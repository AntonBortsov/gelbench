use super::SimState;
use crate::transport;

// One RK4 departure-point step per macrostep for the spatially varying field.
const CISL_NSUB: usize = 1;
// Radius (cells) of the Zalesak flux limiter for fractional displacements.
const ZALESAK_WINDOW: usize = 1;
impl SimState {
    /// Commit outstanding row displacements before an edit that changes velocity or concentration.
    pub fn apply_pending_advection(&mut self) {
        if self.shift_cells.iter().all(|&s| s == 0.0) {
            return;
        }
        for row in 0..self.shift_cells.len() {
            self.commit_row_drift(row);
        }
        self.view_valid = false;
    }

    pub(super) fn commit_row_drift(&mut self, row: usize) {
        let shift = std::mem::take(&mut self.shift_cells[row]);
        if shift == 0.0 {
            return;
        }
        let n = self.params.topology.cells_per_lane.max(1);
        let range = row * n..(row + 1) * n;
        let outflow = Self::commit_displacement(
            &mut self.conc[range.clone()],
            &mut self.adv[range],
            &mut self.faces,
            shift,
            self.dx,
            self.params.open_anode,
        );
        self.exited[row / self.params.topology.subcolumns] += outflow;
    }

    fn commit_displacement(
        c: &mut [f64],
        scratch: &mut [f64],
        faces: &mut [f64],
        shift: f64,
        dx: f64,
        open: bool,
    ) -> f64 {
        let (fraction, mut exited) = Self::rotate_integer_shift(c, shift, open);
        if fraction != 0.0 {
            exited += Self::advect_segment(c, scratch, faces, fraction * dx, 1.0, dx, open);
            c.copy_from_slice(scratch);
        }
        exited
    }

    /// Returns fractional displacement and outflow. The cathode retains mass.
    fn rotate_integer_shift(seg: &mut [f64], shift_cells: f64, open_right: bool) -> (f64, f64) {
        let n = seg.len();
        let k = shift_cells.trunc();
        let frac = shift_cells - k;
        let k = k as isize;
        if k == 0 {
            return (frac, 0.0);
        }
        let ku = (k.unsigned_abs()).min(n);
        let mut exited = 0.0;
        if k > 0 {
            let lost: f64 = seg[n - ku..].iter().sum();
            seg.copy_within(0..n - ku, ku);
            seg[..ku].fill(0.0);
            if open_right {
                exited = lost;
            } else {
                seg[n - 1] += lost;
            }
        } else {
            let lost: f64 = seg[..ku].iter().sum();
            seg.copy_within(ku..n, 0);
            seg[n - ku..].fill(0.0);
            seg[0] += lost;
        }
        (frac, exited)
    }

    /// Materialize outstanding fractional shifts in the observation buffer. Reading the result leaves the transport state unchanged.
    pub(super) fn refresh_view(&mut self) {
        if self.view_valid {
            return;
        }
        self.view_conc.copy_from_slice(&self.conc);
        self.view_exited.fill(0.0);
        if self.shift_cells.iter().any(|&s| s != 0.0) {
            let topo = self.params.topology;
            let open = self.params.open_anode;
            let n = topo.cells_per_lane.max(1);
            let n_pops = self.n_pops_total();
            let k_sub = topo.subcolumns;
            for lane in 0..topo.lanes {
                for pop in 0..n_pops {
                    for k in 0..k_sub {
                        let row = (lane * n_pops + pop) * k_sub + k;
                        let sh = self.shift_cells[row];
                        if sh == 0.0 {
                            continue;
                        }
                        let base = row * n;
                        let gone = Self::advect_segment(
                            &self.view_conc[base..base + n],
                            &mut self.adv[base..base + n],
                            &mut self.faces,
                            sh * self.dx,
                            1.0,
                            self.dx,
                            open,
                        );
                        self.view_exited[row / k_sub] += gone;
                        let (s, e) = (base, base + n);
                        self.view_conc[s..e].copy_from_slice(&self.adv[s..e]);
                    }
                }
            }
        }
        self.view_valid = true;
    }

    /// Subdivide above CFL 0.9 to keep the donor-cell update positive.
    fn advect_segment(
        src: &[f64],
        dst: &mut [f64],
        faces: &mut [f64],
        v: f64,
        dt_adv: f64,
        dx: f64,
        open_right: bool,
    ) -> f64 {
        const CFL_ADVECT_MAX: f64 = 0.9;
        let lam = (v * dt_adv / dx).abs();
        let nsub = ((lam / CFL_ADVECT_MAX).ceil().max(1.0)) as usize;
        if nsub == 1 {
            return Self::advect_segment_once(src, dst, faces, v, dt_adv, dx, open_right);
        }
        let tau = dt_adv / nsub as f64;
        let mut tmp = vec![0.0f64; src.len()];
        let mut exited = Self::advect_segment_once(src, dst, faces, v, tau, dx, open_right);
        for i in 1..nsub {
            if i % 2 == 1 {
                exited += Self::advect_segment_once(dst, &mut tmp, faces, v, tau, dx, open_right);
            } else {
                exited += Self::advect_segment_once(&tmp, dst, faces, v, tau, dx, open_right);
            }
        }
        if nsub.is_multiple_of(2) {
            dst.copy_from_slice(&tmp);
        }
        exited
    }

    /// Take one limited flux step and correct its centroid to the analytic constant-velocity displacement when no mass has left.
    fn advect_segment_once(
        src: &[f64],
        dst: &mut [f64],
        faces: &mut [f64],
        v: f64,
        dt_adv: f64,
        dx: f64,
        open_right: bool,
    ) -> f64 {
        let m: f64 = src.iter().sum();
        let com_before = if m > 0.0 {
            src.iter()
                .enumerate()
                .map(|(i, c)| i as f64 * c)
                .sum::<f64>()
                / m
        } else {
            0.0
        };
        let exited =
            transport::advect_fct_open(src, dst, faces, v, dt_adv, dx, ZALESAK_WINDOW, open_right);
        let m_new: f64 = dst.iter().sum();
        if m_new > 0.0 && m > 0.0 && exited == 0.0 {
            let com_after = dst
                .iter()
                .enumerate()
                .map(|(i, c)| i as f64 * c)
                .sum::<f64>()
                / m_new;
            let target = com_before + v * dt_adv / dx;
            transport::shift_profile_cells(dst, target - com_after);
        }
        exited
    }

    /// Uniform fields allow deferred row shifts. Spatially varying fields require fresh characteristics each step.
    pub(super) fn advect_populations(&mut self, dt: f64) {
        let topo = self.params.topology;
        let n = topo.cells_per_lane.max(1);
        let populations = self.n_pops_total();
        let open = self.params.open_anode;
        for lane in 0..topo.lanes {
            let lane_factor = self.casting_factor(lane) * self.lane_thermal_factor(lane);
            let variable_field = super::electrolyte::field_for_lane(
                self.medium.as_ref(),
                self.field_phi.as_deref(),
                lane,
            )
            .filter(|p| super::electrolyte::field_is_nonuniform(p));
            for pop in 0..populations {
                for k in 0..topo.subcolumns {
                    let row = (lane * populations + pop) * topo.subcolumns + k;
                    let mut velocity = self.velocities[pop] * self.subcolumn_theta[k] * lane_factor;
                    if let Some(theta) = &self.dynamic_subcolumn_theta {
                        velocity *= theta[k];
                    }
                    let shift = &mut self.shift_cells[row];
                    let range = row * n..(row + 1) * n;
                    match variable_field {
                        None => {
                            *shift += velocity * dt / self.dx;
                            if shift.abs() >= 1.0 {
                                let (fraction, exited) =
                                    Self::rotate_integer_shift(&mut self.conc[range], *shift, open);
                                *shift = fraction;
                                self.exited[row / topo.subcolumns] += exited;
                            }
                        }
                        Some(phi) => {
                            let c = &mut self.conc[range.clone()];
                            let mut outflow = 0.0;
                            if *shift != 0.0 {
                                outflow = Self::commit_displacement(
                                    c,
                                    &mut self.adv[range],
                                    &mut self.faces,
                                    std::mem::take(shift),
                                    self.dx,
                                    open,
                                );
                            }
                            let local_velocity =
                                |x| velocity * transport::interp_on(&self.xc, phi, x);
                            let (next, exited) = transport::cisl_step_open(
                                c,
                                &self.xe,
                                &local_velocity,
                                dt,
                                CISL_NSUB,
                                open,
                            );
                            c.copy_from_slice(&next);
                            self.exited[row / topo.subcolumns] += outflow + exited;
                        }
                    }
                }
            }
        }
    }
}
