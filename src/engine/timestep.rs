//! Fixed-order operator splitting. Electrolyte properties and mobility use the
//! temperature at the start of the macrostep. Electrical heating affects the next one.

use super::{SimFault, SimState};
use crate::{stain, transport};

impl SimState {
    /// Advance one configured macrostep. On medium failure the model clock stops.
    pub fn step(&mut self) {
        if self.fault.is_some() {
            return;
        }
        self.view_valid = false;
        let dt = f64::from(self.params.macro_dt_s);
        let bath_temp_k = self.update_thermal_feedback();
        let power = self.advance_circuit(dt, bath_temp_k);
        if let Err(fault) = self.advance_medium(dt) {
            self.fault = Some(fault);
            return;
        }
        self.advance_stain(dt);
        self.diffuse_populations(dt);
        self.advect_populations(dt);
        if let (Some(bath), Some((total_w, gel_w))) = (&mut self.bath_thermal, power) {
            bath.advance(dt, total_w, gel_w);
        }
        self.t += dt;
        self.steps += 1;
    }

    fn update_thermal_feedback(&mut self) -> f64 {
        let (feedback_temp_k, medium_feedback_temp_k) = match &self.bath_thermal {
            Some(t) => (
                t.feedback_temperature_k(),
                t.feedback_mean_gel_temperature_k(),
            ),
            None => (self.params.gel.temp_k, self.params.gel.temp_k),
        };
        if let Some(t) = &self.bath_thermal {
            self.dynamic_subcolumn_theta = Some(t.theta_for_subcolumns());
            if let Some(y) = &self.lane_y_m {
                self.dynamic_lane_theta = Some(t.theta_for_positions(y));
            }
        }
        if let Some(m) = self.medium.as_mut() {
            m.set_temperature_k(medium_feedback_temp_k);
        }
        feedback_temp_k
    }

    /// Use the same gel resistivity for the ammeter, gel power and population drift.
    fn advance_circuit(&mut self, dt: f64, bath_temp_k: f64) -> Option<(f64, f64)> {
        let chamber = self.chamber.as_mut()?;
        let rho = self.medium.as_ref().map(|m| m.gel_resistivity_ohm_m());
        let op = chamber.step_with_gel(dt, bath_temp_k, rho);
        let gel_power = chamber.gel_power_w_at(&op, bath_temp_k, rho);
        self.apply_field_kinetics_at(op.e_gel_v_per_cm, bath_temp_k);
        if let Some(medium) = &mut self.medium {
            medium.set_u_gel(op.e_gel_v_per_cm * f64::from(self.params.gel_length_cm));
        }
        Some((op.p_w.max(0.0), gel_power))
    }

    fn advance_medium(&mut self, dt: f64) -> Result<(), SimFault> {
        if let Some(medium) = &mut self.medium {
            medium.step(dt).map_err(SimFault::Medium)?;
        }
        Ok(())
    }

    fn advance_stain(&mut self, dt: f64) {
        let Some(staining) = &mut self.stain else {
            return;
        };
        if !matches!(staining.mode, stain::StainMode::InGel { .. }) {
            return;
        }
        // Free ethidium is cationic. An open circuit stops its drift as well as DNA drift.
        let n = self.params.topology.cells_per_lane.max(1);
        let mut face_velocities = vec![0.0; n + 1];
        for (lane, concentration) in staining.free_m.iter_mut().enumerate() {
            let phi = super::electrolyte::field_for_lane(
                self.medium.as_ref(),
                self.field_phi.as_deref(),
                lane,
            );
            for (face, velocity) in face_velocities.iter_mut().enumerate() {
                let scale = match phi {
                    None => 1.0,
                    Some(p) if face == 0 => p[0],
                    Some(p) if face == n => p[n - 1],
                    Some(p) => 0.5 * (p[face - 1] + p[face]),
                };
                *velocity = -stain::MU_ETBR_CM2_V_S * self.kinetics_e_v_cm * scale;
            }
            *concentration = transport::sg_step_reservoir(
                concentration,
                self.dx,
                &face_velocities,
                stain::D_ETBR_CM2_S,
                dt,
                0.0,
                0.0,
            );
        }
    }

    fn diffuse_populations(&mut self, dt_macro: f64) {
        let d_max = self.diffusivities.iter().cloned().fold(0.0f64, f64::max);

        let (n_dif, dt_dif) = transport::substep_count(dt_macro, self.dx, 0.0, d_max);
        if n_dif > 0 && d_max > 0.0 {
            let topo = self.params.topology;
            let n = topo.cells_per_lane.max(1);
            let n_pops = self.n_pops_total();
            let k_sub = topo.subcolumns;
            for _ in 0..n_dif {
                for lane in 0..topo.lanes {
                    for pop in 0..n_pops {
                        for k in 0..k_sub {
                            let base = (lane * n_pops + pop) * k_sub * n + k * n;
                            transport::diffuse_explicit(
                                &self.conc[base..base + n],
                                &mut self.adv[base..base + n],
                                self.diffusivities[pop],
                                dt_dif,
                                self.dx,
                            );
                        }
                    }
                }
                self.conc.copy_from_slice(&self.adv);
            }
        }
    }
}
