use super::SimState;
use crate::{fmath::DetMath, loading};

impl SimState {
    /// Remove all analytes from a lane after accounting for their pending displacement.
    pub fn clear_lane(&mut self, lane: usize) {
        self.apply_pending_advection();
        let topo = self.params.topology;
        let n = topo.cells_per_lane.max(1);
        let n_pops = self.n_pops_total();
        for pop in 0..n_pops {
            for k in 0..topo.subcolumns {
                let base = (lane * n_pops + pop) * topo.subcolumns * n + k * n;
                self.conc[base..base + n].fill(0.0);
            }
        }
        self.view_valid = false;
    }

    /// Distribute one band's mass over a raised cosine, normalized across all transverse subcolumns.
    fn seed_band(&mut self, lane: usize, pop: usize, mass_ng: f64, half_w_cells: f64, add: bool) {
        let topo = self.params.topology;
        let n = topo.cells_per_lane.max(1);
        let n_pops = self.n_pops_total();
        let k_sub = topo.subcolumns.max(1);
        let peak = ((self.params.well_position_cm as f64 / self.dx).floor() as usize)
            .clamp(2, n.saturating_sub(4));
        let half_w = half_w_cells.max(1.5);
        let mut profile = vec![0.0f64; n];
        let mut psum = 0.0;
        for (x, pv) in profile.iter_mut().enumerate() {
            *pv = Self::well_profile_height(x, peak, half_w);
            psum += *pv;
        }
        assert!(psum > 0.0);
        let amp = mass_ng / (psum * k_sub as f64);
        for k in 0..k_sub {
            let row = (lane * n_pops + pop) * k_sub + k;
            self.commit_row_drift(row);
            let base = row * n;
            for x in 0..n {
                let v = amp * profile[x];
                if add {
                    self.conc[base + x] += v;
                } else {
                    self.conc[base + x] = v;
                }
            }
        }
        self.view_valid = false;
    }

    /// Seed a lane and account for material retained, lost to flotation or spilled into neighboring wells.
    pub fn load_sample(
        &mut self,
        lane: usize,
        comb: &loading::Comb,
        load: &loading::SampleLoad,
    ) -> loading::LoadReport {
        let n_dna = self.params.population_sizes_bp.len();
        assert_eq!(
            load.masses_ng.len(),
            n_dna,
            "masses_ng обязан быть параллелен population_sizes_bp"
        );
        let topo = self.params.topology;
        assert!(lane < topo.lanes, "дорожки {lane} нет");
        let cap = comb.well_capacity_ul();
        let has_left = lane > 0;
        let has_right = lane + 1 < topo.lanes;
        let (stay, per_side, lost_frac) =
            loading::overflow_split(load.volume_ul, cap, has_left, has_right);
        let f_loss = loading::density_loss_fraction(load);
        let m_total: f64 = load.masses_ng.iter().sum();

        let half_w_cells = 0.5 * comb.tooth_thickness_m * 1.0e2 / self.dx;
        for (pop, &m) in load.masses_ng.iter().enumerate() {
            self.seed_band(lane, pop, m * stay * (1.0 - f_loss), half_w_cells, false);
            if per_side > 0.0 {
                if has_left {
                    self.seed_band(lane - 1, pop, m * per_side, half_w_cells * 1.5, true);
                }
                if has_right {
                    self.seed_band(lane + 1, pop, m * per_side, half_w_cells * 1.5, true);
                }
            }
        }
        let n_sides = has_left as u32 + has_right as u32;
        loading::LoadReport {
            loaded_ng: m_total * stay * (1.0 - f_loss),
            density_loss_ng: m_total * stay * f_loss,
            spilled_to_neighbors_ng: m_total * per_side * n_sides as f64,
            spilled_lost_ng: m_total * lost_frac,
            overload: loading::overload_status(m_total * stay * (1.0 - f_loss), comb.tooth_width_m),
            overflow_fraction: (load.volume_ul - cap).max(0.0) / load.volume_ul,
        }
    }

    pub fn load_dyes(&mut self, lane: usize, comb: &loading::Comb, masses_ng: &[f64]) {
        let n_dna = self.params.population_sizes_bp.len();
        assert_eq!(
            masses_ng.len(),
            self.params.dyes.len(),
            "masses_ng обязан быть параллелен params.dyes"
        );
        assert!(lane < self.params.topology.lanes, "дорожки {lane} нет");
        let half_w_cells = 0.5 * comb.tooth_thickness_m * 1.0e2 / self.dx;
        for (i, &m) in masses_ng.iter().enumerate() {
            self.seed_band(lane, n_dna + i, m, half_w_cells, false);
        }
    }

    pub fn dye_population_index(&self, i: usize) -> usize {
        assert!(i < self.params.dyes.len(), "красителя {i} нет");
        self.params.population_sizes_bp.len() + i
    }

    pub(super) fn well_profile_height(cell: usize, center: usize, half_width_cells: f64) -> f64 {
        let distance = ((cell as f64 - center as f64).abs() / half_width_cells).min(1.0);
        0.5 * (1.0 + (std::f64::consts::PI * distance).dcos())
    }
}
