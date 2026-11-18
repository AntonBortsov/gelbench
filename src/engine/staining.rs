use super::SimState;
use crate::{dyes, mobility2, stain};

impl SimState {
    /// Indicator color of bromophenol blue at its mass centroid. Returns no color if the dye is absent or its band is empty.
    pub fn bpb_color(&mut self, lane: usize) -> Option<dyes::BpbColor> {
        let di = self
            .params
            .dyes
            .iter()
            .position(|d| *d == dyes::DyeKind::BromophenolBlue)?;
        let pop = self.params.population_sizes_bp.len() + di;
        let topo = self.params.topology;
        let n = topo.cells_per_lane.max(1);
        let n_pops = self.n_pops_total();
        let k_sub = topo.subcolumns;
        let view = self.concentrations();
        let mut m0 = 0.0;
        let mut mi = 0.0;
        for k in 0..k_sub {
            let base = (lane * n_pops + pop) * k_sub * n + k * n;
            for x in 0..n {
                m0 += view[base + x];
                mi += view[base + x] * x as f64;
            }
        }
        if m0 <= 0.0 {
            return None;
        }
        let cell = (mi / m0).round().clamp(0.0, (n - 1) as f64) as usize;
        let ph = match &self.medium {
            Some(m) => m.ph_profile(lane)[cell],
            None => 8.3,
        };
        Some(dyes::bpb_color_at_ph(ph))
    }

    /// Set the stain mode after manual kinetics overrides. In-gel ethidium bromide changes DNA mobility during transport.
    pub fn enable_staining(&mut self, mode: stain::StainMode) {
        let topo = self.params.topology;
        let n = topo.cells_per_lane.max(1);
        self.stain = Some(stain::StainState::new(mode, topo.lanes, n));
        if let stain::StainMode::InGel { c0_ug_ml } = mode {
            let free_etbr_m = stain::ug_ml_to_molar(c0_ug_ml);
            let e_now = self.kinetics_e_v_cm;
            for p in 0..self.params.population_sizes_bp.len() {
                if self.manual_kinetics[p] {
                    let f = mobility2::MobilityV2::etbr_geometry_factor(
                        self.params.population_sizes_bp[p],
                        self.params
                            .population_forms
                            .get(p)
                            .copied()
                            .unwrap_or_default(),
                        self.params.gel.percent,
                        e_now.abs(),
                        free_etbr_m,
                    );
                    self.velocities[p] *= f;
                    self.manual_v_ref[p] *= f;
                } else {
                    let (v, d) = Self::kinetics_for(
                        &self.params,
                        Some(mode),
                        p,
                        e_now,
                        self.kinetics_temp_k,
                    );
                    self.velocities[p] = v;
                    self.diffusivities[p] = d;
                }
            }
        }
        self.view_valid = false;
    }

    /// Free ethidium bromide concentration along a lane (mol/L).
    pub fn stain_profile(&self, lane: usize) -> Option<&[f64]> {
        self.stain.as_ref().map(|s| s.free_m[lane].as_slice())
    }

    /// DNA fluorescence from bound stain, summed across populations and subcolumns.
    pub fn stain_band_intensity(&mut self, lane: usize) -> Option<Vec<f64>> {
        let topo = self.params.topology;
        let n = topo.cells_per_lane.max(1);
        let n_dna = self.params.population_sizes_bp.len();
        let n_pops = self.n_pops_total();
        let k_sub = topo.subcolumns;
        self.stain.as_ref()?;
        self.refresh_view();
        let st = self.stain.as_ref()?;
        let view = &self.view_conc;
        let mut out = vec![0.0; n];
        for (x, o) in out.iter_mut().enumerate() {
            let nu = stain::nu_bound(st.free_m[lane][x]);
            let mut c_sum = 0.0;
            for pop in 0..n_dna {
                for k in 0..k_sub {
                    let base = (lane * n_pops + pop) * k_sub * n + k * n;
                    c_sum += view[base + x];
                }
            }
            *o = c_sum * nu * stain::FLUOR_ENHANCEMENT;
        }
        Some(out)
    }

    /// Free stain concentration relative to the initial bath concentration, used for the fluorescence background.
    pub fn stain_haze(&self, lane: usize) -> Option<Vec<f64>> {
        let st = self.stain.as_ref()?;
        let c0 = stain::ug_ml_to_molar(st.mode.c0_ug_ml()).max(1e-30);
        Some(st.free_m[lane].iter().map(|&e| e / c0).collect())
    }
}
