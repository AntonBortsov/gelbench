use super::SimState;
use crate::chamber;

/// Instrument readings at a macrostep boundary. Units follow the field names.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observation {
    /// Elapsed model time (s).
    pub t: f64,
    pub steps: u64,
    /// Mass currently on the gel, in loading units.
    pub total_mass: f64,
    pub max_conc: f64,
    /// Chamber current (A). Instrument fields are zero when no chamber is attached.
    pub current_a: f64,
    pub e_gel_v_per_cm: f64,
    pub power_w: f64,
    pub buffer_level_mm: f64,
    /// Maximum gas coverage fraction of the two electrodes.
    pub theta_max: f64,
    pub gel_delta_t_est_k: f64,
    pub melt_risk: bool,
    /// Cumulative anode outflow, in the same units as total_mass.
    pub exited_mass: f64,
}
impl SimState {
    /// Mass, centroid (cm) and standard deviation (cm) on the observation grid, summing over subcolumns.
    pub(crate) fn band_moments(&mut self, lane: usize, pop: usize) -> (f64, f64, f64) {
        let topo = self.params.topology;
        assert!(lane < topo.lanes, "дорожки {lane} нет");
        assert!(pop < self.n_pops_total(), "популяции {pop} нет");
        let n = topo.cells_per_lane.max(1);
        let k_sub = topo.subcolumns;
        let n_pops = self.n_pops_total();
        let dx = self.dx;
        self.refresh_view();
        let (mut m0, mut mx, mut mxx) = (0.0f64, 0.0f64, 0.0f64);
        for k in 0..k_sub {
            let base = ((lane * n_pops + pop) * k_sub + k) * n;
            for i in 0..n {
                let c = self.view_conc[base + i];
                let x = (i as f64 + 0.5) * dx;
                m0 += c;
                mx += c * x;
                mxx += c * x * x;
            }
        }
        if m0 <= 0.0 {
            return (0.0, f64::NAN, f64::NAN);
        }
        let com = mx / m0;
        let var = (mxx / m0 - com * com).max(0.0);
        (m0, com, var.sqrt())
    }

    /// Centroid from the cathode in cm, or NaN for an empty band.
    pub fn band_com_cm(&mut self, lane: usize, pop: usize) -> f64 {
        self.band_moments(lane, pop).1
    }

    /// Standard deviation in cm, or NaN for an empty band.
    pub fn band_sigma_cm(&mut self, lane: usize, pop: usize) -> f64 {
        self.band_moments(lane, pop).2
    }

    /// Mass retained on the gel in the units used to load it. An empty band has zero mass.
    pub fn band_mass(&mut self, lane: usize, pop: usize) -> f64 {
        self.band_moments(lane, pop).0
    }

    /// Mass lost from one band through the anode, including the fractional displacement of the current observation.
    pub fn exited_mass(&mut self, lane: usize, pop: usize) -> f64 {
        assert!(lane < self.params.topology.lanes, "дорожки {lane} нет");
        assert!(pop < self.n_pops_total(), "популяции {pop} нет");
        self.refresh_view();
        let i = lane * self.n_pops_total() + pop;
        self.exited[i] + self.view_exited[i]
    }

    pub fn exited_total(&mut self) -> f64 {
        self.refresh_view();
        self.exited.iter().sum::<f64>() + self.view_exited.iter().sum::<f64>()
    }

    /// Materialize concentrations in lane, population, subcolumn, cell order. Reading never changes the next simulated state.
    pub fn concentrations(&mut self) -> &[f64] {
        self.refresh_view();
        &self.view_conc
    }

    pub fn observe(&mut self) -> Observation {
        self.refresh_view();
        let mut total_mass = 0.0f64;
        let mut max_conc = 0.0f64;
        for &v in &self.view_conc {
            total_mass += v;
            if v > max_conc {
                max_conc = v;
            }
        }
        let (current_a, e_gel_v_per_cm, power_w, buffer_level_mm, theta_max, static_dt_est) =
            match &self.chamber {
                Some(ch) => (
                    ch.last.i_a,
                    ch.last.e_gel_v_per_cm,
                    ch.last.p_w,
                    ch.level_m * 1.0e3,
                    ch.theta_anode.max(ch.theta_cathode),
                    ch.melt_delta_t_estimate_k(&ch.last),
                ),
                None => (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            };
        let (dt_est, melt_risk) = match &self.bath_thermal {
            Some(t) => (
                t.local_delta_t_k().iter().copied().fold(0.0, f64::max),
                t.max_gel_temp_k() >= chamber::MELT_RISK_K,
            ),
            None => (
                static_dt_est,
                self.params.gel.temp_k + static_dt_est >= chamber::MELT_RISK_K,
            ),
        };
        Observation {
            t: self.t,
            steps: self.steps,
            total_mass,
            max_conc,
            current_a,
            e_gel_v_per_cm,
            power_w,
            buffer_level_mm,
            theta_max,
            gel_delta_t_est_k: dt_est,
            melt_risk,
            exited_mass: self.exited.iter().sum::<f64>() + self.view_exited.iter().sum::<f64>(),
        }
    }
}
