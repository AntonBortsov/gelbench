use super::SimState;
use crate::{fmath::DetMath, rng};

/// A single column resolves differential displacement below 0.1 band widths. Otherwise target 0.25 between adjacent columns.
/// The profile already spans the physical width supplied by the caller.
pub fn select_subcolumn_count(
    theta: &[f64],
    _width_m: f64,
    band_speed_cm_s: f64,
    run_s: f64,
    sigma_band_cm: f64,
) -> usize {
    let band_speed_cm_s = band_speed_cm_s.abs();
    let k_max = theta.len();
    assert!(k_max >= 1);
    let th_c = theta[theta.len() / 2];
    let max_disp_diff =
        band_speed_cm_s * run_s * theta.iter().map(|t| (t - th_c).abs()).fold(0.0, f64::max);
    if max_disp_diff < 0.1 * sigma_band_cm {
        return 1;
    }
    let sample = |kk: usize, k: usize| -> f64 {
        let frac = (kk as f64 + 0.5) / k as f64;
        let pos = frac * (theta.len() - 1) as f64;
        let i = pos.floor() as usize;
        let j = (i + 1).min(theta.len() - 1);
        let f = pos - i as f64;
        theta[i] * (1.0 - f) + theta[j] * f
    };
    for k in 2..=k_max {
        let mut worst = 0.0f64;
        for kk in 0..k - 1 {
            let d = (sample(kk + 1, k) - sample(kk, k)).abs();
            worst = worst.max(d);
        }
        if band_speed_cm_s * run_s * worst <= 0.25 * sigma_band_cm {
            return k;
        }
    }
    k_max
}
impl SimState {
    /// Set the mobility multiplier for each transverse subcolumn. Materialize pending displacement before switching the profile.
    pub fn set_subcolumn_theta(&mut self, theta: &[f64]) {
        assert!(
            theta.len() == self.params.topology.subcolumns,
            "theta.len() = {} против K = {}",
            theta.len(),
            self.params.topology.subcolumns
        );
        self.apply_pending_advection();
        self.subcolumn_theta.clear();
        self.subcolumn_theta.extend_from_slice(theta);
        self.view_valid = false;
    }

    /// Apply a seeded two-harmonic variation in lane mobility. This is a qualitative casting artifact. Zero restores a flat gel. Waviness and an active medium loop are mutually exclusive.
    pub fn set_casting_waviness(&mut self, amplitude: f64, seed: u64) {
        const WAVELENGTH_LANES: f64 = 6.18;
        const SECOND_HARMONIC: f64 = 0.35;
        assert!(
            (0.0..0.5).contains(&amplitude),
            "амплитуда волнистости {amplitude} вне [0; 0.5)"
        );
        assert!(
            self.medium.is_none(),
            "№25 совместно с петлёй среды v1 не поддерживает (PHYSICS §32)"
        );
        self.apply_pending_advection();
        if amplitude == 0.0 {
            self.casting_theta_lane = None;
            self.view_valid = false;
            return;
        }
        let lanes = self.params.topology.lanes;
        let mut r = rng::Rng::new(seed);
        let tau = core::f64::consts::TAU;
        let ph1 = r.next_f64() * tau;
        let ph2 = r.next_f64() * tau;
        let w = tau / WAVELENGTH_LANES;
        let norm = 1.0 + SECOND_HARMONIC;
        self.casting_theta_lane = Some(
            (0..lanes)
                .map(|l| {
                    let x = l as f64;
                    let s = ((w * x + ph1).dsin() + SECOND_HARMONIC * (2.0 * w * x + ph2).dsin())
                        / norm;
                    1.0 + amplitude * s
                })
                .collect(),
        );
        self.view_valid = false;
    }

    pub fn casting_theta(&self) -> Option<&[f64]> {
        self.casting_theta_lane.as_deref()
    }

    pub(super) fn casting_factor(&self, lane: usize) -> f64 {
        match &self.casting_theta_lane {
            Some(c) => c[lane],
            None => 1.0,
        }
    }

    /// Set transverse lane centers in meters from the gel center. An empty slice disables lane-wise thermal feedback.
    pub fn set_lane_positions_m(&mut self, y_m: &[f64]) {
        assert!(
            y_m.is_empty() || y_m.len() == self.params.topology.lanes,
            "положений {} против дорожек {}",
            y_m.len(),
            self.params.topology.lanes
        );
        assert!(
            y_m.iter().all(|y| y.is_finite()),
            "положения дорожек обязаны быть конечными"
        );
        self.apply_pending_advection();
        if y_m.is_empty() {
            self.lane_y_m = None;
            self.dynamic_lane_theta = None;
        } else {
            self.lane_y_m = Some(y_m.to_vec());
            self.dynamic_lane_theta = self
                .bath_thermal
                .as_ref()
                .map(|t| t.theta_for_positions(y_m));
        }
        self.view_valid = false;
    }

    pub fn lane_positions_m(&self) -> Option<&[f64]> {
        self.lane_y_m.as_deref()
    }

    /// Lane mobility multipliers from the current bath-temperature profile.
    pub fn lane_theta(&self) -> Option<&[f64]> {
        self.dynamic_lane_theta.as_deref()
    }

    pub(super) fn lane_thermal_factor(&self, lane: usize) -> f64 {
        match &self.dynamic_lane_theta {
            Some(t) => t[lane],
            None => 1.0,
        }
    }
}
