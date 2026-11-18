//! Migration offsets from the transverse temperature profile and a 5 kb reference band.

use super::heat_transfer::bisect_result;
use super::{SteadyProfileError, ThermalConfig};

impl ThermalConfig {
    /// Расстояние от бокового края геля до крайней дорожки, м.
    pub const LANE_INSET_M: f64 = 0.010;

    /// Разность пути центральной и крайней дорожек, мм. Опорный путь при 5 В/см масштабируется по полю и перегреву.
    pub fn spread_mm_inset(
        &self,
        path_5cm_cm: f64,
        lane_inset_m: f64,
    ) -> Result<f64, SteadyProfileError> {
        let kt = crate::mobility::KT_PER_CELSIUS;
        let (t, _h) = self.steady_profile()?;
        let n = t.len();
        let dy = self.width_m / (n - 1) as f64;
        let idx = |y_m: f64| -> f64 {
            let pos = (y_m + self.width_m / 2.0) / dy;
            let i = pos.floor() as usize;
            if i + 1 >= n {
                return t[n - 1];
            }
            let frac = pos - i as f64;
            t[i] * (1.0 - frac) + t[i + 1] * frac
        };
        let t_center = idx(0.0);
        let t_edge = idx(self.width_m / 2.0 - lane_inset_m);
        let path_cm = path_5cm_cm * (self.e_v_per_cm / 5.0);
        Ok(path_cm * 10.0 * kt * (t_center - t_edge))
    }

    /// Перепад пути между центром геля и дорожкой в 1 см от бокового края.
    pub fn spread_mm(&self, path_5cm_cm: f64) -> Result<f64, SteadyProfileError> {
        self.spread_mm_inset(path_5cm_cm, Self::LANE_INSET_M)
    }

    /// Опорный путь за 45 минут при заданном поле, см. Порог изгиба применим, пока полоса остаётся в геле.
    pub fn band_path_cm(&self, e_v_per_cm: f64, path_5cm_cm: f64) -> f64 {
        path_5cm_cm * (e_v_per_cm / 5.0)
    }

    /// Перепад температуры между центром и дорожкой в 1 см от края.
    pub fn drop_center_to_lane(&self, t: &[f64]) -> f64 {
        let n = t.len();
        let dy = self.width_m / (n - 1) as f64;
        let lane = n - 1 - (Self::LANE_INSET_M / dy).round() as usize;
        t[n / 2] - t[lane]
    }

    /// Поле, при котором разброс достигает 1 мм. `None`, если корень вне диапазона 0.3–100 В/см.
    pub fn smile_threshold_e(&self, path_5cm_cm: f64) -> Result<Option<f64>, SteadyProfileError> {
        self.smile_threshold_e_inset(path_5cm_cm, Self::LANE_INSET_M)
    }

    /// Порог изгиба при сравнении центра с самим торцом геля.
    pub fn smile_threshold_e_metric_b(
        &self,
        path_5cm_cm: f64,
    ) -> Result<Option<f64>, SteadyProfileError> {
        self.smile_threshold_e_inset(path_5cm_cm, 0.0)
    }

    fn smile_threshold_e_inset(
        &self,
        path_5cm_cm: f64,
        lane_inset_m: f64,
    ) -> Result<Option<f64>, SteadyProfileError> {
        let mut eval = |e| {
            let cfg = ThermalConfig {
                e_v_per_cm: e,
                ..*self
            };
            Ok(cfg.spread_mm_inset(path_5cm_cm, lane_inset_m)? - 1.0)
        };
        let mut lo = 0.3;
        let mut flo = eval(lo)?;
        if flo == 0.0 {
            return Ok(Some(lo));
        }
        while lo < 100.0 {
            let hi = (lo * 1.15_f64).max(lo + 0.1).min(100.0);
            let fhi = eval(hi)?;
            if fhi == 0.0 {
                return Ok(Some(hi));
            }
            if flo.signum() != fhi.signum() {
                return bisect_result(&mut eval, lo, hi);
            }
            lo = hi;
            flo = fhi;
        }
        Ok(None)
    }
}

/// Опорный путь модели за 45 минут при 5 В/см для 5 кб в 1% агарозе (TAE), см.
/// Масштабирование по полю сохраняет эту опорную концентрацию геля.
pub fn reference_band_path_cm() -> f64 {
    const REF_BP: f64 = 5000.0;
    const REF_PERCENT: f64 = 1.0;
    const REF_E_V_CM: f64 = 5.0;
    const RUN_S: f64 = 45.0 * 60.0;
    crate::mobility2::MobilityV2::mobility_cm2_v_s(
        REF_BP,
        REF_PERCENT,
        REF_E_V_CM,
        crate::mobility::T_REF_K,
    ) * REF_E_V_CM
        * RUN_S
}

/// Температурный множитель подвижности в центре каждой из K поперечных полос.
pub fn theta_from_delta_t(cfg: &ThermalConfig, d_t: &[f64], k: usize) -> Vec<f64> {
    let kt = crate::mobility::KT_PER_CELSIUS;
    (0..k)
        .map(|kk| {
            let y = -cfg.width_m / 2.0 + (kk as f64 + 0.5) * (cfg.width_m / k as f64);
            1.0 + kt * local_delta_t_at(cfg, d_t, y)
        })
        .collect()
}

/// Линейная интерполяция ΔT на расстоянии y (м) от центра. За пределами геля берётся крайнее значение.
pub fn local_delta_t_at(cfg: &ThermalConfig, d_t: &[f64], y_m: f64) -> f64 {
    let dy = cfg.width_m / (d_t.len() - 1) as f64;
    let pos = ((y_m + cfg.width_m / 2.0) / dy).clamp(0.0, (d_t.len() - 1) as f64);
    let i = (pos.floor() as usize).min(d_t.len() - 2);
    let frac = pos - i as f64;
    d_t[i] * (1.0 - frac) + d_t[i + 1] * frac
}

/// Стационарные тепловые множители для `SimState::set_subcolumn_theta`.
pub fn theta_profile(cfg: &ThermalConfig, k: usize) -> Result<Vec<f64>, SteadyProfileError> {
    let (t, _h) = cfg.steady_profile()?;
    Ok(theta_from_delta_t(cfg, &t, k))
}
