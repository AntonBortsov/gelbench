//! Heat transfer at the upper and lower gel faces.

use super::G;
use crate::fmath::DetMath;

/// Среда для корреляции Нуссельта–Рэлея.
#[derive(Debug, Clone, Copy)]
pub struct Medium {
    /// β = (1/T) для газов; табличное для жидкостей, 1/К.
    pub beta: f64,
    /// Кинематическая вязкость, м²/с.
    pub nu: f64,
    /// Температуропроводность, м²/с.
    pub alpha: f64,
    /// Теплопроводность, Вт/(м·К).
    pub k: f64,
}

/// Воздух, 293 K (литература: CRC Handbook / Incropera гл. 9).
pub const AIR: Medium = Medium {
    beta: 1.0 / 293.0,
    nu: 1.57e-5,
    alpha: 2.22e-5,
    k: 0.026,
};

/// Вода, 293 K (литература: CRC Handbook / Incropera гл. 9).
pub const WATER: Medium = Medium {
    beta: 2.07e-4,
    nu: 1.004e-6,
    alpha: 1.43e-7,
    k: 0.598,
};

/// Корреляция McAdams для свободной конвекции, Nu = C·Ra^(1/4) (Incropera & DeWitt, гл. 9, точность ±30%).
#[inline]
pub fn h_natural_convection(med: &Medium, delta_t: f64, l_c: f64, c_coef: f64) -> f64 {
    let ra = G * med.beta * delta_t.max(1e-6) * l_c.powi(3) / (med.nu * med.alpha);
    c_coef * ra.dpowf(0.25) * med.k / l_c
}

/// Число Рэлея для замкнутого горизонтального слоя толщиной `depth_m`.
#[inline]
pub fn rayleigh_horizontal_layer(med: &Medium, delta_t: f64, depth_m: f64) -> f64 {
    assert!(depth_m.is_finite() && depth_m > 0.0);
    G * med.beta * delta_t.max(0.0) * depth_m.powi(3) / (med.nu * med.alpha)
}

/// Корреляция Hollands–Unny–Raithby–Konicek для замкнутого слоя, нагретого снизу.
/// При Ra ≤ 1708 тепло переносится только проводимостью (Nu = 1).
#[inline]
pub fn nusselt_horizontal_layer(ra: f64) -> f64 {
    if !ra.is_finite() || ra <= 0.0 {
        return 1.0;
    }
    1.0 + 1.44 * (1.0 - 1708.0 / ra).max(0.0) + ((ra / 5830.0).dpowf(1.0 / 3.0) - 1.0).max(0.0)
}

/// Поддон и внешняя теплоотдача действуют последовательно. Для перемешанной ванны внешний h бесконечен.
#[inline]
pub(super) fn tray_series_heat_transfer(
    external_h_w_m2_k: f64,
    tray_thickness_m: f64,
    tray_k_w_m_k: f64,
) -> f64 {
    assert!(external_h_w_m2_k > 0.0);
    assert!(tray_thickness_m.is_finite() && tray_thickness_m >= 0.0);
    assert!(tray_k_w_m_k > 0.0);
    if tray_thickness_m == 0.0 || tray_k_w_m_k == f64::INFINITY {
        return external_h_w_m2_k;
    }
    let tray_h = tray_k_w_m_k / tray_thickness_m;
    if external_h_w_m2_k == f64::INFINITY {
        return tray_h;
    }
    1.0 / (1.0 / tray_h + 1.0 / external_h_w_m2_k)
}

/// Бисекция на отрезке. Возвращает `None`, если концы не заключают корень.
pub fn bisect(mut f: impl FnMut(f64) -> f64, lo: f64, hi: f64) -> Option<f64> {
    let mut lo = lo;
    let mut hi = hi;
    let mut flo = f(lo);
    let fhi = f(hi);
    if flo == 0.0 {
        return Some(lo);
    }
    if fhi == 0.0 {
        return Some(hi);
    }
    if flo.is_nan() || fhi.is_nan() || flo.signum() == fhi.signum() {
        return None;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            return Some(mid); // достигнут предел точности f64
        }
        let fm = f(mid);
        if fm == 0.0 {
            return Some(mid);
        }
        if fm.signum() == flo.signum() {
            lo = mid;
            flo = fm;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

/// Различает отсутствие корня в интервале и ошибку вычисления теплового профиля.
pub(super) fn bisect_result<E>(
    mut f: impl FnMut(f64) -> Result<f64, E>,
    lo: f64,
    hi: f64,
) -> Result<Option<f64>, E> {
    let mut lo = lo;
    let mut hi = hi;
    let mut flo = f(lo)?;
    let fhi = f(hi)?;
    if flo == 0.0 {
        return Ok(Some(lo));
    }
    if fhi == 0.0 {
        return Ok(Some(hi));
    }
    if flo.is_nan() || fhi.is_nan() || flo.signum() == fhi.signum() {
        return Ok(None);
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            return Ok(Some(mid));
        }
        let fm = f(mid)?;
        if fm == 0.0 {
            return Ok(Some(mid));
        }
        if fm.signum() == flo.signum() {
            lo = mid;
            flo = fm;
        } else {
            hi = mid;
        }
    }
    Ok(Some(0.5 * (lo + hi)))
}
