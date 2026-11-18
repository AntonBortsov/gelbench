//! Свободный EtBr дрейфует к катоду. В приближении разбавленной ДНК связывание в полосах не истощает краситель.
//! Яркость зависит от локального связывания, геометрия мигрирующей ДНК — от начальной концентрации EtBr.
//! Изотерма МакГи–фон Хиппеля учитывает исключение соседних сайтов (J Mol Biol 86:469, 1974).
//!
//! ```text
//! ν/L = K·(1−nν)·[(1−nν)/(1−(n−1)ν)]^(n−1)
//! ```
//!
//! ν — число связанных молекул на пару оснований, L — свободный EtBr (М).
//! Параметры взяты для умеренной соли (n ≈ 2–3, K ≈ 1–2·10⁵ М⁻¹). Зависимость K от ионной силы не моделируется.

use crate::fmath::DetMath;
use serde::{Deserialize, Serialize};

/// Молярная масса EtBr, г/моль.
pub const MW_ETBR_G_MOL: f64 = 394.3;
/// Сайт-эксклюзия МакГи–фон Хиппеля, пар оснований на молекулу.
pub const N_SITE_BP: f64 = 2.5;
/// Константа ассоциации, М⁻¹.
pub const K_ASSOC_M: f64 = 1.5e5;
/// Раскрутка дуплекса на интеркаляцию, градусов.
pub const UNWIND_DEG: f64 = 26.0;
/// Пар оснований на виток B-формы.
pub const BP_PER_TURN: f64 = 10.5;
/// Подвижность одновалентного катиона, см²/(В·с), оценена из D по Нернсту–Эйнштейну.
pub const MU_ETBR_CM2_V_S: f64 = 1.9e-4;
/// Диффузия свободного EtBr, см²/с (Стокс–Эйнштейн, 394 Да).
pub const D_ETBR_CM2_S: f64 = 4.9e-6;
/// Усиление флуоресценции при интеркаляции.
pub const FLUOR_ENHANCEMENT: f64 = 25.0;
pub fn ug_ml_to_molar(ug_ml: f64) -> f64 {
    ug_ml * 1.0e-3 / MW_ETBR_G_MOL
}

fn isotherm_rhs(nu: f64) -> f64 {
    let a = 1.0 - N_SITE_BP * nu;
    if a <= 0.0 {
        return 0.0;
    }
    let b = 1.0 - (N_SITE_BP - 1.0) * nu;
    K_ASSOC_M * a * (a / b).dpowf(N_SITE_BP - 1.0)
}

/// Число связанных молекул на пару оснований при свободном EtBr в моль/л. Физический интервал 0 ≤ ν < 1/n.
pub fn nu_bound(free_m: f64) -> f64 {
    if free_m <= 0.0 {
        return 0.0;
    }
    let mut lo = 0.0_f64;
    let mut hi = 1.0 / N_SITE_BP - 1e-12;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if free_m * isotherm_rhs(mid) > mid {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Свободный EtBr (М), дающий заданное ν (обратная изотерма, явная).
pub fn free_for_nu(nu: f64) -> f64 {
    assert!(nu > 0.0 && nu < 1.0 / N_SITE_BP);
    nu / isotherm_rhs(nu)
}

/// Предельное дополнительное трение при интеркаляции сверх удлинения контура.
/// Феноменологическая подгонка по наклонам для λ/HindIII (Sigmon & Larcom 1996, Table 1).
/// Форма насыщения ниже 0.25 мкг/мл — экстраполяция.
pub const ETBR_FRICTION_MAX: f64 = 0.1971;
/// Полунасыщение по числу связанных молекул на пару оснований.
pub const ETBR_FRICTION_NU_HALF: f64 = 0.00339;

/// Относительная подвижность с поправкой на трение при заданной доле связывания.
pub fn friction_factor_with(nu: f64, r_max: f64, nu_half: f64) -> f64 {
    if nu <= 0.0 {
        return 1.0;
    }
    1.0 - r_max * nu / (nu + nu_half)
}

pub fn friction_factor(nu: f64) -> f64 {
    friction_factor_with(nu, ETBR_FRICTION_MAX, ETBR_FRICTION_NU_HALF)
}

/// Одна интеркаляция удлиняет контур на шаг B-ДНК (0.34 нм), поэтому его относительная длина равна 1 + ν.
pub fn contour_length_factor(nu: f64) -> f64 {
    assert!(nu.is_finite() && (0.0..1.0 / N_SITE_BP).contains(&nu));
    1.0 + nu
}

/// Свободный EtBr (мкг/мл), при котором раскрутка компенсирует начальную отрицательную суперспиральность.
pub fn inversion_free_ug_ml(sigma0: f64) -> f64 {
    assert!(
        sigma0 < 0.0,
        "инверсия — про отрицательную суперспиральность"
    );
    let nu_star = -sigma0 * 360.0 / (UNWIND_DEG * BP_PER_TURN);
    free_for_nu(nu_star) * MW_ETBR_G_MOL * 1.0e3
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum StainMode {
    /// Краситель внесён в гель, резервуары исходно чистые. Уход EtBr к катоду обедняет анодную часть.
    InGel { c0_ug_ml: f64 },
    /// После прогона в чистом геле равномерная ванна красителя определяет только наблюдаемое свечение.
    PostStain { c0_ug_ml: f64 },
}

impl StainMode {
    pub fn c0_ug_ml(&self) -> f64 {
        match self {
            StainMode::InGel { c0_ug_ml } | StainMode::PostStain { c0_ug_ml } => *c0_ug_ml,
        }
    }
}

#[derive(Debug, Clone)]
pub struct StainState {
    pub mode: StainMode,
    /// Свободный EtBr по дорожкам, моль/л. В режиме постокраски поле постоянно.
    pub free_m: Vec<Vec<f64>>,
}

impl StainState {
    pub fn new(mode: StainMode, lanes: usize, n: usize) -> Self {
        let c0 = ug_ml_to_molar(mode.c0_ug_ml());
        Self {
            mode,
            free_m: vec![vec![c0; n]; lanes],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isotherm_roundtrip_and_limits() {
        // Прямая и обратная изотермы согласованы; разбавленный предел
        // ν ≈ K·L при L → 0.
        for &nu in &[0.01, 0.05, 0.111, 0.2, 0.35] {
            let l = free_for_nu(nu);
            let back = nu_bound(l);
            assert!((back - nu).abs() < 1e-9, "ν {nu}: назад {back}");
        }
        let l_small = 1e-9;
        let nu_small = nu_bound(l_small);
        assert!(
            ((nu_small / (K_ASSOC_M * l_small)) - 1.0).abs() < 0.01,
            "разбавленный предел: ν/KL = {}",
            nu_small / (K_ASSOC_M * l_small)
        );
    }

    #[test]
    fn contour_extension_at_half_ug_ml_has_no_fitted_slowdown_constant() {
        // Э6: при 0.5 мкг/мл единственный вход в подвижность — измеримая
        // доля интеркаляции. Скорость из неё предсказывает mobility2.
        let nu = nu_bound(ug_ml_to_molar(0.5));
        assert!(
            (contour_length_factor(nu) - 1.111).abs() < 0.005,
            "L_c(EtBr)/L_c(0) = {} при ν = {nu:.3}",
            contour_length_factor(nu)
        );
    }

    #[test]
    fn anchor_13_inversion_window_is_predicted() {
        // №13 — ПРЕДСКАЗАНИЕ (не подгонка): литературные n, K, 26°/10.5
        // дают окно свободного EtBr 0.24–0.37 мкг/мл для σ₀ −0.05…−0.07.
        let lo = inversion_free_ug_ml(-0.05);
        let hi = inversion_free_ug_ml(-0.07);
        assert!(
            (0.22..=0.27).contains(&lo),
            "σ₀=−0.05: инверсия при {lo:.3} мкг/мл (якорь ~0.24)"
        );
        assert!(
            (0.34..=0.40).contains(&hi),
            "σ₀=−0.07: инверсия при {hi:.3} мкг/мл (якорь ~0.37)"
        );
        // Стандартные 0.5 мкг/мл — за точкой перехода для всего диапазона.
        assert!(0.5 > hi, "0.5 мкг/мл обязан быть за инверсией");
    }

    #[test]
    fn nu_at_standard_gel_concentration() {
        // Санитарная точка: 0.5 мкг/мл = 1.27 мкМ → ν ≈ 0.111 (покрытие
        // 28% от насыщения 1/n) — числа, на которые ссылаются доки.
        let nu = nu_bound(ug_ml_to_molar(0.5));
        assert!((nu - 0.111).abs() < 0.005, "ν = {nu:.4}");
    }
}
