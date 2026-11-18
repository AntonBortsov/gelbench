//! Электрическая цепь камеры при стабилизации напряжения, тока или мощности.
//! Гель и слой буфера над ним соединены параллельно, два электродных моста —
//! последовательно с этой парой. Поле в геле определяется падением напряжения
//! на параллельном участке. Геометрия и проводимость буфера задают сопротивления.

use crate::fmath::DetMath;

/// Приближённая относительная проводимость 1%-го агарозного геля (буфер = 1).
pub const KAPPA_GEL_OVER_BUFFER: f64 = 0.97;

/// Температурный коэффициент проводимости, 1/°C. Для водных растворов типичен диапазон 0,015–0,03.
pub const K_KAPPA_PER_C: f64 = 0.02;

/// Линейная температурная поправка к проводимости буфера при 25 °C. Температура задаётся в кельвинах.
#[inline]
pub fn kappa_at(kappa_25_s_m: f64, t_kelvin: f64) -> f64 {
    kappa_25_s_m * (1.0 + K_KAPPA_PER_C * (t_kelvin - crate::mobility::T_REF_K))
}

/// Геометрия цепи ячейки, м.
#[derive(Debug, Clone, Copy)]
pub struct CellGeometry {
    /// Длина геля вдоль направления миграции.
    pub l_x_m: f64,
    pub width_m: f64,
    pub t_gel_m: f64,
    /// Толщина слоя буфера НАД гелем (шунт).
    pub d_overlay_m: f64,
    /// Длина буферного моста электрод↔торец геля (с каждой стороны).
    pub l_bridge_m: f64,
    pub h_bridge_m: f64,
}

impl CellGeometry {
    /// Сопротивление геля, Ом.
    pub fn r_gel(&self, kappa_gel_s_m: f64) -> f64 {
        self.l_x_m / (kappa_gel_s_m * self.width_m * self.t_gel_m)
    }

    /// Сопротивление последовательных ячеек, Ом. Интеграл 1/σ берётся по значениям в центрах.
    /// Для полного геля n·dx должно равняться его длине. Срез профиля даёт сопротивление участка.
    pub fn r_gel_profile(&self, sigma_s_m: &[f64], dx_m: f64) -> f64 {
        assert!(!sigma_s_m.is_empty(), "пустой профиль σ");
        assert!(dx_m > 0.0);
        let integral: f64 = sigma_s_m
            .iter()
            .map(|s| {
                assert!(*s > 0.0, "σ ≤ 0 в профиле: обеднение среды");
                1.0 / s
            })
            .sum::<f64>()
            * dx_m;
        integral / (self.width_m * self.t_gel_m)
    }

    /// Сопротивление слоя буфера над гелем, Ом.
    pub fn r_overlay(&self, kappa_buf_s_m: f64) -> f64 {
        self.l_x_m / (kappa_buf_s_m * self.width_m * self.d_overlay_m)
    }

    /// Параллельная пара гель ∥ шунт, Ом.
    pub fn r_span(&self, kappa_buf_s_m: f64) -> f64 {
        let r_gel = self.r_gel(kappa_buf_s_m * KAPPA_GEL_OVER_BUFFER);
        let r_ov = self.r_overlay(kappa_buf_s_m);
        r_gel * r_ov / (r_gel + r_ov)
    }

    /// Сопротивление параллельной пары по заданным R геля и проводимости буфера над ним.
    pub fn r_span_with_gel(&self, r_gel_ohm: f64, kappa_buf_s_m: f64) -> f64 {
        assert!(r_gel_ohm > 0.0, "R_gel из профиля обязан быть положителен");
        let r_ov = self.r_overlay(kappa_buf_s_m);
        r_gel_ohm * r_ov / (r_gel_ohm + r_ov)
    }

    /// Сопротивление по геометрии камеры и эффективному удельному сопротивлению геля, Ом.
    pub fn r_gel_from_resistivity(&self, rho_ohm_m: f64) -> f64 {
        assert!(rho_ohm_m > 0.0, "удельное сопротивление обязано быть > 0");
        self.l_x_m * rho_ohm_m / (self.width_m * self.t_gel_m)
    }

    /// Один буферный мост, Ом (в цепи их два, последовательно).
    pub fn r_bridge(&self, kappa_buf_s_m: f64) -> f64 {
        self.l_bridge_m / (kappa_buf_s_m * self.width_m * self.h_bridge_m)
    }

    /// Полное сопротивление ячейки, Ом.
    pub fn r_total(&self, kappa_buf_s_m: f64) -> f64 {
        self.r_span(kappa_buf_s_m) + 2.0 * self.r_bridge(kappa_buf_s_m)
    }

    /// Доля тока в буферном шунте. Последовательные мосты не участвуют в делении тока.
    pub fn shunt_fraction(&self, kappa_buf_s_m: f64) -> f64 {
        let r_gel = self.r_gel(kappa_buf_s_m * KAPPA_GEL_OVER_BUFFER);
        let r_ov = self.r_overlay(kappa_buf_s_m);
        r_gel / (r_gel + r_ov)
    }
}

/// Приближённая геометрия типовой мини-камеры класса Mini-Sub.
pub const MINI_SUB_REFERENCE: CellGeometry = CellGeometry {
    l_x_m: 0.10,
    width_m: 0.10,
    t_gel_m: 0.006,
    d_overlay_m: 0.003,
    l_bridge_m: 0.010,
    h_bridge_m: 0.005,
};

/// Уставка источника питания.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SourceMode {
    /// Стабилизация напряжения, В.
    ConstantVoltage(f64),
    /// Стабилизация тока, А.
    ConstantCurrent(f64),
    /// Стабилизация мощности, Вт.
    ConstantPower(f64),
}

/// Напряжение и ток на клеммах, суммарная мощность и поле на гелевом участке.
#[derive(Debug, Clone, Copy)]
pub struct OperatingPoint {
    pub u_v: f64,
    pub i_a: f64,
    pub p_w: f64,
    pub e_gel_v_per_m: f64,
    pub e_gel_v_per_cm: f64,
}

/// Решает цепь для однородного геля при заданных κ₂₅ буфера и температуре.
pub fn solve(
    mode: SourceMode,
    geom: &CellGeometry,
    kappa_25_s_m: f64,
    t_kelvin: f64,
) -> OperatingPoint {
    let kappa = kappa_at(kappa_25_s_m, t_kelvin);
    solve_with_r_gel(mode, geom, kappa, r_gel_of(geom, kappa))
}

fn r_gel_of(geom: &CellGeometry, kappa: f64) -> f64 {
    geom.r_gel(kappa * KAPPA_GEL_OVER_BUFFER)
}

/// Заданное сопротивление геля может быть получено из профиля σ(x).
/// Шунт и мосты используют проводимость буфера при текущей температуре.
pub fn solve_with_r_gel(
    mode: SourceMode,
    geom: &CellGeometry,
    kappa_buf_s_m: f64,
    r_gel_ohm: f64,
) -> OperatingPoint {
    let r_ov = geom.r_overlay(kappa_buf_s_m);
    let r_span = r_gel_ohm * r_ov / (r_gel_ohm + r_ov);
    let r_total = r_span + 2.0 * geom.r_bridge(kappa_buf_s_m);
    let u_v = match mode {
        SourceMode::ConstantVoltage(u) => u,
        SourceMode::ConstantCurrent(i) => i * r_total,
        SourceMode::ConstantPower(p) => (p * r_total).sqrt(),
    };
    let i_a = u_v / r_total;
    let divider = r_span / r_total;
    let e_gel_v_per_m = u_v * divider / geom.l_x_m;
    OperatingPoint {
        u_v,
        i_a,
        p_w: u_v * i_a,
        e_gel_v_per_m,
        e_gel_v_per_cm: e_gel_v_per_m / 100.0,
    }
}

/// ЭДС жидкостного соединения действует последовательно в гелевой ветви.
/// Через буферный шунт течёт ток Gш·U, через гель — Gг·(U − Uj), где U — падение на пролёте.
pub fn solve_with_r_gel_and_junction(
    mode: SourceMode,
    geom: &CellGeometry,
    kappa_buf_s_m: f64,
    r_gel_ohm: f64,
    u_junction_v: f64,
) -> OperatingPoint {
    assert!(r_gel_ohm > 0.0, "R_геля ≤ 0: обеднение среды");
    let g_gel = 1.0 / r_gel_ohm;
    let g_ov = 1.0 / geom.r_overlay(kappa_buf_s_m);
    let g_sum = g_gel + g_ov;
    let r_bridges = 2.0 * geom.r_bridge(kappa_buf_s_m);

    let (u_span, i_a) = match mode {
        SourceMode::ConstantVoltage(u) => {
            let u_span = (u + r_bridges * g_gel * u_junction_v) / (r_bridges * g_sum + 1.0);
            let i_a = u_span * g_sum - g_gel * u_junction_v;
            (u_span, i_a)
        }
        SourceMode::ConstantCurrent(i) => ((i + g_gel * u_junction_v) / g_sum, i),
        SourceMode::ConstantPower(p) => {
            // U_клемм(I) = r_bridges·I + (I + G_г·U_j)/g_sum; P = I·U_клемм.
            let a_c = r_bridges + 1.0 / g_sum;
            let b_c = g_gel * u_junction_v / g_sum;
            let i_a = (-b_c + (b_c * b_c + 4.0 * a_c * p).sqrt()) / (2.0 * a_c);
            let u_span = (i_a + g_gel * u_junction_v) / g_sum;
            (u_span, i_a)
        }
    };
    let u_v = r_bridges * i_a + u_span;
    // Падение на гелевой ветви без ЭДС — это I_гель·R_гель; среднее поле.
    let i_gel = (u_span - u_junction_v) * g_gel;
    let e_gel_v_per_m = i_gel * r_gel_ohm / geom.l_x_m;
    OperatingPoint {
        u_v,
        i_a,
        p_w: u_v * i_a,
        e_gel_v_per_m,
        e_gel_v_per_cm: e_gel_v_per_m / 100.0,
    }
}

/// Проводимость буфера при температуре измерения, См/м, по току и напряжению на клеммах.
/// Все сопротивления пропорциональны 1/κ. Точность оценки ограничена точностью геометрии.
pub fn kappa_from_reading(u_v: f64, i_a: f64, geom: &CellGeometry, t_kelvin: f64) -> f64 {
    let _ = t_kelvin;
    let c_geom = geom.r_total(1.0) * 1.0;
    i_a * c_geom / u_v
}

/// Расстояние от дна мостового канала до нижней поверхности геля, м. Приближённое значение для пресетов.
pub const TRAY_SURFACE_PLACEHOLDER_M: f64 = 0.003;

/// Порог толщины проводящей плёнки, м.
const MIN_WET_FILM_M: f64 = 1e-6;

/// Контакт определяется слоем буфера над гелем и погружением мостов.
/// При высыхании шунта ток идёт через гель, при обрыве моста ток прекращается.
/// О контактных мостах см. Kaplan, Russo & Wilcox, Anal. Biochem. (1976).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LevelContact {
    Wet { d_overlay_m: f64, h_bridge_m: f64 },
    DryShunt { h_bridge_m: f64 },
    Open,
}

/// Проводящее сечение моста ограничено высотой его горловины.
pub fn contact_at_level(level_m: f64, tray_surface_m: f64, geom: &CellGeometry) -> LevelContact {
    assert!(
        level_m.is_finite() && tray_surface_m.is_finite(),
        "уровень и высота поддона обязаны быть конечными"
    );
    if level_m <= 0.0 {
        return LevelContact::Open;
    }
    let gel_top_m = tray_surface_m + geom.t_gel_m;
    let h_bridge_m = level_m.min(geom.h_bridge_m);
    let d_overlay_m = level_m - gel_top_m;
    if d_overlay_m >= MIN_WET_FILM_M {
        LevelContact::Wet {
            d_overlay_m,
            h_bridge_m,
        }
    } else {
        LevelContact::DryShunt { h_bridge_m }
    }
}

impl CellGeometry {
    /// Геометрия при фактическом уровне буфера. При разомкнутой цепи возвращает `None`.
    pub fn with_level(&self, level_m: f64, tray_surface_m: f64) -> Option<CellGeometry> {
        match contact_at_level(level_m, tray_surface_m, self) {
            LevelContact::Wet {
                d_overlay_m,
                h_bridge_m,
            } => Some(CellGeometry {
                d_overlay_m,
                h_bridge_m,
                ..*self
            }),
            // Сухой шунт представляем большим конечным сопротивлением.
            LevelContact::DryShunt { h_bridge_m } => Some(CellGeometry {
                d_overlay_m: MIN_WET_FILM_M * 1e-3,
                h_bridge_m,
                ..*self
            }),
            LevelContact::Open => None,
        }
    }
}

/// При обрыве постоянное напряжение сохраняется на клеммах, ток и поле в геле равны нулю.
/// Для стабилизации I или P обрыв даёт нулевые показания. Предел напряжения источника не моделируется.
pub fn solve_at_level(
    mode: SourceMode,
    geom: &CellGeometry,
    kappa_25_s_m: f64,
    t_kelvin: f64,
    level_m: f64,
    tray_surface_m: f64,
) -> OperatingPoint {
    match geom.with_level(level_m, tray_surface_m) {
        Some(g) => solve(mode, &g, kappa_25_s_m, t_kelvin),
        None => {
            let u_v = match mode {
                SourceMode::ConstantVoltage(u) => u,
                _ => 0.0,
            };
            OperatingPoint {
                u_v,
                i_a: 0.0,
                p_w: 0.0,
                e_gel_v_per_m: 0.0,
                e_gel_v_per_cm: 0.0,
            }
        }
    }
}

/// Постоянная Фарадея, Кл/моль (CODATA 2018).
pub const FARADAY_C_PER_MOL: f64 = 96485.33212;

const R_GAS_J_PER_MOL_K: f64 = 8.314462618;
const P_ATM_PA: f64 = 101325.0;

/// Объём газа по закону Фарадея при 25 °C и 1 атм, мл. Один O₂ требует четыре электрона,
/// H₂ — два. Принят единичный выход по току на Pt (PRX Energy 4, 013011, 2025).
pub fn faraday_gas_volume_ml(charge_c: f64, electrons_per_molecule: f64) -> f64 {
    assert!(charge_c >= 0.0, "заряд отрицательный");
    assert!(electrons_per_molecule > 0.0, "z обязано быть положительным");
    let n_mol = charge_c / (electrons_per_molecule * FARADAY_C_PER_MOL);
    n_mol * R_GAS_J_PER_MOL_K * crate::mobility::T_REF_K / P_ATM_PA * 1e6
}

/// Проводимость пузырьковой дисперсии по Бруггману. Для φ ≤ 0,3 альтернативные модели
/// различаются примерно на 7% или меньше (Sci Rep 15, s41598-025-25842-3; arXiv:1304.2701).
#[inline]
pub fn kappa_with_void_fraction(kappa_s_m: f64, void_fraction: f64) -> f64 {
    assert!(
        (0.0..=1.0).contains(&void_fraction),
        "объёмная доля газа φ вне [0; 1]: {void_fraction}"
    );
    kappa_s_m * (1.0 - void_fraction).dpowf(1.5)
}

/// Объёмная доля газа при покрытии электрода θ. Приближение φ = 2θ/3.
#[inline]
pub fn void_fraction_from_coverage(theta: f64) -> f64 {
    assert!(
        (0.0..=1.0).contains(&theta),
        "покрытие Θ вне [0; 1]: {theta}"
    );
    2.0f64 / 3.0 * theta
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::buffers::KAPPA_TAE_S_M;

    const T_ROOM: f64 = 298.0;

    #[test]
    fn r_gel_profile_matches_uniform_scalar() {
        // Однородный профиль σ(x): r_gel_profile обязан совпасть со скалярным
        // r_gel до машинного округления — точка входа M3 непрерывно продолжает
        // M2.
        let kappa = KAPPA_TAE_S_M * 0.97;
        let n = 64;
        let dx = MINI_SUB_REFERENCE.l_x_m / n as f64;
        let sigma = vec![kappa; n];
        let r_prof = MINI_SUB_REFERENCE.r_gel_profile(&sigma, dx);
        let r_scal = MINI_SUB_REFERENCE.r_gel(kappa);
        assert!(
            (r_prof - r_scal).abs() / r_scal < 1e-12,
            "r_gel_profile {r_prof} против r_gel {r_scal}"
        );
    }

    #[test]
    fn solve_with_r_gel_reproduces_m2_solve_on_uniform_medium() {
        // solve_with_r_gel с R из однородного профиля (σ_геля = 0.97·κ)
        // воспроизводит старый solve() по току и полю.
        let kappa = KAPPA_TAE_S_M;
        let op_old = solve(
            SourceMode::ConstantVoltage(150.0),
            &MINI_SUB_REFERENCE,
            kappa,
            T_ROOM,
        );
        let n = 32;
        let dx = MINI_SUB_REFERENCE.l_x_m / n as f64;
        let sigma = vec![kappa_at(kappa, T_ROOM) * KAPPA_GEL_OVER_BUFFER; n];
        let op_new = solve_with_r_gel(
            SourceMode::ConstantVoltage(150.0),
            &MINI_SUB_REFERENCE,
            kappa_at(kappa, T_ROOM),
            MINI_SUB_REFERENCE.r_gel_profile(&sigma, dx),
        );
        assert!((op_old.u_v - op_new.u_v).abs() < 1e-9);
        assert!((op_old.i_a - op_new.i_a).abs() / op_old.i_a < 1e-12);
        assert!((op_old.e_gel_v_per_cm - op_new.e_gel_v_per_cm).abs() < 1e-10);
    }

    #[test]
    fn limits_no_shunt_and_no_bridge() {
        // Без слоя над гелем (d_over → 0 ⟹ R_overlay → ∞) и без мостов
        // (L_br → 0): вся цепь — чистый R_gel, E_gel = U/L_x точно.
        let bare = CellGeometry {
            d_overlay_m: 1e-12,
            l_bridge_m: 1e-12,
            ..MINI_SUB_REFERENCE
        };
        let op = solve(
            SourceMode::ConstantVoltage(100.0),
            &bare,
            KAPPA_TAE_S_M,
            T_ROOM,
        );
        let e_expected = 100.0 / bare.l_x_m;
        assert!(
            (op.e_gel_v_per_m - e_expected).abs() / e_expected < 1e-9,
            "E_gel = {} против U/L_x = {e_expected}",
            op.e_gel_v_per_m
        );
        // Полная геометрия (мосты + шунт): делитель и шунт занижают E_gel
        // заметно ниже U/L_x — именно поэтому якорь №1 фиксирует E_gel,
        // а не напряжение на клеммах.
        let op_full = solve(
            SourceMode::ConstantVoltage(100.0),
            &MINI_SUB_REFERENCE,
            KAPPA_TAE_S_M,
            T_ROOM,
        );
        assert!(
            op_full.e_gel_v_per_m < e_expected * 0.99,
            "шунт+мосты обязаны занижать E_gel, получено {}",
            op_full.e_gel_v_per_m
        );
    }

    #[test]
    fn shunt_fraction_matches_review_3_4() {
        // Ревью §3.4 (гель 6 мм): доля тока мимо геля 26% при слое 2 мм
        // и 63% при 10 мм. Точный расчёт с коэффициентом 0.97 даёт
        // 25.57% и 63.21%; документ печатал целые проценты, поэтому
        // страж по половине последнего напечатанного разряда (±0.5 п.п.)
        // накрывает и систематику 0.97.
        let base = CellGeometry {
            t_gel_m: 0.006,
            ..MINI_SUB_REFERENCE
        };
        for (d_over_m, doc_share) in [(0.002, 0.26), (0.010, 0.63)] {
            let g = CellGeometry {
                d_overlay_m: d_over_m,
                ..base
            };
            let share = g.shunt_fraction(KAPPA_TAE_S_M);
            assert!(
                (share - doc_share).abs() < 0.005,
                "d_over = {} мм: доля шунта {:.4}, документ {}",
                d_over_m * 1e3,
                share,
                doc_share
            );
            // И точная алгебра: share = d/(d + 0.97·t_g).
            let exact = d_over_m / (d_over_m + KAPPA_GEL_OVER_BUFFER * base.t_gel_m);
            assert!((share - exact).abs() < 1e-12);
        }
    }

    #[test]
    fn modes_agree_at_same_current() {
        // При одном рабочем токе все три режима обязаны дать одну и ту же
        // рабочую точку. Алгебра без аппроксимаций — допуск машинный.
        let i0 = 0.05_f64;
        let i_mode = solve(
            SourceMode::ConstantCurrent(i0),
            &MINI_SUB_REFERENCE,
            KAPPA_TAE_S_M,
            T_ROOM,
        );
        let u_mode = solve(
            SourceMode::ConstantVoltage(i_mode.u_v),
            &MINI_SUB_REFERENCE,
            KAPPA_TAE_S_M,
            T_ROOM,
        );
        let p_mode = solve(
            SourceMode::ConstantPower(i_mode.p_w),
            &MINI_SUB_REFERENCE,
            KAPPA_TAE_S_M,
            T_ROOM,
        );
        for (name, a, b) in [("U", &u_mode, &i_mode), ("P", &p_mode, &i_mode)] {
            assert!(
                (a.i_a - b.i_a).abs() / b.i_a < 1e-12,
                "{name}-режим: ток разошёлся"
            );
            assert!(
                (a.e_gel_v_per_cm - b.e_gel_v_per_cm).abs() / b.e_gel_v_per_cm < 1e-12,
                "{name}-режим: E_gel разошёлся"
            );
            assert!(
                (a.p_w - b.p_w).abs() / b.p_w < 1e-12,
                "{name}-режим: мощность разошлась"
            );
        }
    }

    #[test]
    fn feedback_signs_6a_6b_6c_by_finite_difference() {
        // Якоря №6a/6b/6c — знак обратной связи мощность↔температура,
        // конечной разностью по T на реализованных функциях:
        // U-режим dQ/dT > 0 (разгон), I-режим dQ/dT < 0 (устойчиво),
        // P-режим dQ/dT = 0 (нейтрально).
        let dt = 10.0;
        let q = |mode: SourceMode| {
            let cold = solve(mode, &MINI_SUB_REFERENCE, KAPPA_TAE_S_M, 298.0);
            let hot = solve(mode, &MINI_SUB_REFERENCE, KAPPA_TAE_S_M, 298.0 + dt);
            (hot.p_w - cold.p_w) / dt
        };
        let dq_u = q(SourceMode::ConstantVoltage(100.0));
        let dq_i = q(SourceMode::ConstantCurrent(0.05));
        let dq_p = q(SourceMode::ConstantPower(5.0));
        assert!(dq_u > 0.0, "U-режим: dQ/dT = {dq_u} — ожидался разгон (>0)");
        assert!(
            dq_i < 0.0,
            "I-режим: dQ/dT = {dq_i} — ожидалась устойчивость (<0)"
        );
        assert!(
            dq_p.abs() < 1e-12,
            "P-режим: dQ/dT = {dq_p} — мощность обязана быть константой"
        );
    }

    #[test]
    fn anchor7_mini_sub_current_order_of_magnitude() {
        // Якорь №7 (лабораторный фольклор, оксфордский протокол Bio-Rad
        // Mini-Sub): 1× TAE, 100 В → ток 50–100 мА. Это ВАЛИДАЦИЯ, не
        // калибровка: геометрия — заглушка MINI_SUB_REFERENCE, κ — верхняя
        // оценка (предельные λ° без поправки на ионную силу), поэтому
        // предсказание ~110 мА выше окна ЗАВЕДОМО. Подгонять запрещено;
        // после измерения κ окно ужесточается до 50–100 мА. Сейчас
        // проверяем порядок величины: [40; 250] мА.
        let op = solve(
            SourceMode::ConstantVoltage(100.0),
            &MINI_SUB_REFERENCE,
            KAPPA_TAE_S_M,
            T_ROOM,
        );
        let ma = op.i_a * 1e3;
        println!("якорь №7: I(100 В, TAE, мини-камера) = {ma:.1} мА");
        assert!(
            (40.0..=250.0).contains(&ma),
            "якорь №7: I = {ma:.1} мА вне окна порядка величины [40; 250] мА"
        );
    }

    #[test]
    fn inverse_kappa_roundtrip_exact() {
        // Обратная задача: κ, восстановленная из посчитанного тока,
        // совпадает с исходной до машинного нуля.
        for t in [298.0, 303.0, 308.0] {
            let op = solve(
                SourceMode::ConstantVoltage(90.0),
                &MINI_SUB_REFERENCE,
                KAPPA_TAE_S_M,
                t,
            );
            let restored = kappa_from_reading(op.u_v, op.i_a, &MINI_SUB_REFERENCE, t);
            let expected = kappa_at(KAPPA_TAE_S_M, t);
            assert!(
                (restored - expected).abs() / expected < 1e-12,
                "T={t}: восстановлено {restored}, исходно {expected}"
            );
        }
    }

    #[test]
    fn u_mode_field_independent_of_temperature() {
        // Инвариант, делающий порог разгона замкнутой формулой: все
        // сопротивления ∝ 1/κ, поэтому делитель R_span/R_total от κ не
        // зависит — в U-режиме E_gel не зависит от температуры. Именно
        // поэтому q_v(T) = κ(T)·E² растёт при фиксированном поле.
        let e_ref = solve(
            SourceMode::ConstantVoltage(100.0),
            &MINI_SUB_REFERENCE,
            KAPPA_TAE_S_M,
            298.0,
        )
        .e_gel_v_per_cm;
        for t in [298.0, 320.0, 350.0] {
            let op = solve(
                SourceMode::ConstantVoltage(100.0),
                &MINI_SUB_REFERENCE,
                KAPPA_TAE_S_M,
                t,
            );
            assert!(
                (op.e_gel_v_per_cm - e_ref).abs() < 1e-12,
                "T={t}: E_gel {} против {e_ref}",
                op.e_gel_v_per_cm
            );
        }
    }

    #[test]
    fn anchor1_field_is_reachable_and_consistent() {
        // Якорь №1 не сломан: цепь позволяет получить E_gel = 5 В/см,
        // и прямая функция воспроизводит поле из найденного U точно.
        // Сам путь 2.98 см держит интеграционный тест m1_anchors
        // (anchor1_front_position_5kb) — он остаётся зелёным.
        let target_e_v_per_m = 500.0;
        let kappa = kappa_at(KAPPA_TAE_S_M, T_ROOM);
        let divider = MINI_SUB_REFERENCE.r_span(kappa) / MINI_SUB_REFERENCE.r_total(kappa);
        let u_needed = target_e_v_per_m * MINI_SUB_REFERENCE.l_x_m / divider;
        let op = solve(
            SourceMode::ConstantVoltage(u_needed),
            &MINI_SUB_REFERENCE,
            KAPPA_TAE_S_M,
            T_ROOM,
        );
        assert!(
            (op.e_gel_v_per_m - target_e_v_per_m).abs() / target_e_v_per_m < 1e-12,
            "E_gel = {} против цели {target_e_v_per_m}",
            op.e_gel_v_per_m
        );
    }
}
