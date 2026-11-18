//! Состояние камеры связывает уровень буфера и покрытие электродов с полем в геле.
//! Пузыри меняют сопротивление мостов по закону Бруггмана. Покрытие растёт с прошедшим
//! зарядом, а отрыв пузырей моделируется пуассоновским процессом. Параметры роста
//! и отрыва для горизонтальной платиновой проволоки не калиброваны.

use crate::circuit::{
    contact_at_level, kappa_at, void_fraction_from_coverage, CellGeometry, LevelContact,
    OperatingPoint, SourceMode,
};
use crate::fmath::DetMath;
use crate::rng::Rng;

/// Mini-Sub Cell GT (Bio-Rad M1704400B, Bulletin 6181). Длина мостов выведена
/// из расстояния между электродами и длины геля. Высота горловины моста оценочная.
pub const MINI_SUB_GT: CellGeometry = CellGeometry {
    l_x_m: 0.07,
    width_m: 0.07,
    t_gel_m: 0.006,
    d_overlay_m: 0.003,
    l_bridge_m: 0.0375,
    h_bridge_m: 0.005,
};

/// Объём буфера базы Mini-Sub Cell GT, л (каталог: 265–320 мл, ~270).
pub const MINI_SUB_GT_BUFFER_L: f64 = 0.270;

/// Прирост покрытия на кулон прошедшего заряда.
pub const THETA_GROWTH_PER_COULOMB: f64 = 0.6;
/// Частота отрыва пузырей, 1/с.
pub const DETACH_RATE_PER_S: f64 = 0.5;
/// Доля покрытия, остающаяся после отрыва.
pub const DETACH_KEEP: f64 = 0.3;

#[derive(Debug, Clone, Copy)]
pub struct ChamberSpec {
    /// Номинальная геометрия (горловина мостов полная).
    pub geometry: CellGeometry,
    /// Высота поддона над дном мостового канала, м.
    pub tray_surface_m: f64,
    /// Уставка на электродах. Падение напряжения на геле определяется геометрией цепи.
    pub source: SourceMode,
    /// Проводимость буфера при 25 °C, См/м.
    pub kappa_25_s_m: f64,
    /// Стартовый уровень буфера ℓ, м (от дна мостового канала).
    pub level_m: f64,
    pub bubble_seed: u64,
}

#[derive(Debug, Clone)]
pub struct ChamberState {
    pub spec: ChamberSpec,
    /// Уровень от дна мостового канала, м. Изменяется внешним вызовом при доливе или испарении.
    pub level_m: f64,
    /// Доля площади электрода под пузырями.
    pub theta_anode: f64,
    pub theta_cathode: f64,
    rng: Rng,
    /// Полный заряд через камеру, Кл.
    pub charge_c: f64,
    pub last: OperatingPoint,
}

impl ChamberState {
    pub fn new(spec: ChamberSpec) -> Self {
        let mut me = Self {
            level_m: spec.level_m,
            theta_anode: 0.0,
            theta_cathode: 0.0,
            rng: Rng::new(spec.bubble_seed),
            charge_c: 0.0,
            last: OperatingPoint {
                u_v: 0.0,
                i_a: 0.0,
                p_w: 0.0,
                e_gel_v_per_m: 0.0,
                e_gel_v_per_cm: 0.0,
            },
            spec,
        };
        me.last = me.solve_now(298.15);
        me
    }

    /// Решение цепи для текущих уровня буфера и покрытий электродов.
    pub fn solve_now(&self, t_kelvin: f64) -> OperatingPoint {
        self.solve_now_with_gel(t_kelvin, None)
    }

    /// Заданное удельное сопротивление геля (Ом·м) берётся при рабочей температуре.
    /// `None` использует проводимость однородного буфера с поправкой на агарозу.
    pub fn solve_now_with_gel(
        &self,
        t_kelvin: f64,
        gel_resistivity_ohm_m: Option<f64>,
    ) -> OperatingPoint {
        let geom = match self
            .spec
            .geometry
            .with_level(self.level_m, self.spec.tray_surface_m)
        {
            Some(g) => g,
            None => {
                let u_v = match self.spec.source {
                    SourceMode::ConstantVoltage(u) => u,
                    _ => 0.0,
                };
                return OperatingPoint {
                    u_v,
                    i_a: 0.0,
                    p_w: 0.0,
                    e_gel_v_per_m: 0.0,
                    e_gel_v_per_cm: 0.0,
                };
            }
        };
        let kappa = kappa_at(self.spec.kappa_25_s_m, t_kelvin);
        // Удельное сопротивление из петли среды уже учитывает температуру.
        let r_span = match gel_resistivity_ohm_m {
            Some(rho_at_t) => geom.r_span_with_gel(geom.r_gel_from_resistivity(rho_at_t), kappa),
            None => geom.r_span(kappa),
        };
        let m_a = (1.0 - void_fraction_from_coverage(self.theta_anode)).dpowf(1.5);
        let m_c = (1.0 - void_fraction_from_coverage(self.theta_cathode)).dpowf(1.5);
        let r_br = geom.r_bridge(kappa);
        let r_total = r_span + r_br / m_a + r_br / m_c;
        let (u_v, i_a) = match self.spec.source {
            SourceMode::ConstantVoltage(u) => (u, u / r_total),
            SourceMode::ConstantCurrent(i) => (i * r_total, i),
            SourceMode::ConstantPower(p) => {
                let i = (p / r_total).sqrt();
                (i * r_total, i)
            }
        };
        let u_span = i_a * r_span;
        let e_gel_v_per_m = u_span / geom.l_x_m;
        OperatingPoint {
            u_v,
            i_a,
            p_w: u_v * i_a,
            e_gel_v_per_m,
            e_gel_v_per_cm: e_gel_v_per_m / 100.0,
        }
    }

    /// Вычисляет рабочую точку в начале шага, затем обновляет покрытия по прошедшему заряду.
    pub fn step(&mut self, dt_s: f64, t_kelvin: f64) -> OperatingPoint {
        self.step_with_gel(dt_s, t_kelvin, None)
    }

    /// Шаг с выбором удельного сопротивления геля, как в [`Self::solve_now_with_gel`].
    pub fn step_with_gel(
        &mut self,
        dt_s: f64,
        t_kelvin: f64,
        gel_resistivity_ohm_m: Option<f64>,
    ) -> OperatingPoint {
        let op = self.solve_now_with_gel(t_kelvin, gel_resistivity_ohm_m);
        let dq = op.i_a * dt_s;
        self.charge_c += dq;
        // Анод первым получает случайную выборку, что фиксирует порядок для воспроизводимых прогонов.
        self.theta_anode = (self.theta_anode + THETA_GROWTH_PER_COULOMB * dq).min(1.0);
        self.theta_cathode = (self.theta_cathode + THETA_GROWTH_PER_COULOMB * dq).min(1.0);
        let p_detach = (DETACH_RATE_PER_S * dt_s).min(1.0);
        if self.rng.next_f64() < p_detach {
            self.theta_anode *= DETACH_KEEP;
        }
        if self.rng.next_f64() < p_detach {
            self.theta_cathode *= DETACH_KEEP;
        }
        self.last = op;
        op
    }

    pub fn level_contact(&self) -> LevelContact {
        contact_at_level(self.level_m, self.spec.tray_surface_m, &self.spec.geometry)
    }

    /// Мощность гелевой ветви в ваттах для локального теплового профиля.
    /// Баланс тепла всей ванны использует полную мощность.
    pub fn gel_power_w_at(
        &self,
        op: &OperatingPoint,
        t_kelvin: f64,
        gel_resistivity_ohm_m: Option<f64>,
    ) -> f64 {
        let Some(g) = self
            .spec
            .geometry
            .with_level(self.level_m, self.spec.tray_surface_m)
        else {
            return 0.0;
        };
        let r_gel = match gel_resistivity_ohm_m {
            Some(rho_at_t) => g.r_gel_from_resistivity(rho_at_t),
            None => g.r_gel(
                kappa_at(self.spec.kappa_25_s_m, t_kelvin) * crate::circuit::KAPPA_GEL_OVER_BUFFER,
            ),
        };
        let u_span = op.e_gel_v_per_m * g.l_x_m;
        (u_span * u_span / r_gel).max(0.0).min(op.p_w.max(0.0))
    }

    /// Оценка перегрева середины геля, К. Мощность пролёта распределена по объёму пластины
    /// с двусторонней естественной конвекцией. Коэффициент теплоотдачи оценочный.
    pub fn melt_delta_t_estimate_k(&self, op: &OperatingPoint) -> f64 {
        let g = match self
            .spec
            .geometry
            .with_level(self.level_m, self.spec.tray_surface_m)
        {
            Some(g) => g,
            None => return 0.0,
        };
        let kappa = kappa_at(self.spec.kappa_25_s_m, 298.15);
        let r_span = g.r_span(kappa);
        let volume = g.l_x_m * g.width_m * g.t_gel_m;
        let q_v = op.i_a * op.i_a * r_span / volume;
        let t = g.t_gel_m;
        const H_EST_W_M2K: f64 = 15.0; // оценка (вода/воздух между)
        q_v * (t * t / (8.0 * crate::thermal::K_GEL) + t / (2.0 * H_EST_W_M2K))
    }
}

/// Порог риска размягчения 65 °C (Thermo Tip 8). Обычная агароза плавится при 85–95 °C.
pub const MELT_RISK_K: f64 = 338.15;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuit::MINI_SUB_REFERENCE;

    const KAPPA: f64 = crate::buffers::KAPPA_TAE_S_M;

    fn spec(level_m: f64, seed: u64) -> ChamberSpec {
        ChamberSpec {
            geometry: MINI_SUB_REFERENCE,
            tray_surface_m: 0.003,
            source: SourceMode::ConstantVoltage(100.0),
            kappa_25_s_m: KAPPA,
            level_m,
            bubble_seed: seed,
        }
    }

    #[test]
    fn bubbles_jitter_inside_anchor19_window_and_deterministic() {
        // №19 в динамике: ток гуляет в [0.45·I₀; I₀], бит-в-бит по сиду.
        let mut a = ChamberState::new(spec(0.012, 7));
        let mut b = ChamberState::new(spec(0.012, 7));
        let i0 = a.solve_now(298.15).i_a;
        let mut min_ratio = f64::INFINITY;
        for _ in 0..2000 {
            let opa = a.step(0.5, 298.15);
            let opb = b.step(0.5, 298.15);
            assert_eq!(opa.i_a.to_bits(), opb.i_a.to_bits(), "сид не детерминирует");
            let r = opa.i_a / i0;
            assert!(r <= 1.0 + 1e-12, "ток выше номинала: {r}");
            assert!(r >= 0.45, "ток ниже пола №19: {r:.3}");
            min_ratio = min_ratio.min(r);
        }
        assert!(min_ratio < 0.98, "джиттер не виден: min {min_ratio:.3}");
        // Другой сид — другая траектория (но то же окно).
        let mut c = ChamberState::new(spec(0.012, 8));
        let mut differs = false;
        for _ in 0..50 {
            let oa = ChamberState::new(spec(0.012, 7)).step(0.5, 298.15);
            let oc = c.step(0.5, 298.15);
            differs |= oa.i_a.to_bits() != oc.i_a.to_bits();
        }
        assert!(differs, "разные сиды дали одну траекторию");
    }

    #[test]
    fn open_circuit_gives_zero_current_and_refill_restores() {
        // №16 на состоянии: ℓ ≤ 0 → I = 0; долив возвращает рабочую точку.
        let mut ch = ChamberState::new(spec(0.0, 1));
        assert_eq!(ch.step(0.5, 298.15).i_a, 0.0);
        ch.level_m = 0.012;
        let op = ch.step(0.5, 298.15);
        assert!(op.i_a > 0.05, "долив не восстановил ток: {}", op.i_a);
    }

    #[test]
    fn vendor_preset_dimensions_are_catalog() {
        // §1.1: межэлектродный путь = гель + два моста = 14.5 см.
        let total = MINI_SUB_GT.l_x_m + 2.0 * MINI_SUB_GT.l_bridge_m;
        assert!((total - 0.145).abs() < 1e-12);
        assert_eq!(MINI_SUB_GT.t_gel_m, 0.006);
    }
}
