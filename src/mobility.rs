//! Опорные значения свободной подвижности и её температурного коэффициента при 25 °C.
//! Историческая модель [`MobilityModel`] воспроизводит таблицу PHYSICS.md §3.3. В прогонах используется [`crate::mobility2`].
//! Здесь линейный график Фергюсона задаёт μ/μ₀ = exp(−K(N)·T%) (Ferguson 1964, Rodbard & Chrambach 1970).
//! Просеивание сменяется рептацией и высокомолекулярным плато (Ogston 1958, Lerman & Frisch 1982, Stellwagen 2009).

use crate::fmath::DetMath;
use crate::materials::agarose::rg_nm;

/// Свободнорешёточная подвижность дсДНК, см²/(В·с), 25 °C.
/// Литература: Stellwagen, Gelfi & Righetti, Biopolymers 42:687–703, 1997.
pub const MU0_CM2_V_S: f64 = 3.68e-4;

/// Температурный коэффициент подвижности, доля на °C. Литературный диапазон 1.8–2.5 %/°C (Stellwagen 2009).
pub const KT_PER_CELSIUS: f64 = 0.02;

/// Опорная температура для μ₀ по Stellwagen (1997) и для поправки проводимости камеры.
pub const T_REF_K: f64 = 298.15;

/// Калибровочные отношения μ/μ₀ для 1% агарозы, сгруппированные по длине ДНК.
pub const TARGET_MU100: f64 = 0.98;
pub const TARGET_MU300: f64 = 0.92;
/// Значение внутри диапазона 0.55–0.75 μ₀ для плато (Stellwagen 2009).
pub const TARGET_MU5000: f64 = 0.60;
pub const TARGET_MU23000: f64 = 0.533;

/// Предел коэффициента торможения для длинной ДНК в постоянном поле (Stellwagen 2009).
pub const K_MAX: f64 = 0.65;

/// Узел «просеивание → рептационный переход», п.о. Подгонка.
pub const SEG_SMALL_MAX_BP: f64 = 1000.0;
/// Узел «переход → плато», п.о. Подгонка.
pub const SEG_BRIDGE_MAX_BP: f64 = 5000.0;
/// Размер-якорь наклона плато, п.о. Подгонка.
pub const PLATEAU_ANCHOR_BP: f64 = 23000.0;

/// Опорная таблица (п.о., μ/μ₀) для 1% агарозы, PHYSICS.md §3.3.
pub const DOC_MOBILITY_TABLE: [(f64, f64); 9] = [
    (100.0, 0.980),
    (300.0, 0.920),
    (500.0, 0.850),
    (1000.0, 0.681),
    (3000.0, 0.627),
    (5000.0, 0.600),
    (10000.0, 0.568),
    (23000.0, 0.533),
    (100000.0, 0.522),
];

/// Коэффициенты кусочной кривой торможения, полученные из калибровочных отношений μ/μ₀.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MobilityModel {
    k_l: f64,
    k_q: f64,
    gamma: f64,
    k1: f64,
    k2: f64,
    s_r: f64,
}

impl MobilityModel {
    /// Определяет два коэффициента просеивания и непрерывно присоединяет переход и плато.
    pub fn calibrated() -> Self {
        let u = |n_bp: f64| rg_nm(n_bp) / rg_nm(SEG_SMALL_MAX_BP);
        let u1 = u(100.0);
        let u3 = u(300.0);

        let k100 = -TARGET_MU100.dln();
        let k300 = -TARGET_MU300.dln();

        // [[u1, u1²], [u3, u3²]] · [k_l, k_q]ᵀ = [K100, K300]; метод Крамера.
        let det = u1 * u3 * (u3 - u1);
        debug_assert!(det > 1e-6, "калибровочная система вырождена");
        let k_l = (k100 * u3 * u3 - k300 * u1 * u1) / det;
        let k_q = (u1 * k300 - u3 * k100) / det;
        debug_assert!(
            k_l > 0.0 && k_q > 0.0,
            "калибровка просеивания выродилась: k_l={k_l}, k_q={k_q}"
        );

        let k1 = k_l + k_q;
        let k2 = -TARGET_MU5000.dln();
        let gamma = (k2 / k1).dln() / (SEG_BRIDGE_MAX_BP / SEG_SMALL_MAX_BP).dln();

        let k23 = -TARGET_MU23000.dln();
        let s_r = (k23 - k2) / (PLATEAU_ANCHOR_BP / SEG_BRIDGE_MAX_BP).dln();

        debug_assert!(gamma > 0.0 && s_r > 0.0, "рептационные сегменты убывают");

        Self {
            k_l,
            k_q,
            gamma,
            k1,
            k2,
            s_r,
        }
    }

    /// Порядок коэффициентов: k_l, k_q, γ, K₁, K₂, s_r.
    pub fn coeffs(&self) -> [f64; 6] {
        [self.k_l, self.k_q, self.gamma, self.k1, self.k2, self.s_r]
    }

    /// Коэффициент торможения на 1% агарозы. Все три ветви непрерывны в точках стыка.
    pub fn retardation(&self, n_bp: f64) -> f64 {
        debug_assert!(n_bp > 0.0, "размер популяции должен быть положительным");
        if n_bp <= SEG_SMALL_MAX_BP {
            let u = rg_nm(n_bp) / rg_nm(SEG_SMALL_MAX_BP);
            self.k_l * u + self.k_q * u * u
        } else if n_bp <= SEG_BRIDGE_MAX_BP {
            self.k1 * (n_bp / SEG_SMALL_MAX_BP).dpowf(self.gamma)
        } else {
            (self.k2 + self.s_r * (n_bp / SEG_BRIDGE_MAX_BP).dln()).min(K_MAX)
        }
    }

    /// μ/μ0 при данной концентрации геля (в процентах).
    pub fn mobility_ratio(&self, n_bp: f64, percent_gel: f64) -> f64 {
        (-self.retardation(n_bp) * percent_gel).dexp()
    }

    /// Абсолютная подвижность, см²/(В·с), с температурной поправкой θ(T).
    pub fn mobility_cm2_v_s(&self, n_bp: f64, percent_gel: f64, temp_k: f64) -> f64 {
        let theta = 1.0 + KT_PER_CELSIUS * (temp_k - T_REF_K);
        MU0_CM2_V_S * self.mobility_ratio(n_bp, percent_gel) * theta.max(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn calibration_targets_hit_exactly() {
        // Цели, по которым решался калибратор, выполняются с машинной точностью.
        let m = MobilityModel::calibrated();
        for &(n, target) in &[
            (100.0, TARGET_MU100),
            (300.0, TARGET_MU300),
            (5000.0, TARGET_MU5000),
            (23000.0, TARGET_MU23000),
        ] {
            let got = m.mobility_ratio(n, 1.0);
            assert!(
                (got - target).abs() < 1e-12,
                "N={n}: {got:.12} против цели {target}"
            );
        }
    }

    #[test]
    fn doc_table_matches_engine() {
        // Правило документ↔код: таблица PHYSICS.md §3.3 — округлённые значения
        // движка; допуск 1.5e-3 покрывает округление до трёх знаков.
        let m = MobilityModel::calibrated();
        for &(n, mu_doc) in &DOC_MOBILITY_TABLE {
            let mu_engine = m.mobility_ratio(n, 1.0);
            assert!(
                (mu_engine - mu_doc).abs() < 1.5e-3,
                "μ({n}): движок {mu_engine:.5}, документ {mu_doc}"
            );
        }
    }

    #[test]
    fn retardation_monotone_nonnegative() {
        // Монотонный неубывающий рост K(N) на лог-сетке 50…200 000 п.о.,
        // строгие увеличения хотя бы до выхода на плато; K ≥ 0 всюду.
        let m = MobilityModel::calibrated();
        let mut prev = -1.0;
        let mut strict_growth_seen = false;
        let steps = 400;
        for i in 0..=steps {
            let n = 50.0 * (200000.0f64 / 50.0).dpowf(i as f64 / steps as f64);
            let k = m.retardation(n);
            assert!(k >= 0.0, "K({n}) = {k} отрицателен");
            assert!(k >= prev - 1e-15, "K({n}) = {k} меньше предыдущего {prev}");
            if k > prev + 1e-12 {
                strict_growth_seen = true;
            }
            prev = k;
        }
        assert!(strict_growth_seen, "K(N) должна расти хотя бы где-то");
    }

    #[test]
    fn mobility_never_exceeds_free_solution() {
        // Физический потолок μ ≤ μ0 на сетке размеров и концентраций.
        let m = MobilityModel::calibrated();
        for i in 0..=60 {
            let n = 50.0 * (200000.0f64 / 50.0).dpowf(i as f64 / 60.0);
            for &t in &[0.5, 0.8, 1.0, 1.4, 2.0] {
                let r = m.mobility_ratio(n, t);
                assert!(
                    r > 0.0 && r <= 1.0 + 1e-12,
                    "μ/μ0({n}, {t}%) = {r} вне (0, 1]"
                );
            }
        }
    }

    #[test]
    fn joints_continuous() {
        // Стыки режимов непрерывны по построению — проверяем численно.
        let m = MobilityModel::calibrated();
        let eps = 1e-6;
        for knot in [SEG_SMALL_MAX_BP, SEG_BRIDGE_MAX_BP] {
            let a = m.retardation(knot - eps);
            let b = m.retardation(knot);
            let c = m.retardation(knot + eps);
            assert!(
                (b - a).abs() < 1e-3 && (c - b).abs() < 1e-3,
                "стык {knot} рваный: {a:.6} / {b:.6} / {c:.6}"
            );
        }
    }

    #[test]
    fn temperature_correction_linear_and_positive() {
        // θ(T) линейна с коэффициентом KT_PER_CELSIUS и не обращает μ в ноль
        // в рабочем диапазоне температур.
        let m = MobilityModel::calibrated();
        let base = m.mobility_cm2_v_s(1000.0, 1.0, T_REF_K);
        let warmer = m.mobility_cm2_v_s(1000.0, 1.0, T_REF_K + 10.0);
        assert!(((warmer / base) - (1.0 + 10.0 * KT_PER_CELSIUS)).abs() < 1e-12);
        assert!(m.mobility_cm2_v_s(1000.0, 1.0, T_REF_K - 20.0) > 0.0);
    }
}
