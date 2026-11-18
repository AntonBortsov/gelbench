//! Transverse heat-transfer geometry, boundary conditions and buffer conductivity feedback.

use super::heat_transfer::tray_series_heat_transfer;
use super::{
    bisect, h_natural_convection, nusselt_horizontal_layer, rayleigh_horizontal_layer, Medium,
    DEFAULT_BUFFER_OVERLAY_M, DEFAULT_TRAY_K_W_M_K, DEFAULT_TRAY_THICKNESS_M, WATER,
};

/// Краевое условие на торцах пластины (y = ±W/2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeKind {
    /// Пришпилен: T(±W/2) = T_amb (торцы в буфере при температуре среды).
    Pinned,
    /// Конвективная граница Робина с разностью первого порядка по внешней нормали.
    Convective,
}

/// Причина отказа стационарного решения.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SteadyProfileError {
    KappaFeedbackRunaway { loop_gain: f64 },
    NonFiniteKappaFeedback { iteration: usize },
    IterationLimit { iterations: usize },
}

impl std::fmt::Display for SteadyProfileError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::KappaFeedbackRunaway { loop_gain } => write!(
                f,
                "стационар κ(T) не существует: коэффициент петли {loop_gain:.6} >= 1"
            ),
            Self::NonFiniteKappaFeedback { iteration } => write!(
                f,
                "κ(T)-итерация дала нечисловое значение на шаге {iteration}"
            ),
            Self::IterationLimit { iterations } => {
                write!(f, "κ(T)-итерация не сошлась за {iterations} шагов")
            }
        }
    }
}

impl std::error::Error for SteadyProfileError {}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoolingModel {
    /// Одинаковая корреляция свободной конвекции на обеих сторонах воздушной кассеты.
    FreeConvection { c_coef: f64 },
    /// Погружённый гель: сверху замкнутый слой буфера, снизу поддон.
    ImmersedHorizontal {
        buffer_depth_m: f64,
        tray_thickness_m: f64,
        tray_k_w_m_k: f64,
    },
    /// Конвекция над верхней стороной геля и под поддоном. Для горячей горизонтальной пластины
    /// McAdams даёт Nu = C·Ra^(1/4) с разными C для восходящего и нисходящего теплоотвода (Incropera, гл. 9).
    OpenHorizontal {
        c_up: f64,
        c_down: f64,
        tray_thickness_m: f64,
        tray_k_w_m_k: f64,
    },
}

/// Коэффициент McAdams для верхней поверхности горячей пластины.
pub const NU_HORIZONTAL_UP: f64 = 0.54;
/// Коэффициент McAdams для нижней поверхности горячей пластины.
pub const NU_HORIZONTAL_DOWN: f64 = 0.27;

/// Тепловая модель и граничные условия описаны в PHYSICS.md §11–14.
#[derive(Debug, Clone, Copy)]
pub struct ThermalConfig {
    pub medium: Medium,
    /// Проводимость буфера при 25 °C, См/м. По умолчанию используется равновесная оценка для TAE.
    pub kappa_s_m: f64,
    pub e_v_per_cm: f64,
    /// Множитель теплоотдачи для анализа чувствительности (номинал 1, диапазон ±30%).
    pub h_scale: f64,
    pub cooling: CoolingModel,
    pub edge: EdgeKind,
    pub width_m: f64,
    pub t_gel_m: f64,
    /// Число узлов поперёк геля.
    pub n_points: usize,
    /// Характерный размер конвекции, м. Для вертикальной пластины это высота, для горизонтальной — площадь/периметр.
    pub l_c_m: f64,
    /// Пересчитывать джоулев источник по проводимости при среднем перегреве геля.
    /// В стационаре положительная обратная связь может исключить физическое решение.
    pub kappa_feedback: bool,
}

impl Default for ThermalConfig {
    fn default() -> Self {
        Self {
            medium: WATER,
            kappa_s_m: crate::buffers::KAPPA_TAE_S_M,
            e_v_per_cm: 5.0,
            h_scale: 1.0,
            cooling: CoolingModel::ImmersedHorizontal {
                buffer_depth_m: DEFAULT_BUFFER_OVERLAY_M,
                tray_thickness_m: DEFAULT_TRAY_THICKNESS_M,
                tray_k_w_m_k: DEFAULT_TRAY_K_W_M_K,
            },
            edge: EdgeKind::Convective,
            width_m: 0.10,
            t_gel_m: 0.006,
            n_points: 801,
            l_c_m: 0.08,
            kappa_feedback: false,
        }
    }
}

impl ThermalConfig {
    /// Джоулев источник, Вт/м³: q_v = κ·E² (E переведена из В/см в В/м).
    pub fn q_v(&self) -> f64 {
        self.kappa_s_m * (self.e_v_per_cm * 100.0).powi(2)
    }

    /// Коэффициент стока через плоскости геля, Вт/(м³·К): s = 2h/t_g.
    pub fn sink_coeff(&self, h_mean: f64) -> f64 {
        2.0 * h_mean / self.t_gel_m
    }

    /// Замыкает корреляцию теплоотдачи по оценке среднего перегрева q_v·t_g/(2h).
    pub fn solve_h_selfconsistent(&self) -> f64 {
        self.solve_h_for_q_v(self.q_v())
    }

    /// Использует фактическую плотность тепловыделения в геле с учётом шунта и профиля проводимости.
    pub fn solve_h_for_q_v(&self, q_v: f64) -> f64 {
        assert!(q_v.is_finite() && q_v >= 0.0);
        bisect(
            |h| {
                let delta_t = q_v * self.t_gel_m / (2.0 * h);
                self.mean_surface_h(delta_t) * self.h_scale - h
            },
            1e-3,
            1e5,
        )
        // Вне диапазона бисекции используем её верхнюю границу.
        .unwrap_or(1.0e5)
    }

    /// Возвращает h сверху и снизу (Вт/(м²·К)), число Рэлея и число Нуссельта слоя. Для воздуха оба числа слоя равны NaN.
    pub fn surface_heat_transfer(&self, delta_t: f64) -> (f64, f64, f64, f64) {
        match self.cooling {
            CoolingModel::FreeConvection { c_coef } => {
                let h = h_natural_convection(&self.medium, delta_t, self.l_c_m, c_coef);
                (h, h, f64::NAN, f64::NAN)
            }
            CoolingModel::OpenHorizontal {
                c_up,
                c_down,
                tray_thickness_m,
                tray_k_w_m_k,
            } => {
                let up = h_natural_convection(&self.medium, delta_t, self.l_c_m, c_up);
                let bare_down = h_natural_convection(&self.medium, delta_t, self.l_c_m, c_down);
                let down = tray_series_heat_transfer(bare_down, tray_thickness_m, tray_k_w_m_k);
                (up, down, f64::NAN, f64::NAN)
            }
            CoolingModel::ImmersedHorizontal {
                buffer_depth_m,
                tray_thickness_m,
                tray_k_w_m_k,
            } => {
                assert!(buffer_depth_m.is_finite() && buffer_depth_m > 0.0);
                assert!(tray_thickness_m.is_finite() && tray_thickness_m > 0.0);
                assert!(tray_k_w_m_k.is_finite() && tray_k_w_m_k > 0.0);
                let ra = rayleigh_horizontal_layer(&self.medium, delta_t, buffer_depth_m);
                let nu = nusselt_horizontal_layer(ra);
                let h_top = nu * self.medium.k / buffer_depth_m;
                // Температура внешней стороны поддона совпадает с температурой перемешанной ванны.
                let h_bottom =
                    tray_series_heat_transfer(f64::INFINITY, tray_thickness_m, tray_k_w_m_k);
                (h_top, h_bottom, ra, nu)
            }
        }
    }

    /// Среднее h сохраняет суммарный поток через обе поверхности.
    pub fn mean_surface_h(&self, delta_t: f64) -> f64 {
        let (h_top, h_bottom, _, _) = self.surface_heat_transfer(delta_t);
        0.5 * (h_top + h_bottom)
    }
}
