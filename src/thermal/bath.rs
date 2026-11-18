//! Bath enthalpy and transverse gel heating. Temperatures feed the next transport step.

use super::{
    local_delta_t_at, theta_from_delta_t, CoolingModel, EdgeKind, ThermalConfig, AIR,
    BUFFER_HEAT_CAPACITY_J_L_K, DEFAULT_BUFFER_OVERLAY_M, DEFAULT_TRAY_K_W_M_K,
    DEFAULT_TRAY_THICKNESS_M, K_GEL, MIN_NODES_PER_THERMAL_LENGTH, NU_HORIZONTAL_DOWN,
    NU_HORIZONTAL_UP, OVERHEAT_WARNING_K, RHO_CP, THERMAL_MODEL_MAX_K, THERMAL_MODEL_MIN_K, WATER,
};
use crate::fmath::DetMath;

/// Объём буфера здесь относится только к жидкости. Теплоёмкость геля задаёт геометрия.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BathThermalConfig {
    pub ambient_k: f64,
    pub buffer_volume_l: f64,
    /// Суммарная теплопередача всей камеры наружу, Вт/К.
    pub ua_w_per_k: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum ThermalStatus {
    Normal = 0,
    Overheat = 1,
    MeltRisk = 2,
    OutOfDomain = 3,
}

/// Полная мощность камеры входит в баланс ванны. Доля мощности в геле задаёт поперечный профиль.
#[derive(Debug, Clone)]
pub struct BathThermalState {
    pub config: BathThermalConfig,
    pub initial_temp_k: f64,
    pub bath_temp_k: f64,
    pub heat_capacity_j_per_k: f64,
    /// В сосредоточенной модели гель и жидкость имеют одну среднюю температуру.
    gel_heat_capacity_j_per_k: f64,
    pub energy_in_j: f64,
    /// Энтальпия долитого или слитого буфера относительно начальной температуры ванны.
    pub buffer_exchange_enthalpy_j: f64,
    /// Интеграл теплопотерь в среду. Отрицателен, когда среда греет ванну.
    pub energy_to_ambient_j: f64,
    pub last_total_power_w: f64,
    pub last_gel_power_w: f64,
    pub last_outside_gel_power_w: f64,
    local_delta_t_k: Vec<f64>,
    local_cfg: ThermalConfig,
    subcolumns: usize,
}

impl BathThermalState {
    /// Для минимальной тепловой длины ℓ = √(k·t / 2h) берётся верхняя оценка водяного теплоотвода.
    fn profile_points(width_m: f64, t_gel_m: f64) -> usize {
        const H_UPPER_BOUND_W_M2K: f64 = 400.0;
        let ell_m = (K_GEL * t_gel_m / (2.0 * H_UPPER_BOUND_W_M2K)).sqrt();
        let want = (MIN_NODES_PER_THERMAL_LENGTH * width_m / ell_m).ceil() as usize;
        (want | 1).clamp(129, 1601)
    }

    pub fn new(
        config: BathThermalConfig,
        initial_temp_k: f64,
        geometry: crate::circuit::CellGeometry,
        subcolumns: usize,
    ) -> Self {
        assert!(config.ambient_k.is_finite() && config.ambient_k > 0.0);
        assert!(config.buffer_volume_l.is_finite() && config.buffer_volume_l > 0.0);
        assert!(config.ua_w_per_k.is_finite() && config.ua_w_per_k >= 0.0);
        assert!(initial_temp_k.is_finite() && initial_temp_k > 0.0);
        let gel_volume_m3 = geometry.l_x_m * geometry.width_m * geometry.t_gel_m;
        let gel_heat_capacity_j_per_k = RHO_CP * gel_volume_m3;
        let heat_capacity_j_per_k =
            BUFFER_HEAT_CAPACITY_J_L_K * config.buffer_volume_l + gel_heat_capacity_j_per_k;
        // Сетка теплового слоя независима от числа полос. Редкая сетка завышает перепад центр–край.
        let n_points = Self::profile_points(geometry.width_m, geometry.t_gel_m);
        let local_cfg = ThermalConfig {
            medium: WATER,
            kappa_s_m: 1.0,
            e_v_per_cm: 0.0,
            h_scale: 1.0,
            cooling: CoolingModel::ImmersedHorizontal {
                buffer_depth_m: DEFAULT_BUFFER_OVERLAY_M,
                tray_thickness_m: DEFAULT_TRAY_THICKNESS_M,
                tray_k_w_m_k: DEFAULT_TRAY_K_W_M_K,
            },
            edge: EdgeKind::Convective,
            width_m: geometry.width_m,
            t_gel_m: geometry.t_gel_m,
            n_points,
            l_c_m: geometry.l_x_m,
            kappa_feedback: false,
        };
        Self {
            config,
            initial_temp_k,
            bath_temp_k: initial_temp_k,
            heat_capacity_j_per_k,
            gel_heat_capacity_j_per_k,
            energy_in_j: 0.0,
            buffer_exchange_enthalpy_j: 0.0,
            energy_to_ambient_j: 0.0,
            last_total_power_w: 0.0,
            last_gel_power_w: 0.0,
            last_outside_gel_power_w: 0.0,
            local_delta_t_k: vec![0.0; n_points],
            local_cfg,
            subcolumns: subcolumns.max(1),
        }
    }

    /// Толщина слоя над гелем определяется уровнем буфера и геометрией поддона.
    pub fn new_in_chamber(
        config: BathThermalConfig,
        initial_temp_k: f64,
        geometry: crate::circuit::CellGeometry,
        subcolumns: usize,
        tray_surface_m: f64,
        level_m: f64,
    ) -> Self {
        let mut state = Self::new(config, initial_temp_k, geometry, subcolumns);
        let gel_top_m = tray_surface_m + geometry.t_gel_m;
        state.set_immersion_depth_m((level_m - gel_top_m).max(0.0));
        state
    }

    /// Покрытый гель отдаёт тепло через слой буфера и поддон, оголённый — в воздух и через поддон.
    pub fn set_immersion_depth_m(&mut self, depth_m: f64) {
        assert!(depth_m.is_finite() && depth_m >= 0.0);
        if depth_m > 0.0 {
            self.local_cfg.medium = WATER;
            self.local_cfg.cooling = CoolingModel::ImmersedHorizontal {
                buffer_depth_m: depth_m,
                tray_thickness_m: DEFAULT_TRAY_THICKNESS_M,
                tray_k_w_m_k: DEFAULT_TRAY_K_W_M_K,
            };
        } else {
            self.local_cfg.medium = AIR;
            self.local_cfg.cooling = CoolingModel::OpenHorizontal {
                c_up: NU_HORIZONTAL_UP,
                c_down: NU_HORIZONTAL_DOWN,
                tray_thickness_m: DEFAULT_TRAY_THICKNESS_M,
                tray_k_w_m_k: DEFAULT_TRAY_K_W_M_K,
            };
        }
    }

    pub fn gel_heat_capacity_j_per_k(&self) -> f64 {
        self.gel_heat_capacity_j_per_k
    }

    pub fn buffer_heat_capacity_j_per_k(&self) -> f64 {
        BUFFER_HEAT_CAPACITY_J_L_K * self.config.buffer_volume_l
    }

    /// Долив мгновенно смешивается с буфером и гелем. Геометрический уровень задаётся отдельно.
    pub fn add_buffer(&mut self, volume_l: f64, temperature_k: f64) {
        assert!(volume_l.is_finite() && volume_l >= 0.0);
        assert!(temperature_k.is_finite() && temperature_k > 0.0);
        if volume_l == 0.0 {
            return;
        }
        let added_capacity = BUFFER_HEAT_CAPACITY_J_L_K * volume_l;
        let old_capacity = self.heat_capacity_j_per_k;
        self.bath_temp_k = (old_capacity * self.bath_temp_k + added_capacity * temperature_k)
            / (old_capacity + added_capacity);
        self.heat_capacity_j_per_k = old_capacity + added_capacity;
        self.config.buffer_volume_l += volume_l;
        self.buffer_exchange_enthalpy_j += added_capacity * (temperature_k - self.initial_temp_k);
    }

    /// Слив уносит жидкость при температуре ванны.
    pub fn drain_buffer(&mut self, volume_l: f64) {
        assert!(volume_l.is_finite() && volume_l >= 0.0);
        assert!(volume_l <= self.config.buffer_volume_l);
        if volume_l == 0.0 {
            return;
        }
        let removed_capacity = BUFFER_HEAT_CAPACITY_J_L_K * volume_l;
        self.buffer_exchange_enthalpy_j -=
            removed_capacity * (self.bath_temp_k - self.initial_temp_k);
        self.heat_capacity_j_per_k -= removed_capacity;
        self.config.buffer_volume_l -= volume_l;
    }

    /// Аналитический баланс за шаг при постоянной электрической мощности.
    pub fn advance(&mut self, dt_s: f64, total_power_w: f64, gel_power_w: f64) {
        assert!(dt_s.is_finite() && dt_s > 0.0);
        assert!(total_power_w.is_finite() && total_power_w >= 0.0);
        assert!(gel_power_w.is_finite() && gel_power_w >= 0.0);
        let gel_power_w = gel_power_w.min(total_power_w);
        let old_t = self.bath_temp_k;
        let c = self.heat_capacity_j_per_k;
        let ua = self.config.ua_w_per_k;
        let new_t = if ua == 0.0 {
            old_t + total_power_w * dt_s / c
        } else {
            let equilibrium_k = self.config.ambient_k + total_power_w / ua;
            equilibrium_k + (old_t - equilibrium_k) * (-ua * dt_s / c).dexp()
        };
        self.bath_temp_k = if new_t.is_finite() { new_t } else { old_t };
        let input_step_j = total_power_w * dt_s;
        let stored_step_j = c * (self.bath_temp_k - old_t);
        self.energy_in_j += input_step_j;
        // Теплоотдачу находим из баланса, чтобы погрешность отдельной квадратуры не накапливалась.
        self.energy_to_ambient_j += input_step_j - stored_step_j;
        self.last_total_power_w = total_power_w;
        self.last_gel_power_w = gel_power_w;
        self.last_outside_gel_power_w = total_power_w - gel_power_w;

        let gel_volume_m3 = self.local_cfg.width_m * self.local_cfg.t_gel_m * self.local_cfg.l_c_m;
        let q_v = if gel_volume_m3 > 0.0 {
            gel_power_w / gel_volume_m3
        } else {
            0.0
        };
        let h = self.local_cfg.solve_h_for_q_v(q_v);
        self.local_delta_t_k = self
            .local_cfg
            .transient_step(&self.local_delta_t_k, dt_s, h, q_v);
        for dt in &mut self.local_delta_t_k {
            if !dt.is_finite() || *dt < 0.0 {
                *dt = 0.0;
            }
        }
    }

    /// Обратная связь ограничена диапазоном свойств. Энергетический баланс хранит фактическую температуру.
    pub fn feedback_temperature_k(&self) -> f64 {
        self.bath_temp_k
            .clamp(THERMAL_MODEL_MIN_K, THERMAL_MODEL_MAX_K)
    }

    /// Одномерный продольный перенос использует средний поперечный перегрев геля.
    pub fn feedback_mean_gel_temperature_k(&self) -> f64 {
        let mean_local_delta_k =
            self.local_delta_t_k.iter().sum::<f64>() / self.local_delta_t_k.len() as f64;
        (self.bath_temp_k + mean_local_delta_k).clamp(THERMAL_MODEL_MIN_K, THERMAL_MODEL_MAX_K)
    }

    pub fn local_delta_t_k(&self) -> &[f64] {
        &self.local_delta_t_k
    }

    pub fn max_gel_temp_k(&self) -> f64 {
        self.bath_temp_k + self.local_delta_t_k.iter().copied().fold(0.0, f64::max)
    }

    pub fn stored_energy_j(&self) -> f64 {
        self.heat_capacity_j_per_k * (self.bath_temp_k - self.initial_temp_k)
    }

    pub fn balance_residual_j(&self) -> f64 {
        self.energy_in_j + self.buffer_exchange_enthalpy_j
            - self.energy_to_ambient_j
            - self.stored_energy_j()
    }

    pub fn theta_for_subcolumns(&self) -> Vec<f64> {
        let local = theta_from_delta_t(&self.local_cfg, &self.local_delta_t_k, self.subcolumns);
        self.compose_local_with_bath(local)
    }

    /// Координаты дорожек заданы в метрах от середины геля. За пределами ширины берётся температура края.
    pub fn theta_for_positions(&self, y_m: &[f64]) -> Vec<f64> {
        let kt = crate::mobility::KT_PER_CELSIUS;
        let local: Vec<f64> = y_m
            .iter()
            .map(|&y| 1.0 + kt * local_delta_t_at(&self.local_cfg, &self.local_delta_t_k, y))
            .collect();
        self.compose_local_with_bath(local)
    }

    fn compose_local_with_bath(&self, local: Vec<f64>) -> Vec<f64> {
        let kt = crate::mobility::KT_PER_CELSIUS;
        let bath_k = self.feedback_temperature_k();
        let bath_theta = (1.0 + kt * (bath_k - crate::mobility::T_REF_K)).max(f64::MIN_POSITIVE);
        local
            .into_iter()
            .map(|local_theta| {
                let local_delta_k = (local_theta - 1.0) / kt;
                let total_k =
                    (bath_k + local_delta_k).clamp(THERMAL_MODEL_MIN_K, THERMAL_MODEL_MAX_K);
                let total_theta = 1.0 + kt * (total_k - crate::mobility::T_REF_K);
                // Кинетика уже включает поправку ванны. Деление исключает двойной учёт её температуры.
                total_theta / bath_theta
            })
            .collect()
    }

    pub fn status(&self) -> ThermalStatus {
        let max_t = self.max_gel_temp_k();
        if self.bath_temp_k < THERMAL_MODEL_MIN_K || max_t >= THERMAL_MODEL_MAX_K {
            ThermalStatus::OutOfDomain
        } else if max_t >= crate::chamber::MELT_RISK_K {
            ThermalStatus::MeltRisk
        } else if max_t >= OVERHEAT_WARNING_K {
            ThermalStatus::Overheat
        } else {
            ThermalStatus::Normal
        }
    }
}
