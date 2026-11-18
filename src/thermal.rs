//! Joule heating of a submerged or exposed gel with transverse cooling and a mixed bath.

mod bath;
mod config;
mod heat_transfer;
mod lane_temperature;
mod steady;
#[cfg(test)]
mod tests;
mod transient;

pub use bath::{BathThermalConfig, BathThermalState, ThermalStatus};
pub use config::{
    CoolingModel, EdgeKind, SteadyProfileError, ThermalConfig, NU_HORIZONTAL_DOWN, NU_HORIZONTAL_UP,
};
pub use heat_transfer::{
    bisect, h_natural_convection, nusselt_horizontal_layer, rayleigh_horizontal_layer, Medium, AIR,
    WATER,
};
pub use lane_temperature::{
    local_delta_t_at, reference_band_path_cm, theta_from_delta_t, theta_profile,
};

/// Gravitational acceleration (m/s²).
pub const G: f64 = 9.81;

/// Thermal conductivity of water-rich agarose (W/m/K).
pub const K_GEL: f64 = 0.60;

/// Volumetric heat capacity of the gel (J/m³/K).
pub const RHO_CP: f64 = 4.2e6;

/// Target number of nodes per characteristic thermal length.
pub const MIN_NODES_PER_THERMAL_LENGTH: f64 = 20.0;

/// Volumetric heat capacity of water-like buffer (J/L/K).
pub const BUFFER_HEAT_CAPACITY_J_L_K: f64 = 4_180.0;

/// Warning threshold (K).
pub const OVERHEAT_WARNING_K: f64 = 318.15;
/// Lower temperature bound of the thermal property model (K).
pub const THERMAL_MODEL_MIN_K: f64 = 273.15;
/// Upper temperature bound of the thermal property model (K).
pub const THERMAL_MODEL_MAX_K: f64 = 373.15;

pub const DEFAULT_TRAY_THICKNESS_M: f64 = 0.003;
pub const DEFAULT_TRAY_K_W_M_K: f64 = 0.20;
pub const DEFAULT_BUFFER_OVERLAY_M: f64 = 0.003;
