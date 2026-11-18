use super::SimState;
use crate::{chamber, thermal};

/// Rejected liquid transfer. Neither volume nor level is changed on failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BufferActionError {
    ChamberNotEnabled,
    DynamicHeatingNotEnabled,
    InvalidVolume,
    InvalidTemperature,
    InvalidLevel,
}

impl std::fmt::Display for BufferActionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let text = match self {
            Self::ChamberNotEnabled => "камера не включена",
            Self::DynamicHeatingNotEnabled => "динамический нагрев не включён",
            Self::InvalidVolume => "недопустимый объём буфера",
            Self::InvalidTemperature => "температура долива должна быть в [250; 330] K",
            Self::InvalidLevel => "уровень должен быть в [0; 0.2] м",
        };
        f.write_str(text)
    }
}

impl std::error::Error for BufferActionError {}
impl SimState {
    /// Install the chamber and use its gel field for both drift kinetics and electrolyte transport.
    pub fn enable_chamber(&mut self, spec: chamber::ChamberSpec) {
        assert!(
            self.bath_thermal.is_none(),
            "enable_chamber: сначала отключите dynamic heating"
        );
        let ch = chamber::ChamberState::new(spec);
        let e0 = ch.last.e_gel_v_per_cm;
        self.chamber = Some(ch);
        self.apply_field_kinetics_at(e0, self.params.gel.temp_k);
        if let Some(m) = self.medium.as_mut() {
            m.set_u_gel(e0 * self.params.gel_length_cm as f64);
        }
        self.view_valid = false;
    }

    /// Attach a thermal bath using the chamber geometry and current buffer level.
    pub fn enable_dynamic_heating(&mut self, config: thermal::BathThermalConfig) {
        let ch = self
            .chamber
            .as_ref()
            .expect("enable_dynamic_heating: сначала enable_chamber");
        let state = thermal::BathThermalState::new_in_chamber(
            config,
            self.params.gel.temp_k,
            ch.spec.geometry,
            self.params.topology.subcolumns,
            ch.spec.tray_surface_m,
            ch.level_m,
        );
        self.dynamic_subcolumn_theta = Some(state.theta_for_subcolumns());
        self.dynamic_lane_theta = self.lane_y_m.as_ref().map(|y| state.theta_for_positions(y));
        self.bath_thermal = Some(state);
        self.view_valid = false;
    }

    /// Return to the specified isothermal gel temperature and discard thermal multipliers.
    pub fn disable_dynamic_heating(&mut self) {
        self.apply_pending_advection();
        self.bath_thermal = None;
        self.dynamic_subcolumn_theta = None;
        self.dynamic_lane_theta = None;
        self.apply_field_kinetics_at(self.kinetics_e_v_cm, self.params.gel.temp_k);
        if let Some(m) = self.medium.as_mut() {
            m.set_temperature_k(self.params.gel.temp_k);
        }
        self.view_valid = false;
    }

    pub fn thermal_bath(&self) -> Option<&thermal::BathThermalState> {
        self.bath_thermal.as_ref()
    }

    /// Changing the liquid level also changes the immersed depth used for heat transfer.
    pub fn set_buffer_level(&mut self, level_m: f64) {
        let ch = self
            .chamber
            .as_mut()
            .expect("set_buffer_level: сначала enable_chamber");
        ch.level_m = level_m;
        let gel_top_m = ch.spec.tray_surface_m + ch.spec.geometry.t_gel_m;
        if let Some(thermal) = self.bath_thermal.as_mut() {
            thermal.set_immersion_depth_m((level_m - gel_top_m).max(0.0));
        }
    }

    /// Mix the added volume (L) into the thermal bath at the specified temperature (K).
    pub fn add_buffer(&mut self, volume_l: f64, temperature_k: f64) {
        self.bath_thermal
            .as_mut()
            .expect("add_buffer: сначала enable_dynamic_heating")
            .add_buffer(volume_l, temperature_k);
        self.view_valid = false;
    }

    pub fn drain_buffer(&mut self, volume_l: f64) {
        self.bath_thermal
            .as_mut()
            .expect("drain_buffer: сначала enable_dynamic_heating")
            .drain_buffer(volume_l);
        self.view_valid = false;
    }

    /// The bath and chamber level are updated together after validation.
    pub fn refill_buffer(
        &mut self,
        volume_l: f64,
        temperature_k: f64,
        level_m: f64,
    ) -> Result<(), BufferActionError> {
        if self.chamber.is_none() {
            return Err(BufferActionError::ChamberNotEnabled);
        }
        let thermal = self
            .bath_thermal
            .as_ref()
            .ok_or(BufferActionError::DynamicHeatingNotEnabled)?;
        if !(0.0..=20.0).contains(&volume_l) || thermal.config.buffer_volume_l + volume_l > 20.0 {
            return Err(BufferActionError::InvalidVolume);
        }
        if !(250.0..=330.0).contains(&temperature_k) {
            return Err(BufferActionError::InvalidTemperature);
        }
        if !(0.0..=0.2).contains(&level_m) {
            return Err(BufferActionError::InvalidLevel);
        }
        self.add_buffer(volume_l, temperature_k);
        self.set_buffer_level(level_m);
        Ok(())
    }

    pub fn drain_buffer_to_level(
        &mut self,
        volume_l: f64,
        level_m: f64,
    ) -> Result<(), BufferActionError> {
        if self.chamber.is_none() {
            return Err(BufferActionError::ChamberNotEnabled);
        }
        let thermal = self
            .bath_thermal
            .as_ref()
            .ok_or(BufferActionError::DynamicHeatingNotEnabled)?;
        if !(0.0..=thermal.config.buffer_volume_l).contains(&volume_l) {
            return Err(BufferActionError::InvalidVolume);
        }
        if !(0.0..=0.2).contains(&level_m) {
            return Err(BufferActionError::InvalidLevel);
        }
        self.drain_buffer(volume_l);
        self.set_buffer_level(level_m);
        Ok(())
    }

    pub fn chamber(&self) -> Option<&chamber::ChamberState> {
        self.chamber.as_ref()
    }
}
