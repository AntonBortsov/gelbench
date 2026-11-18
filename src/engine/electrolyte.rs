use super::SimState;
use crate::{medium, medium_loop};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SimFault {
    Medium(medium_loop::MediumFault),
}

impl std::fmt::Display for SimFault {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Medium(m) => write!(f, "{m}"),
        }
    }
}

impl std::error::Error for SimFault {}

const FIELD_UNIFORMITY_TOL: f64 = 1e-3;
impl SimState {
    /// Install a positive, dimensionless multiplier of the nominal field at every cell.
    pub fn set_field_phi(&mut self, phi: &[f64]) {
        let n = self.params.topology.cells_per_lane.max(1);
        assert!(
            phi.len() == n,
            "phi.len() = {} против cells_per_lane = {n}",
            phi.len()
        );
        assert!(
            phi.iter().all(|p| p.is_finite() && *p > 0.0),
            "φ(x) обязана быть конечной и положительной"
        );
        self.apply_pending_advection();
        self.field_phi = Some(phi.to_vec());
        self.view_valid = false;
    }

    /// Initialize every lane with the selected electrolyte. The gel voltage is the nominal field times the gel length until a chamber supplies its operating point.
    pub fn enable_medium(&mut self, spec: medium::MediumSpec, ionic: medium::IonicStrengthModel) {
        assert!(
            self.casting_theta_lane.is_none(),
            "петля среды совместно с №25 v1 не поддерживает (PHYSICS §32)"
        );
        let topo = self.params.topology;
        let n = topo.cells_per_lane.max(1);
        let dx_m = self.dx * 1.0e-2;
        let u_gel_v = self.params.field_v_per_cm as f64 * self.params.gel_length_cm as f64;
        self.medium = Some(medium_loop::MediumLoop::new(
            spec,
            topo.lanes,
            n,
            dx_m,
            u_gel_v,
            self.params.gel.temp_k,
            ionic,
        ));
        self.view_valid = false;
    }

    pub fn well_background(&self) -> Option<medium_loop::WellLoad> {
        self.medium
            .as_ref()
            .map(|m| medium_loop::WellLoad::background(m.spec()))
    }

    /// Set the solution at a well. Use `try_load_well` to handle a composition outside the buffer capacity.
    pub fn load_well(&mut self, lane: usize, load: &medium_loop::WellLoad) {
        self.try_load_well(lane, load)
            .expect("load_well: состав лунки вне буферной ёмкости — используйте try_load_well");
    }

    /// Attempt to replace the solution near the well center. Concentrations outside the speciation range return a medium fault.
    pub fn try_load_well(
        &mut self,
        lane: usize,
        load: &medium_loop::WellLoad,
    ) -> Result<(), SimFault> {
        let center_m = self.params.well_position_cm as f64 * 1.0e-2;
        let m = self
            .medium
            .as_mut()
            .expect("load_well: сначала enable_medium");
        m.try_load_well(lane, center_m, 0.5e-3, 0.3e-3, load)
            .map_err(SimFault::Medium)?;
        self.view_valid = false;
        Ok(())
    }

    /// Conductivity of freshly prepared buffer at 25 °C, in S/m.
    pub fn medium_fresh_kappa_s_m(&self) -> Option<f64> {
        self.medium.as_ref().map(|m| m.fresh_kappa_s_m())
    }

    pub fn medium_diag(&self) -> Option<&medium_loop::MediumDiag> {
        self.medium.as_ref().map(|m| &m.diag)
    }

    pub fn medium_temperature_k(&self) -> Option<f64> {
        self.medium.as_ref().map(|m| m.temperature_k())
    }

    /// Effective gel resistivity inferred from lane conductivity profiles, in Ω·m.
    pub fn medium_gel_resistivity_ohm_m(&self) -> Option<f64> {
        self.medium.as_ref().map(|m| m.gel_resistivity_ohm_m())
    }

    pub fn enable_reservoir_books(&mut self, volume_l: f64, lane_area_m2: f64) {
        self.medium
            .as_mut()
            .expect("enable_reservoir_books: сначала enable_medium")
            .enable_books(volume_l, lane_area_m2);
    }

    pub fn reservoir_books(&self) -> Option<&medium_loop::LoopBooks> {
        self.medium.as_ref().and_then(|m| m.books())
    }

    /// Concentration profile for dissolved family `k` along a lane.
    pub fn medium_family(&self, lane: usize, k: usize) -> Option<&[f64]> {
        self.medium.as_ref().map(|m| m.family(lane, k))
    }

    pub fn medium_ph(&self, lane: usize) -> Option<&[f64]> {
        self.medium.as_ref().map(|m| m.ph_profile(lane))
    }

    /// Local-to-nominal field ratio from the electrolyte solve.
    pub fn lane_field_phi(&self, lane: usize) -> Option<&[f64]> {
        self.medium.as_ref().map(|m| m.phi_lane(lane))
    }
}

/// The electrolyte field takes precedence over an externally supplied profile.
pub(super) fn field_for_lane<'a>(
    medium: Option<&'a medium_loop::MediumLoop>,
    shared: Option<&'a [f64]>,
    lane: usize,
) -> Option<&'a [f64]> {
    medium.map(|m| m.phi_lane(lane)).or(shared)
}

pub(super) fn field_is_nonuniform(profile: &[f64]) -> bool {
    profile
        .iter()
        .any(|&p| (p - 1.0).abs() > FIELD_UNIFORMITY_TOL)
}
