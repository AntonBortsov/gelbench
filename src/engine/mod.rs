mod advection;
mod chamber_control;
mod electrolyte;
mod kinetics;
mod measurements;
mod sample;
mod spatial;
mod staining;
#[cfg(test)]
mod tests;
mod timestep;

pub use chamber_control::BufferActionError;
pub use electrolyte::SimFault;
pub use measurements::Observation;
pub use spatial::select_subcolumn_count;

use crate::{chamber, medium_loop, stain, thermal, MobilityModel, Params};

/// State of one gel run.
pub struct SimState {
    params: Params,
    model: MobilityModel,
    velocities: Vec<f64>,
    subcolumn_theta: Vec<f64>,
    casting_theta_lane: Option<Vec<f64>>,
    diffusivities: Vec<f64>,
    /// State storage is contiguous in lane, population, subcolumn, cell order.
    conc: Vec<f64>,
    /// Reused by both advection and diffusion.
    adv: Vec<f64>,
    /// Observation grid after fractional shifts. Transport advances the unshifted rows.
    view_conc: Vec<f64>,
    view_valid: bool,
    faces: Vec<f64>,
    /// Uncommitted signed displacement (cells), one value per row.
    shift_cells: Vec<f64>,
    /// Shared external field divided by nominal field.
    field_phi: Option<Vec<f64>>,
    medium: Option<medium_loop::MediumLoop>,
    stain: Option<stain::StainState>,
    chamber: Option<chamber::ChamberState>,
    bath_thermal: Option<thermal::BathThermalState>,
    dynamic_subcolumn_theta: Option<Vec<f64>>,
    lane_y_m: Option<Vec<f64>>,
    dynamic_lane_theta: Option<Vec<f64>>,
    fault: Option<SimFault>,
    /// Anode outflow per lane and population, summed over subcolumns.
    exited: Vec<f64>,
    /// Additional outflow from fractional shifts in the observation copy.
    view_exited: Vec<f64>,
    manual_kinetics: Vec<bool>,
    /// Reference velocities at the nominal field for manual overrides (cm/s).
    manual_v_ref: Vec<f64>,
    kinetics_e_v_cm: f64,
    kinetics_temp_k: f64,
    /// Cell edges (cm), including both gel boundaries.
    xe: Vec<f64>,
    /// Cell centers (cm) used to interpolate the local field.
    xc: Vec<f64>,
    well_widened_to_cells: Option<usize>,
    /// Longitudinal cell spacing (cm).
    dx: f64,
    t: f64,
    steps: u64,
}
impl SimState {
    /// Raised-cosine wells start at unit peak height. Loading sets mass in ng. Wells narrower than three cells are widened.
    pub fn new(params: Params) -> Self {
        let topo = params.topology;
        let n_dna = params.population_sizes_bp.len();
        let n_pops = n_dna + params.dyes.len();
        assert!(
            n_pops == topo.populations,
            "population_sizes_bp.len() + dyes.len() = {n_pops} не совпадает с \
             topology.populations = {}",
            topo.populations
        );
        let n = topo.cells_per_lane.max(1);
        let total_cells = topo.total_cells().max(1);
        let dx = params.gel_length_cm as f64 / n as f64;

        let peak =
            ((params.well_position_cm as f64 / dx).floor() as usize).clamp(2, n.saturating_sub(4));
        const WELL_HALF_WIDTH_CM: f64 = 0.05;
        let well_full_cells = 2.0 * WELL_HALF_WIDTH_CM / dx;
        let well_widened_to_cells = if well_full_cells < 3.0 { Some(3) } else { None };
        let half_w_cells = if well_widened_to_cells.is_some() {
            1.5
        } else {
            WELL_HALF_WIDTH_CM / dx
        };
        let mut conc = vec![0.0; total_cells];
        if topo.lanes > 0 && n_pops > 0 && topo.subcolumns > 0 {
            let (first, others) = conc.split_at_mut(n);
            for (x, c) in first.iter_mut().enumerate() {
                *c = Self::well_profile_height(x, peak, half_w_cells);
            }
            for row in others.chunks_exact_mut(n) {
                row.copy_from_slice(first);
            }
        }

        let model = MobilityModel::calibrated();
        let field = params.field_v_per_cm as f64;
        assert!(
            params.population_forms.is_empty()
                || params.population_forms.len() == params.population_sizes_bp.len(),
            "population_forms обязан быть пуст или параллелен population_sizes_bp"
        );
        let mut velocities = Vec::with_capacity(n_pops);
        let mut diffusivities = Vec::with_capacity(n_pops);
        for pop in 0..n_pops {
            let (v, d) = Self::kinetics_for(&params, None, pop, field, params.gel.temp_k);
            velocities.push(v);
            diffusivities.push(d);
        }

        Self {
            adv: vec![0.0; total_cells],
            faces: vec![0.0; n + 1],
            shift_cells: vec![0.0; topo.lanes.max(1) * n_pops.max(1) * topo.subcolumns.max(1)],
            field_phi: None,
            medium: None,
            stain: None,
            chamber: None,
            bath_thermal: None,
            dynamic_subcolumn_theta: None,
            lane_y_m: None,
            dynamic_lane_theta: None,
            fault: None,
            exited: vec![0.0; topo.lanes.max(1) * n_pops.max(1)],
            view_exited: vec![0.0; topo.lanes.max(1) * n_pops.max(1)],
            manual_kinetics: vec![false; n_pops.max(1)],
            manual_v_ref: vec![0.0; n_pops.max(1)],
            kinetics_e_v_cm: params.field_v_per_cm as f64,
            kinetics_temp_k: params.gel.temp_k,
            xe: (0..=n).map(|i| i as f64 * dx).collect(),
            xc: (0..n).map(|i| (i as f64 + 0.5) * dx).collect(),
            well_widened_to_cells,
            view_conc: conc.clone(),
            view_valid: true,
            conc,
            dx,
            model,
            velocities,
            subcolumn_theta: vec![1.0; topo.subcolumns],
            casting_theta_lane: None,
            diffusivities,
            params,
            t: 0.0,
            steps: 0,
        }
    }

    pub fn params(&self) -> &Params {
        &self.params
    }

    pub fn time_s(&self) -> f64 {
        self.t
    }

    pub fn steps(&self) -> u64 {
        self.steps
    }

    /// Failure of the last macrostep, if any. Further steps leave the run stopped.
    pub fn fault(&self) -> Option<SimFault> {
        self.fault
    }

    /// DNA populations occupy the first slots, followed by tracking dyes.
    fn n_pops_total(&self) -> usize {
        self.params.population_sizes_bp.len() + self.params.dyes.len()
    }

    /// Reference mobility model for calibration comparisons. Transport uses the field-dependent model in `mobility2`.
    pub fn model(&self) -> &MobilityModel {
        &self.model
    }

    /// Minimum grid width used for an under-resolved starting well, if applicable.
    pub fn well_widened_to_cells(&self) -> Option<usize> {
        self.well_widened_to_cells
    }
}
