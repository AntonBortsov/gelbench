use super::SimState;
use crate::{mobility, mobility2, stain, Params};

impl SimState {
    /// Return the model velocity (cm/s) and diffusivity (cm²/s) at a specified field and temperature. Manual overrides do not enter this calculation.
    pub fn kinetics_at(&self, pop: usize, e_v_cm: f64, temp_k: f64) -> (f64, f64) {
        assert!(pop < self.n_pops_total(), "популяции {pop} нет");
        Self::kinetics_for(
            &self.params,
            self.stain.as_ref().map(|s| s.mode),
            pop,
            e_v_cm,
            temp_k,
        )
    }

    /// Override one population's velocity (cm/s) and diffusivity (cm²/s). Chamber feedback scales the manual reference velocity with the field.
    pub fn set_population_kinetics(&mut self, pop: usize, v_cm_s: f64, d_cm2_s: f64) {
        self.apply_pending_advection();
        self.velocities[pop] = v_cm_s;
        self.diffusivities[pop] = d_cm2_s;
        self.manual_kinetics[pop] = true;
        self.manual_v_ref[pop] = v_cm_s;
        self.view_valid = false;
    }

    /// Mobility uses the field magnitude. Its sign sets drift direction.
    pub(super) fn kinetics_for(
        params: &Params,
        stain_mode: Option<stain::StainMode>,
        pop: usize,
        e_v_cm: f64,
        temp_k: f64,
    ) -> (f64, f64) {
        let n_dna = params.population_sizes_bp.len();
        let theta = (1.0 + mobility::KT_PER_CELSIUS * (temp_k - mobility::T_REF_K)).max(0.0);
        let (mut v, d);
        if pop < n_dna {
            let size_bp = params.population_sizes_bp[pop];
            let form = params
                .population_forms
                .get(pop)
                .copied()
                .unwrap_or_default();
            let c_etbr = match stain_mode {
                Some(stain::StainMode::InGel { c0_ug_ml }) => stain::ug_ml_to_molar(c0_ug_ml),
                _ => 0.0,
            };
            let e_mag = e_v_cm.abs();
            let ratio = mobility2::MobilityV2::form_mobility_ratio(
                size_bp,
                form,
                params.gel.percent,
                e_mag,
                c_etbr,
            );
            let friction = stain::friction_factor(stain::nu_bound(c_etbr));
            v = mobility::MU0_CM2_V_S * ratio * friction * theta * e_v_cm * params.mobility_scale;
            d = mobility2::MobilityV2::dna_diffusivity_cm2_s(size_bp, params.gel.percent, e_mag)
                * theta
                * (temp_k / mobility::T_REF_K);
        } else {
            let dye = params.dyes[pop - n_dna];
            v = dye.mobility_cm2_v_s(temp_k) * e_v_cm * params.mobility_scale;
            d = dye.diffusivity_cm2_s() * theta * (temp_k / mobility::T_REF_K);
        }
        if params.gel.mu_eeo_cm2_v_s != 0.0 {
            v -= params.gel.mu_eeo_cm2_v_s * theta * e_v_cm;
        }
        (v, d)
    }

    pub(super) fn apply_field_kinetics_at(&mut self, e_v_cm: f64, temp_k: f64) {
        if e_v_cm == self.kinetics_e_v_cm && temp_k == self.kinetics_temp_k {
            return;
        }
        let nominal = self.params.field_v_per_cm as f64;
        let mode = self.stain.as_ref().map(|s| s.mode);
        for pop in 0..self.n_pops_total() {
            if self.manual_kinetics[pop] {
                let scale = if nominal != 0.0 {
                    e_v_cm / nominal
                } else {
                    0.0
                };
                self.velocities[pop] = self.manual_v_ref[pop] * scale;
            } else {
                let (v, d) = Self::kinetics_for(&self.params, mode, pop, e_v_cm, temp_k);
                self.velocities[pop] = v;
                self.diffusivities[pop] = d;
            }
        }
        self.kinetics_e_v_cm = e_v_cm;
        self.kinetics_temp_k = temp_k;
        self.view_valid = false;
    }
}
