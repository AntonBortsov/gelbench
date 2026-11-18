//! Gyration radius and fiber impalement for closed DNA. Intercalation changes contour length
//! and superhelical density. The engine applies the additional dye friction separately.

use super::{
    MobilityV2, K_FREE_CIRC, K_FREE_SC, K_IMPALE, RG_CIRC_RATIO, RG_SC_RATIO, SIGMA_COMPACT,
};
use crate::{fmath::DetMath, materials::agarose::rg_nm};

/// DNA topology. Native supercoiled input uses a negative superhelical density.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize, Default)]
pub enum DnaForm {
    #[default]
    Linear,
    Supercoiled {
        sigma0: f64,
    },
    Relaxed,
}

impl MobilityV2 {
    /// Подбираем длину линейной цепи с тем же Rg и учитываем захват кольца волокнами геля.
    /// Штраф растёт с (Rg/a)² и воспроизводит задержку крупных открытых колец (Mickel 1977).
    fn closed_form_ratio(rg_eff_nm: f64, contour_extension: f64, percent: f64, e_v_cm: f64) -> f64 {
        let a = Self::pore_nm(percent);
        // Контур удлиняется после обращения Rg, чтобы интеркаляция учитывалась один раз.
        let l_eff = Self::bp_at_rg(rg_eff_nm / contour_extension.sqrt()) * contour_extension;
        let base = Self::ratio_core(rg_eff_nm, l_eff, l_eff, percent, e_v_cm);
        let x = rg_eff_nm / a;
        base * (-K_IMPALE * x * x).dexp()
    }

    /// μ/μ₀ при связанной доле красителя ν, вычисленной по свободному EtBr.
    /// У линейной цепи и открытого кольца контур растёт в 1+ν раз. У сверхспирали Rg уже учитывает раскрутку
    /// на 26° при интеркаляции. После полной релаксации положительная сверхспирализация снова уплотняет цепь.
    pub fn form_mobility_ratio(
        n_bp: f64,
        form: DnaForm,
        percent: f64,
        e_v_cm: f64,
        free_etbr_m: f64,
    ) -> f64 {
        let nu = crate::stain::nu_bound(free_etbr_m);
        let extension = crate::stain::contour_length_factor(nu);
        let rg_extension = extension.sqrt();
        match form {
            DnaForm::Linear => Self::retarding_only(
                Self::ratio_core(
                    rg_nm(n_bp) * rg_extension,
                    n_bp * extension,
                    n_bp * extension,
                    percent,
                    e_v_cm,
                ),
                || Self::mobility_ratio(n_bp, percent, e_v_cm),
                extension,
            ),
            DnaForm::Relaxed => Self::retarding_only(
                K_FREE_CIRC
                    * Self::closed_form_ratio(
                        rg_nm(n_bp) * rg_extension * RG_CIRC_RATIO,
                        extension,
                        percent,
                        e_v_cm,
                    ),
                || {
                    K_FREE_CIRC
                        * Self::closed_form_ratio(rg_nm(n_bp) * RG_CIRC_RATIO, 1.0, percent, e_v_cm)
                },
                extension,
            ),
            DnaForm::Supercoiled { sigma0 } => {
                assert!(sigma0 < 0.0, "нативная суперспиральность отрицательна");
                let unwind_per_nu = crate::stain::UNWIND_DEG * crate::stain::BP_PER_TURN / 360.0;
                let sigma_eff = sigma0 + unwind_per_nu * nu;
                let compact = (sigma_eff.abs() / SIGMA_COMPACT).min(1.0);
                let rg_factor = RG_CIRC_RATIO - (RG_CIRC_RATIO - RG_SC_RATIO) * compact;
                K_FREE_SC * Self::closed_form_ratio(rg_nm(n_bp) * rg_factor, 1.0, percent, e_v_cm)
            }
        }
    }

    /// На немонотонном стыке просеивания и рептации удлинение контура может численно ускорить цепь.
    /// Для линейной ДНК и открытого кольца ограничиваем подвижность значением без красителя.
    /// У сверхспирали ускорение после перехода через точку релаксации физически допустимо.
    fn retarding_only(stained: f64, clean: impl FnOnce() -> f64, extension: f64) -> f64 {
        if extension == 1.0 {
            return stained;
        }
        stained.min(clean())
    }

    /// Полный множитель подвижности от EtBr объединяет изменение геометрии и дополнительное трение.
    pub fn etbr_total_factor(
        n_bp: f64,
        form: DnaForm,
        percent: f64,
        e_v_cm: f64,
        free_etbr_m: f64,
    ) -> f64 {
        Self::etbr_geometry_factor(n_bp, form, percent, e_v_cm, free_etbr_m)
            * crate::stain::friction_factor(crate::stain::nu_bound(free_etbr_m))
    }

    /// Отношение подвижностей одной формы с EtBr и без него для заданной пользователем кинетики.
    pub fn etbr_geometry_factor(
        n_bp: f64,
        form: DnaForm,
        percent: f64,
        e_v_cm: f64,
        free_etbr_m: f64,
    ) -> f64 {
        let clean = Self::form_mobility_ratio(n_bp, form, percent, e_v_cm, 0.0);
        let stained = Self::form_mobility_ratio(n_bp, form, percent, e_v_cm, free_etbr_m);
        stained / clean
    }
}
