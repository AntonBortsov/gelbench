//! Log-mobility blend of Ogston sieving and biased reptation with field-dependent orientation and trapping.
//! Fiber-bundle radius controls sieving. The effective void size sets the reptation scale.

use super::{
    MigrationRegime, MobilityV2, V2Params, A_SAT, B_KUHN_NM, E_TRAP_OFF_V_CM, PORE_EXP,
    RHO_AGAROSE_G_CM3,
};
use crate::fmath::DetMath;
use crate::materials::agarose::rg_nm;
use crate::mobility::{KT_PER_CELSIUS, MU0_CM2_V_S, T_REF_K};

impl MobilityV2 {
    /// Размер пустоты для рептации, нм.
    pub fn pore_nm(percent: f64) -> f64 {
        Self::pore_nm_with(&V2Params::CURRENT, percent)
    }

    pub fn pore_nm_with(p: &V2Params, percent: f64) -> f64 {
        p.a1_nm * percent.dpowf(-PORE_EXP)
    }

    /// Длина линейной ДНК (п.о.) с заданным радиусом инерции (нм).
    pub fn bp_at_rg(rg_target_nm: f64) -> f64 {
        let (mut lo, mut hi) = (10.0_f64, 1.0e7_f64);
        for _ in 0..80 {
            let mid = (lo * hi).sqrt();
            if rg_nm(mid) < rg_target_nm {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        (lo * hi).sqrt()
    }

    /// Центр перехода от Огстона к рептации (п.о.), где Rg = X_CROSS·a(C).
    pub fn crossover_bp(percent: f64) -> f64 {
        Self::crossover_bp_with(&V2Params::CURRENT, percent)
    }

    pub fn crossover_bp_with(p: &V2Params, percent: f64) -> f64 {
        Self::bp_at_rg(p.x_cross * Self::pore_nm_with(p, percent))
    }

    /// Приведённое поле ε(E, C).
    pub fn epsilon(e_v_cm: f64, percent: f64) -> f64 {
        Self::epsilon_with(&V2Params::CURRENT, e_v_cm, percent)
    }

    pub fn epsilon_with(p: &V2Params, e_v_cm: f64, percent: f64) -> f64 {
        p.eps0_per_v_cm * e_v_cm * percent.dpowf(-1.4)
    }

    /// Число блобов N(bp, C).
    pub fn blobs(n_bp: f64, percent: f64) -> f64 {
        Self::blobs_with(&V2Params::CURRENT, n_bp, percent)
    }

    pub fn blobs_with(p: &V2Params, n_bp: f64, percent: f64) -> f64 {
        let a = Self::pore_nm_with(p, percent);
        n_bp * crate::materials::agarose::BP_NM * B_KUHN_NM / (a * a)
    }

    /// Гладкое окно энтропийного захвата по числу блобов. Подвижность снижается максимум на долю A_TRAP.
    pub fn trapping_factor(n_blobs: f64, percent: f64, e_v_cm: f64) -> f64 {
        Self::trapping_factor_with(&V2Params::CURRENT, n_blobs, percent, e_v_cm)
    }

    pub fn trapping_factor_with(p: &V2Params, n_blobs: f64, percent: f64, e_v_cm: f64) -> f64 {
        let x = e_v_cm / E_TRAP_OFF_V_CM;
        let g_e = 1.0 / (1.0 + x * x);
        let g_c = percent / (percent + p.c_trap_half);
        let ln_ratio = (n_blobs / p.n_trap).dln();
        let g_n = (-ln_ratio * ln_ratio / (2.0 * p.w_trap * p.w_trap)).dexp();
        1.0 - p.a_trap * g_e * g_c * g_n
    }

    /// μ/μ₀ для заданного Rg. Эффективная длина управляет весами ветвей, фактическая — числом блобов.
    pub(super) fn ratio_core(rg: f64, l_eff_bp: f64, n_bp: f64, percent: f64, e_v_cm: f64) -> f64 {
        Self::ratio_core_with(&V2Params::CURRENT, rg, l_eff_bp, n_bp, percent, e_v_cm)
    }

    pub fn ratio_core_with(
        p: &V2Params,
        rg: f64,
        l_eff_bp: f64,
        n_bp: f64,
        percent: f64,
        e_v_cm: f64,
    ) -> f64 {
        let lx = Self::crossover_bp_with(p, percent);
        Self::ratio_core_with_lx(p, lx, rg, l_eff_bp, n_bp, percent, e_v_cm)
    }

    /// Кроссовер L× зависит только от концентрации и параметров. При подгонке его можно переиспользовать
    /// для всех длин ДНК одной концентрации геля.
    #[allow(clippy::too_many_arguments)]
    pub fn ratio_core_with_lx(
        p: &V2Params,
        lx: f64,
        rg: f64,
        l_eff_bp: f64,
        n_bp: f64,
        percent: f64,
        e_v_cm: f64,
    ) -> f64 {
        // Объёмная доля волокон по массовому балансу агарозы задаёт наклон графика Фергюсона.
        let ln_og = -(percent / (100.0 * RHO_AGAROSE_G_CM3)) * (1.0 + rg / p.r_fiber_nm).powi(2);
        let nb = Self::blobs_with(p, n_bp, percent);
        let eps = Self::epsilon_with(p, e_v_cm, percent).max(1e-12);
        let f = 2.0 * eps / (5.0 + 2.0 * A_SAT * eps);
        let h = nb / (nb + Self::n_star_blobs_with(p, e_v_cm, percent));
        let a = 1.0 / (3.0 * nb);
        let brf = if p.join_p == 2.0 {
            (a.powi(2) + (f * h) * (f * h)).sqrt().min(1.0)
        } else {
            (a.dpowf(p.join_p) + (f * h).dpowf(p.join_p))
                .dpowf(1.0 / p.join_p)
                .min(1.0)
        };
        let w = 1.0 / (1.0 + (-(l_eff_bp.dlog10() - lx.dlog10()) / p.w_dec).dexp());
        let mu = ((1.0 - w) * ln_og + w * brf.dln()).dexp();
        mu * Self::trapping_factor_with(p, nb, percent, e_v_cm)
    }

    pub fn mobility_ratio_with(p: &V2Params, n_bp: f64, percent: f64, e_v_cm: f64) -> f64 {
        Self::ratio_core_with(p, rg_nm(n_bp), n_bp, n_bp, percent, e_v_cm)
    }

    pub fn mobility_ratio_with_lx(
        p: &V2Params,
        lx: f64,
        rg: f64,
        n_bp: f64,
        percent: f64,
        e_v_cm: f64,
    ) -> f64 {
        Self::ratio_core_with_lx(p, lx, rg, n_bp, n_bp, percent, e_v_cm)
    }

    /// Подвижность относительно свободного раствора при опорной температуре.
    pub fn mobility_ratio(n_bp: f64, percent: f64, e_v_cm: f64) -> f64 {
        assert!(n_bp > 0.0 && percent > 0.0 && e_v_cm >= 0.0);
        Self::mobility_ratio_with(&V2Params::CURRENT, n_bp, percent, e_v_cm)
    }

    /// Абсолютная подвижность при заданной температуре, см²/(В·с).
    pub fn mobility_cm2_v_s(n_bp: f64, percent: f64, e_v_cm: f64, temp_k: f64) -> f64 {
        let theta = (1.0 + KT_PER_CELSIUS * (temp_k - T_REF_K)).max(0.0);
        MU0_CM2_V_S * Self::mobility_ratio(n_bp, percent, e_v_cm) * theta
    }

    /// Качественный режим переноса. Окно захвата отмечается при снижении подвижности минимум на 10%.
    pub fn regime(n_bp: f64, percent: f64, e_v_cm: f64) -> MigrationRegime {
        let nb = Self::blobs(n_bp, percent);
        if Self::trapping_factor(nb, percent, e_v_cm) < 0.9 {
            return MigrationRegime::EntropicTrapping;
        }
        if n_bp < Self::crossover_bp(percent) {
            return MigrationRegime::Ogston;
        }
        if nb < Self::n_star_blobs(e_v_cm, percent) {
            MigrationRegime::Reptation
        } else {
            MigrationRegime::OrientedReptation
        }
    }

    /// Число блобов при переходе к ориентированной рептации. Heller (1994) наблюдает N* ∝ E⁻⁰·⁶⁵…E⁻⁰·⁸⁰.
    pub fn n_star_blobs(e_v_cm: f64, percent: f64) -> f64 {
        Self::n_star_blobs_with(&V2Params::CURRENT, e_v_cm, percent)
    }

    pub fn n_star_blobs_with(p: &V2Params, e_v_cm: f64, percent: f64) -> f64 {
        p.k_orient
            / Self::epsilon_with(p, e_v_cm, percent)
                .max(1e-12)
                .dpowf(p.q_orient)
    }
}
