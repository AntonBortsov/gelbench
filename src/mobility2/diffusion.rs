//! Longitudinal diffusion combines a two-exponent thermal fit (Pluen, joined at 6 kb)
//! with field dispersion from Meistermann & Tinland (shortest measured fragment 5.7 kb).

use super::{
    MobilityV2, DISP_A0, DISP_ALPHA, DISP_XC, D_1KB_CM2_S, D_BREAK_BP, D_EXP_REPT, D_EXP_ROUSE,
};
use crate::fmath::DetMath;

impl MobilityV2 {
    /// Longitudinal diffusion coefficient (cm²/s) at the given agarose percent and field in V/cm.
    pub fn dna_diffusivity_cm2_s(n_bp: f64, percent: f64, e_v_cm: f64) -> f64 {
        Self::dna_diffusivity_with(
            D_1KB_CM2_S,
            DISP_A0,
            DISP_ALPHA,
            DISP_XC,
            n_bp,
            percent,
            e_v_cm,
        )
    }
    /// Zero-field diffusivity (cm²/s) using a continuous two-exponent fit.
    pub fn thermal_diffusivity_with(d1: f64, n_bp: f64) -> f64 {
        if n_bp <= D_BREAK_BP {
            d1 * (n_bp / 1000.0).dpowf(-D_EXP_ROUSE)
        } else {
            d1 * (D_BREAK_BP / 1000.0).dpowf(-D_EXP_ROUSE) * (n_bp / D_BREAK_BP).dpowf(-D_EXP_REPT)
        }
    }

    pub fn dna_diffusivity_with(
        d1: f64,
        a0: f64,
        alpha: f64,
        xc: f64,
        n_bp: f64,
        percent: f64,
        e_v_cm: f64,
    ) -> f64 {
        let d0 = Self::thermal_diffusivity_with(d1, n_bp);
        let e = e_v_cm.abs();
        let nb = Self::blobs(n_bp, percent);
        if e <= 0.0 || nb <= 0.0 {
            return d0;
        }
        // Подавление дисперсии ниже одного блоба — экстраполяция за пределы измерений.
        let s = 1.0 - (-nb * nb).dexp();
        let x = e * nb;
        d0 + a0 * percent.dpowf(-alpha) * e.dpowf(1.5) * (1.0 + xc / x).sqrt() * s
    }
}
