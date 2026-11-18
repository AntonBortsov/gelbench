//! Calibrated DNA mobility and diffusion in agarose. Gel concentrations are percentages (w/v),
//! lengths are in base pairs, fields in V/cm and temperatures in kelvin.

mod calibration;
mod diffusion;
mod forms;
mod linear;
#[cfg(test)]
mod tests;

pub use calibration::{
    V2Params, A1_NM, A_SAT, A_TRAP, B_KUHN_NM, C_TRAP_HALF, DISP_A0, DISP_ALPHA, DISP_XC,
    D_1KB_CM2_S, D_BREAK_BP, D_EXP_REPT, D_EXP_ROUSE, EPS0_PER_V_CM, E_TRAP_OFF_V_CM, JOIN_P,
    K_FREE_CIRC, K_FREE_SC, K_IMPALE, K_ORIENT, N_TRAP, PORE_EXP, Q_ORIENT, RG_CIRC_RATIO,
    RG_SC_RATIO, RHO_AGAROSE_G_CM3, R_FIBER_NM, SIGMA_COMPACT, W_DEC, W_TRAP, X_CROSS,
};
pub use forms::DnaForm;

/// Dominant linear-chain migration regime at the specified field and gel concentration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MigrationRegime {
    Ogston,
    Reptation,
    OrientedReptation,
    EntropicTrapping,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MobilityV2;
