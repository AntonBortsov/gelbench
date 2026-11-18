//! Linear migration fits use Holmes & Stellwagen and Rill. Diffusion and longitudinal dispersion
//! use Pluen and Meistermann & Tinland. The trapping and plasmid coefficients are empirical.

/// Characteristic void size of a 1% gel (nm), fitted independently of the sieve fiber radius.
pub const A1_NM: f64 = 333.47;
/// Exponent for the decrease of void size with agarose concentration.
pub const PORE_EXP: f64 = 0.7;
/// Effective fiber bundle radius in the Ogston branch (nm).
pub const R_FIBER_NM: f64 = 5.8714;
/// Agarose density (g/cm³) for the fiber volume fraction.
pub const RHO_AGAROSE_G_CM3: f64 = 1.64;
/// Kuhn length of double-stranded DNA (nm).
pub const B_KUHN_NM: f64 = 100.0;
/// Ratio of free-solution mobility to the high-field BRF plateau.
pub const A_SAT: f64 = 3.0;
/// Reduced field per V/cm in a 1% gel.
pub const EPS0_PER_V_CM: f64 = 0.21031;
/// Orientation crossover scale in units of blobs.
pub const K_ORIENT: f64 = 0.27171;
/// Exponent for movement of the orientation crossover with applied field.
pub const Q_ORIENT: f64 = 0.80;
/// Width of the Ogston-to-reptation crossover in log10 base pairs.
pub const W_DEC: f64 = 0.16464;
/// Ratio of gyration radius to void size at crossover.
pub const X_CROSS: f64 = 0.21567;

/// Upper bound on the fractional mobility penalty from entropic trapping.
pub const A_TRAP: f64 = 0.55;
/// Center of the trapping window in blobs.
pub const N_TRAP: f64 = 4.4402;
/// Width of the trapping window in log blob count.
pub const W_TRAP: f64 = 0.5692;
/// Agarose concentration (%) at half the limiting trapping strength.
pub const C_TRAP_HALF: f64 = 0.77594;
/// Field (V/cm) used to attenuate the trapping term.
pub const E_TRAP_OFF_V_CM: f64 = 2.5;
/// Fitted thermal diffusion coefficient for 1 kb (cm²/s). See Pluen et al. and the dispersion fit.
pub const D_1KB_CM2_S: f64 = 1.1323e-8;
/// Low-mass exponent of thermal diffusion.
pub const D_EXP_ROUSE: f64 = 0.52;
/// High-mass exponent of thermal diffusion.
pub const D_EXP_REPT: f64 = 1.55;
/// Continuous join of the two thermal diffusion branches (bp).
pub const D_BREAK_BP: f64 = 6000.0;
/// Field-induced dispersion scale at 1 V/cm and 1% agarose (cm²/s).
pub const DISP_A0: f64 = 7.4147e-10;
/// Exponent relating field-induced dispersion to gel concentration.
pub const DISP_ALPHA: f64 = 0.56704;
/// Crossover in field strength times blob count (V/cm).
pub const DISP_XC: f64 = 8.4343;

/// Gyration radius of a compact supercoil relative to linear DNA.
pub const RG_SC_RATIO: f64 = 0.55;
/// Effective ring radius relative to linear DNA in the impalement model.
pub const RG_CIRC_RATIO: f64 = 1.25;
/// Magnitude of superhelical density at full model compaction.
pub const SIGMA_COMPACT: f64 = 0.06;
/// Relative free-solution mobility scale for supercoiled material.
pub const K_FREE_SC: f64 = 1.1;
/// Relative free-solution mobility scale for relaxed circles.
pub const K_FREE_CIRC: f64 = 0.95;
/// Dimensionless coefficient of fiber impalement for closed forms.
pub const K_IMPALE: f64 = 0.6;
/// Tunable coefficients for evaluating alternate calibrated migration curves.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct V2Params {
    pub r_fiber_nm: f64,
    pub a1_nm: f64,
    pub x_cross: f64,
    pub eps0_per_v_cm: f64,
    pub k_orient: f64,
    pub q_orient: f64,
    pub w_dec: f64,
    pub a_trap: f64,
    pub n_trap: f64,
    pub w_trap: f64,
    pub c_trap_half: f64,
    pub join_p: f64,
}

impl V2Params {
    pub const CURRENT: Self = Self {
        r_fiber_nm: R_FIBER_NM,
        a1_nm: A1_NM,
        x_cross: X_CROSS,
        eps0_per_v_cm: EPS0_PER_V_CM,
        k_orient: K_ORIENT,
        q_orient: Q_ORIENT,
        w_dec: W_DEC,
        a_trap: A_TRAP,
        n_trap: N_TRAP,
        w_trap: W_TRAP,
        c_trap_half: C_TRAP_HALF,
        join_p: JOIN_P,
    };
}

/// Exponent used to join the two BRF branches. The current fitted model uses their linear sum.
pub const JOIN_P: f64 = 1.0;
