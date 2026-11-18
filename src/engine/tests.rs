use super::{select_subcolumn_count, SimState};
use crate::{dyes::DyeKind, loading, Params, Topology};

fn single_lane() -> SimState {
    SimState::new(Params {
        topology: Topology {
            cells_per_lane: 128,
            lanes: 1,
            subcolumns: 1,
            populations: 2,
        },
        population_sizes_bp: vec![1000.0],
        dyes: vec![DyeKind::BromophenolBlue],
        ..Params::default()
    })
}

#[test]
fn reversing_drift_preserves_transverse_resolution() {
    let theta = [0.8, 1.0, 1.2];
    let columns = |velocity| select_subcolumn_count(&theta, 0.10, velocity, 100.0, 0.10);
    assert!(columns(0.01) > 1);
    assert_eq!(columns(-0.01), columns(0.01));
    assert_eq!(columns(0.0), 1);
}

#[test]
fn transport_uses_the_uniformity_of_the_active_field() {
    use crate::medium::{IonicStrengthModel, TAE_MEDIUM};
    let mut external = single_lane();
    external.set_field_phi(&[1.2; 128]);
    external.enable_medium(TAE_MEDIUM, IonicStrengthModel::OnsagerScreened);
    let mut control = single_lane();
    control.enable_medium(TAE_MEDIUM, IonicStrengthModel::OnsagerScreened);
    for _ in 0..3 {
        external.step();
        control.step();
    }
    assert!(external.time_s() > 0.0);
    assert_eq!(external.concentrations(), control.concentrations());
}

#[test]
fn reloaded_bands_start_at_the_well_after_a_fractional_drift() {
    let mut fresh = single_lane();
    let mut running = single_lane();
    let dye = running.dye_population_index(0);
    running.set_population_kinetics(0, 0.04, 0.0);
    running.set_population_kinetics(dye, 0.04, 0.0);
    running.step();
    assert!(running.band_com_cm(0, 0) > fresh.band_com_cm(0, 0));
    let sample = loading::SampleLoad {
        volume_ul: 10.0,
        masses_ng: vec![10.0],
        densifier: loading::Densifier::GlycerolOrFicoll,
        residual_ethanol: false,
        prep_time_s: 0.0,
    };
    // Reload the dye first to check that this operation alone consumes its own offset.
    for run in [&mut fresh, &mut running] {
        run.load_dyes(0, &loading::Comb::BIO_RAD_20, &[5.0]);
    }
    assert_eq!(fresh.band_com_cm(0, dye), running.band_com_cm(0, dye));
    running.step();
    for run in [&mut fresh, &mut running] {
        run.load_sample(0, &loading::Comb::BIO_RAD_20, &sample);
    }
    assert_eq!(fresh.band_com_cm(0, 0), running.band_com_cm(0, 0));
    assert_eq!(fresh.band_mass(0, 0), running.band_mass(0, 0));
}
