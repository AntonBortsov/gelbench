//! Deterministic agarose gel electrophoresis.
//!
//! Create a [`SimState`] from [`Params`], then call [`SimState::step`] and [`SimState::observe`].

pub mod buffers;
pub mod chamber;
pub mod circuit;
pub mod dyes;
mod engine;
pub mod fmath;
pub mod loading;
pub mod materials;
pub mod medium;
pub mod medium_loop;
pub mod mobility;
pub mod mobility2;
pub mod params;
pub mod reservoir;
pub mod rng;
pub mod stain;
pub mod thermal;
pub mod transport;

pub use engine::{select_subcolumn_count, BufferActionError, Observation, SimFault, SimState};
pub use mobility::MobilityModel;
pub use params::{GelParams, Params, Topology};
