//! Параметры эксперимента и пространственной сетки.

use serde::{Deserialize, Serialize};

/// Размеры сетки концентраций.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct Topology {
    /// Ячеек вдоль направления миграции.
    pub cells_per_lane: usize,
    pub lanes: usize,
    /// Поперечных подстолбцов на дорожку.
    pub subcolumns: usize,
    pub populations: usize,
}

impl Default for Topology {
    fn default() -> Self {
        Self {
            cells_per_lane: 1024,
            lanes: 12,
            subcolumns: 16,
            populations: 8,
        }
    }
}

impl Topology {
    /// Число элементов массива концентраций.
    pub fn total_cells(&self) -> usize {
        self.lanes * self.populations * self.subcolumns * self.cells_per_lane
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Params {
    pub topology: Topology,
    /// Номинальное поле, В/см.
    pub field_v_per_cm: f32,
    /// Длина геля, см.
    pub gel_length_cm: f32,
    /// Длительность макрошага, с.
    pub macro_dt_s: f32,
    /// Расстояние от катодного края до центра лунки, см.
    pub well_position_cm: f32,
    pub gel: GelParams,
    /// Размеры ДНК в парах оснований. Число фрагментов и красителей в сумме равно `topology.populations`.
    pub population_sizes_bp: Vec<f64>,
    /// Красители загрузки следуют за ДНК-популяциями.
    pub dyes: Vec<crate::dyes::DyeKind>,
    /// Формы в порядке `population_sizes_bp`. Пустой список означает линейную ДНК.
    pub population_forms: Vec<crate::mobility2::DnaForm>,
    /// Открытый анод выпускает полосы в буфер. При закрытом аноде масса остаётся в крайней ячейке.
    /// Катодная граница всегда закрыта.
    pub open_anode: bool,
    /// Поправка к подвижности ДНК и красителей для конкретной партии геля и буфера.
    pub mobility_scale: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct GelParams {
    /// Агароза, % массы на объём.
    pub percent: f64,
    /// Исходная температура, K.
    pub temp_k: f64,
    /// Электроосмотическая подвижность, см²/(В·с). Поток к катоду замедляет все полосы.
    /// Паспортный показатель агарозы −mr = µ_EEO/(µ_EEO + µ_БСА), обычно 0.05–0.13 для низко-EEO марок.
    pub mu_eeo_cm2_v_s: f64,
}

impl Default for GelParams {
    fn default() -> Self {
        Self {
            percent: 1.0,
            // Опорная температура модели подвижности.
            temp_k: 298.15,
            mu_eeo_cm2_v_s: 0.0,
        }
    }
}

impl Default for Params {
    fn default() -> Self {
        Self {
            topology: Topology::default(),
            field_v_per_cm: 5.0,
            gel_length_cm: 10.0,
            macro_dt_s: 0.5,
            well_position_cm: 0.25,
            gel: GelParams::default(),
            population_sizes_bp: vec![
                100.0, 300.0, 500.0, 1000.0, 3000.0, 5000.0, 10000.0, 23000.0,
            ],
            dyes: Vec::new(),
            population_forms: Vec::new(),
            open_anode: false,
            mobility_scale: 1.0,
        }
    }
}
