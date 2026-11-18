//! Загрузка задаёт начальную массу полосы в нанограммах на ячейку.
//! Утечка до включения поля и перелив в соседние дорожки учитываются раздельно.
//! Пороги перегрузки взяты из рекомендаций Thermo/NCI/Unicorn (0,1–0,2 мкг/мм ширины,
//! 500 нг на полосу). Они классифицируют риск шлейфа, его форма здесь не моделируется.

use crate::fmath::DetMath;
use serde::{Deserialize, Serialize};

/// Геометрия зубца задаёт объём лунки и начальную протяжённость полосы вдоль поля.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Comb {
    pub teeth: usize,
    /// Ширина зубца поперёк поля, м.
    pub tooth_width_m: f64,
    /// Толщина зубца вдоль поля, м.
    pub tooth_thickness_m: f64,
    pub tooth_depth_m: f64,
}

impl Comb {
    /// Гребень Bio-Rad на 20 лунок, каталожная ёмкость 18,2 мкл (по геометрии 18,15 мкл).
    pub const BIO_RAD_20: Comb = Comb {
        teeth: 20,
        tooth_width_m: 4.84e-3,
        tooth_thickness_m: 0.75e-3,
        tooth_depth_m: 5.0e-3,
    };

    /// Ёмкость лунки V_cap = b·t·d, мкл.
    pub fn well_capacity_ul(&self) -> f64 {
        self.tooth_width_m * self.tooth_thickness_m * self.tooth_depth_m * 1.0e9
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Densifier {
    /// 6× краситель с глицерином/сахарозой/Ficoll — проба тяжелее буфера.
    GlycerolOrFicoll,
    /// Проба без уплотнителя близка по плотности к буферу.
    None,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SampleLoad {
    pub volume_ul: f64,
    /// Массы по ДНК-популяциям, нг (параллелен population_sizes_bp).
    pub masses_ng: Vec<f64>,
    pub densifier: Densifier,
    /// Остаточный этанол после очистки может вызвать всплытие даже при наличии глицерина (QIAGEN FAQ 205).
    pub residual_ethanol: bool,
    /// Время утечки между загрузкой и включением поля, с.
    pub prep_time_s: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OverloadStatus {
    /// ≤ 0.2 мкг/мм ширины лунки и ≤ 500 нг/полосу.
    Ok,
    /// Масса превышает рекомендацию 0,2 мкг/мм ширины.
    Heavy,
    /// Масса превышает порог риска шлейфа (500 нг).
    SmearRisk,
}

/// Баланс массы пробы после загрузки, нг.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadReport {
    /// Масса в дорожке после утечки из лунки.
    pub loaded_ng: f64,
    /// Потеря при всплытии пробы.
    pub density_loss_ng: f64,
    /// Суммарный перелив в соседние дорожки.
    pub spilled_to_neighbors_ng: f64,
    /// Перелив в буфер над гелем.
    pub spilled_lost_ng: f64,
    /// Статус перегрузки по фактически доехавшей массе.
    pub overload: OverloadStatus,
    /// Доля объёма сверх ёмкости.
    pub overflow_fraction: f64,
}

/// Плотность буфера, г/см³ (разбавленный водный буфер).
pub const RHO_BUFFER: f64 = 1.000;
/// Оценка плотности пробы с уплотнителем, г/см³ (около 5% глицерина в лунке).
pub const RHO_WITH_DENSIFIER: f64 = 1.015;
/// Оценка плотности пробы без уплотнителя, г/см³.
pub const RHO_PLAIN: f64 = 0.995;
/// Некалиброванная поправка к плотности пробы при наличии остаточного этанола, г/см³.
pub const RHO_ETHANOL_SHIFT: f64 = -0.050;
/// Некалиброванный масштаб времени для утечки, с.
pub const TAU0_LEAK_S: f64 = 10.0;

/// Плотность пробы, г/см³.
pub fn sample_density(load: &SampleLoad) -> f64 {
    let base = match load.densifier {
        Densifier::GlycerolOrFicoll => RHO_WITH_DENSIFIER,
        Densifier::None => RHO_PLAIN,
    };
    base + if load.residual_ethanol {
        RHO_ETHANOL_SHIFT
    } else {
        0.0
    }
}

/// Доля всплывшей пробы до старта поля по экспоненциальной модели утечки.
pub fn density_loss_fraction(load: &SampleLoad) -> f64 {
    let d_rho = sample_density(load) - RHO_BUFFER;
    if d_rho >= 0.0 {
        return 0.0;
    }
    let tau = TAU0_LEAK_S / (d_rho.abs() / RHO_BUFFER);
    1.0 - (-load.prep_time_s / tau).dexp()
}

/// Статус перегрузки по массе `total_ng` в лунке шириной `width_m`.
pub fn overload_status(total_ng: f64, width_m: f64) -> OverloadStatus {
    let per_band_limit_ng = 500.0;
    let per_mm_heavy_ng = 200.0; // 0.2 мкг/мм — верх рабочей рекомендации
    if total_ng > per_band_limit_ng {
        OverloadStatus::SmearRisk
    } else if total_ng > per_mm_heavy_ng * (width_m * 1.0e3) {
        OverloadStatus::Heavy
    } else {
        OverloadStatus::Ok
    }
}

/// Каждый сосед получает четверть избытка. Остаток уходит в буфер над гелем.
/// Возвращает доли полной массы в лунке, на одного существующего соседа и в буфере.
pub fn overflow_split(
    volume_ul: f64,
    capacity_ul: f64,
    has_left: bool,
    has_right: bool,
) -> (f64, f64, f64) {
    assert!(volume_ul > 0.0 && capacity_ul > 0.0);
    if volume_ul <= capacity_ul {
        return (1.0, 0.0, 0.0);
    }
    let f_over = (volume_ul - capacity_ul) / volume_ul;
    let stay = 1.0 - f_over;
    let per_side = 0.25 * f_over;
    let n_sides = has_left as u32 + has_right as u32;
    let lost = f_over - per_side * n_sides as f64;
    (stay, per_side, lost)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn anchor20_capacity_matches_bio_rad_catalog() {
        // V_cap = b·t·d: 4.84 × 0.75 × 5.0 мм = 18.15 мкл против
        // каталожных 18.2 (M4-PLAN §4.1 — формула подтверждена до мкл).
        let v = Comb::BIO_RAD_20.well_capacity_ul();
        assert!((v - 18.15).abs() < 0.01, "V_cap = {v:.3}");
        assert!((v - 18.2).abs() / 18.2 < 0.005, "каталог разошёлся: {v:.3}");
    }

    #[test]
    fn overflow_split_conserves_and_reports() {
        // Баланс долей = 1 точно; перелив делится по правилу §11-Р.1.
        for &(v, has_l, has_r) in &[
            (25.0, true, true),
            (40.0, true, false),
            (19.0, false, false),
        ] {
            let (stay, per_side, lost) = overflow_split(v, 18.15, has_l, has_r);
            let n = has_l as u32 + has_r as u32;
            let total = stay + per_side * n as f64 + lost;
            assert!((total - 1.0).abs() < 1e-15, "доли не сходятся: {total}");
            assert!(stay > 0.0 && lost >= 0.0 && per_side >= 0.0);
        }
        // Без перелива — всё в лунке.
        assert_eq!(overflow_split(10.0, 18.15, true, true), (1.0, 0.0, 0.0));
        // Крайняя лунка теряет долю несуществующего соседа.
        let (_, ps_both, lost_both) = overflow_split(30.0, 18.15, true, true);
        let (_, ps_edge, lost_edge) = overflow_split(30.0, 18.15, true, false);
        assert_eq!(ps_both, ps_edge);
        assert!(lost_edge > lost_both);
    }

    #[test]
    fn anchor21_three_density_cases() {
        // №21: уплотнитель → ноль потерь; без него — утечка растёт со
        // временем; этанол топит даже глицерин; монотонность по |Δρ|.
        let mk = |dens, eth, t| SampleLoad {
            volume_ul: 10.0,
            masses_ng: vec![100.0],
            densifier: dens,
            residual_ethanol: eth,
            prep_time_s: t,
        };
        assert_eq!(
            density_loss_fraction(&mk(Densifier::GlycerolOrFicoll, false, 600.0)),
            0.0,
            "глицерин обязан держать пробу"
        );
        let plain_1 = density_loss_fraction(&mk(Densifier::None, false, 60.0));
        let plain_5 = density_loss_fraction(&mk(Densifier::None, false, 300.0));
        assert!(
            plain_1 > 0.0 && plain_5 > plain_1,
            "утечка не растёт со временем"
        );
        let eth_gly = density_loss_fraction(&mk(Densifier::GlycerolOrFicoll, true, 120.0));
        assert!(
            eth_gly > 0.3,
            "этанол обязан топить даже глицериновую пробу: {eth_gly:.3}"
        );
        // Монотонность по |Δρ|: этанол без глицерина легче ⇒ течёт быстрее.
        let eth_plain = density_loss_fraction(&mk(Densifier::None, true, 120.0));
        assert!(eth_plain > eth_gly);
    }

    #[test]
    fn overload_thresholds_from_literature() {
        let w = Comb::BIO_RAD_20.tooth_width_m; // 4.84 мм → Heavy c ~968 нг
        assert_eq!(overload_status(100.0, w), OverloadStatus::Ok);
        assert_eq!(overload_status(501.0, w), OverloadStatus::SmearRisk);
        // Узкая лунка 1 мм: Heavy раньше SmearRisk.
        assert_eq!(overload_status(300.0, 1.0e-3), OverloadStatus::Heavy);
        assert_eq!(overload_status(150.0, 1.0e-3), OverloadStatus::Ok);
    }
}
