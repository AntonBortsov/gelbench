//! Подвижности бромфенолового синего и ксиленцианола заданы по совместной миграции с ДНК
//! в 1% агарозе при 5 В/см и 25 °C (соответственно около 300 п.о. и 5 кб). С температурой они
//! меняются по вязкостной поправке ДНК. Состав буфера и просеивание красителей модель не учитывает.
//! При 0.5% агарозе расчётные эквиваленты (около 0.5 и 25 кб) расходятся со справочными (около 1 и 10 кб).

use serde::{Deserialize, Serialize};

use crate::mobility::{KT_PER_CELSIUS, T_REF_K};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DyeKind {
    /// Бромфеноловый синий: 670 Да, z ≈ −2 при pH геля.
    BromophenolBlue,
    /// Ксиленцианол FF: 539 Да, z ≈ −2.
    XyleneCyanol,
}

pub const MU_BPB_CM2_V_S: f64 = 3.1336e-4;
pub const MU_XC_CM2_V_S: f64 = 1.2386e-4;

/// Оценки Стокса–Эйнштейна для радиусов 0.57 и 0.54 нм, см²/с. Затруднение диффузии в геле опущено.
pub const D_BPB_CM2_S: f64 = 4.3e-6;
pub const D_XC_CM2_S: f64 = 4.6e-6;

/// Размеры ДНК для калибровки по совместной миграции при 1% агарозы.
pub const BPB_BP_EQUIVALENT: f64 = 300.0;
pub const XC_BP_EQUIVALENT: f64 = 5000.0;

impl DyeKind {
    pub fn mobility_cm2_v_s(&self, temp_k: f64) -> f64 {
        let theta = (1.0 + KT_PER_CELSIUS * (temp_k - T_REF_K)).max(0.0);
        match self {
            DyeKind::BromophenolBlue => MU_BPB_CM2_V_S * theta,
            DyeKind::XyleneCyanol => MU_XC_CM2_V_S * theta,
        }
    }

    pub fn diffusivity_cm2_s(&self) -> f64 {
        match self {
            DyeKind::BromophenolBlue => D_BPB_CM2_S,
            DyeKind::XyleneCyanol => D_XC_CM2_S,
        }
    }
}

/// Цветовое состояние БФБ по локальному pH.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BpbColor {
    Yellow,
    Green,
    Blue,
}

/// БФБ желтеет ниже pH 3 и синеет выше 4.6. Жёлтая полоса указывает на закисление в месте наблюдения.
pub fn bpb_color_at_ph(ph: f64) -> BpbColor {
    if ph < 3.0 {
        BpbColor::Yellow
    } else if ph <= 4.6 {
        BpbColor::Green
    } else {
        BpbColor::Blue
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dye_mobilities_match_comigration_anchor() {
        // Константы красителей ВЫВЕДЕНЫ из модели ДНК в точке якоря
        // (1 %, 5 В/см, 25 °C) по фольклору сомиграции. Страж связывает
        // их с источником: если модель ДНК изменят, тест упадёт и
        // константы придётся пересчитать ОСОЗНАННО, а не молча.
        let mu300 =
            crate::mobility2::MobilityV2::mobility_cm2_v_s(BPB_BP_EQUIVALENT, 1.0, 5.0, T_REF_K);
        let mu5000 =
            crate::mobility2::MobilityV2::mobility_cm2_v_s(XC_BP_EQUIVALENT, 1.0, 5.0, T_REF_K);
        assert!(
            (MU_BPB_CM2_V_S / mu300 - 1.0).abs() < 0.005,
            "µ_БФБ разошлась с моделью(300 бп): {MU_BPB_CM2_V_S:.5e}              против {mu300:.5e}"
        );
        assert!(
            (MU_XC_CM2_V_S / mu5000 - 1.0).abs() < 0.005,
            "µ_XC разошлась с моделью(5 кб): {MU_XC_CM2_V_S:.5e}              против {mu5000:.5e}"
        );
    }

    #[test]
    fn dye_equivalents_follow_vendor_table() {
        // ПРЕДСКАЗАНИЕ, а не калибровка: подвижность красителя постоянна,
        // µ ДНК падает с концентрацией геля, значит «эквивалент в п.о.»
        // обязан уменьшаться с C. Сверяем с таблицами производителей
        // (окна широкие — таблицы сами дают порядок величины).
        let equiv = |mu_dye: f64, percent: f64| {
            let (mut lo, mut hi) = (20.0_f64, 3.0e5_f64);
            for _ in 0..80 {
                let m = (lo * hi).sqrt();
                let mu = crate::mobility2::MobilityV2::mobility_cm2_v_s(m, percent, 5.0, T_REF_K);
                if mu > mu_dye {
                    lo = m;
                } else {
                    hi = m;
                }
            }
            (lo * hi).sqrt()
        };
        for (percent, bpb_lo, bpb_hi, xc_lo, xc_hi) in [
            (1.0, 280.0, 320.0, 4500.0, 5500.0), // точка якоря
            (1.5, 150.0, 300.0, 1500.0, 3000.0), // таблица ~200 бп / ~2 кб
            (2.0, 110.0, 240.0, 900.0, 2000.0),  // таблица ~150 бп / ~1.5 кб
        ] {
            let b = equiv(MU_BPB_CM2_V_S, percent);
            let x = equiv(MU_XC_CM2_V_S, percent);
            assert!(
                (bpb_lo..=bpb_hi).contains(&b),
                "БФБ при {percent} %: эквивалент {b:.0} бп вне таблицы"
            );
            assert!(
                (xc_lo..=xc_hi).contains(&x),
                "XC при {percent} %: эквивалент {x:.0} бп вне таблицы"
            );
        }
        // Монотонность по C — само существо предсказания.
        assert!(equiv(MU_XC_CM2_V_S, 0.7) > equiv(MU_XC_CM2_V_S, 1.0));
    }

    #[test]
    fn bpb_color_thresholds() {
        assert_eq!(bpb_color_at_ph(8.3), BpbColor::Blue);
        assert_eq!(bpb_color_at_ph(4.7), BpbColor::Blue);
        assert_eq!(bpb_color_at_ph(4.6), BpbColor::Green);
        assert_eq!(bpb_color_at_ph(3.5), BpbColor::Green);
        assert_eq!(bpb_color_at_ph(2.9), BpbColor::Yellow);
    }

    #[test]
    fn dye_order_bpb_ahead_of_xc() {
        // БФБ всегда впереди XC (меньше «эквивалент» — быстрее): порядок
        // полос на геле, который знает каждый студент.
        assert!(MU_BPB_CM2_V_S > MU_XC_CM2_V_S * 1.3);
    }
}
