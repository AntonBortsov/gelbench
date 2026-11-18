//! Электролиз воды в двух перемешанных электродных отсеках без пространственного разрешения.
//! Анод накапливает H⁺, катод — OH⁻: 10,364 мкмоль на кулон на каждом электроде.
//! Противоионы приходят через буферный мост. Обмен массой между гелем и отсеком
//! в балансе не учитывается (типичные объёмы — около 1 и 150 мл соответственно).

use crate::medium::{MediumSpec, NodeComposition, StrongIons};

pub use crate::medium::FARADAY as FARADAY_C_MOL;

/// Состояние одного перемешанного электродного отсека.
#[derive(Debug, Clone)]
pub struct ReservoirBook {
    /// Состав ящика: семейства буфера + накопленный титрант.
    pub node: NodeComposition,
    /// Сильные ионы пробы в резервуаре.
    pub strong: StrongIons,
    pub volume_l: f64,
    /// Накопленный модуль заряда в кулонах.
    pub charge_c: f64,
    ph_warm: f64,
}

impl ReservoirBook {
    pub fn fresh(spec: &MediumSpec, volume_l: f64) -> Self {
        assert!(volume_l > 0.0, "объём книги обязан быть положительным");
        let node = spec.fresh_node();
        let ph0 = node.speciate().ph;
        Self {
            node,
            strong: StrongIons::ZERO,
            volume_l,
            charge_c: 0.0,
            ph_warm: ph0,
        }
    }

    /// Анодный заряд `q_c` Кл (> 0): +q/(F·V) моль/л сильной кислоты.
    pub fn add_anodic_charge(&mut self, q_c: f64) {
        assert!(q_c >= 0.0);
        self.charge_c += q_c;
        self.node.titrant_strong_acid_m += q_c / (FARADAY_C_MOL * self.volume_l);
    }

    /// Катодный заряд, Кл. Щёлочь уменьшает концентрацию кислотного титранта.
    pub fn add_cathodic_charge(&mut self, q_c: f64) {
        assert!(q_c >= 0.0);
        self.charge_c += q_c;
        self.node.titrant_strong_acid_m -= q_c / (FARADAY_C_MOL * self.volume_l);
    }

    /// pH перемешанного отсека. При выходе корня за диапазон решателя метод возвращает
    /// 0 или 14 по знаку титранта. На границе это признак истощения буфера.
    pub fn ph(&mut self) -> f64 {
        match self.node.try_speciate_warm_with(&self.strong, self.ph_warm) {
            Ok(eq) => {
                self.ph_warm = eq.ph;
                eq.ph
            }
            Err(_) => {
                let edge = if self.node.titrant_strong_acid_m > 0.0 {
                    0.0
                } else {
                    14.0
                };
                self.ph_warm = edge;
                edge
            }
        }
    }

    /// Остаток баланса заряда при текущем pH, моль экв/л.
    pub fn electroneutrality_residual_m(&mut self) -> f64 {
        let ph = self.ph();
        self.node.charge_balance_with(ph, &self.strong)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn overtitrated_book_saturates_instead_of_panicking() {
        // Дебаг ревизии V12: перетитрованная книга роняла ВЕСЬ движок
        // (`.expect()` в `ph()`), хотя §24 требует от этого состояния быть
        // отчётом, а не падением. Достижимо без экзотики: маленький
        // резервуар при штатном токе выжигается за минуты.
        let mut anode = ReservoirBook::fresh(&crate::medium::TAE_MEDIUM, 0.002);
        let mut cathode = ReservoirBook::fresh(&crate::medium::TAE_MEDIUM, 0.002);
        // Заряд заведомо за скобкой pH [0; 14].
        anode.add_anodic_charge(2000.0);
        cathode.add_cathodic_charge(2000.0);
        let (pa, pc) = (anode.ph(), cathode.ph());
        assert!(pa.is_finite() && pc.is_finite(), "pH книги не конечен");
        assert!((0.0..=14.0).contains(&pa), "pH анода вне скобки: {pa}");
        assert!((0.0..=14.0).contains(&pc), "pH катода вне скобки: {pc}");
        assert!(pa < pc, "анод обязан быть кислее катода: {pa} / {pc}");
        // На штатном объёме поведение прежнее — корень настоящий.
        let mut normal = ReservoirBook::fresh(&crate::medium::TAE_MEDIUM, 0.150);
        normal.add_anodic_charge(100.0);
        let ph = normal.ph();
        assert!(
            (5.0..9.0).contains(&ph),
            "штатная книга обязана давать НАСТОЯЩИЙ корень, а не край: {ph}"
        );
    }

    use super::*;
    use crate::medium::{TAE_MEDIUM, TBE_MEDIUM};

    /// Минуты до порога pH при токе `i_a` и объёме книги `v_l` — шагами
    /// по секунде через ФАРАДЕЕВСКИЙ путь (не замкнутой формулой): тест
    /// проверяет ту же арифметику, что работает в петле.
    fn minutes_to_ph(spec: &MediumSpec, v_l: f64, i_a: f64, target_ph: f64, anodic: bool) -> f64 {
        let mut b = ReservoirBook::fresh(spec, v_l);
        let dt = 1.0_f64;
        let mut t = 0.0;
        loop {
            let q = i_a * dt;
            if anodic {
                b.add_anodic_charge(q);
            } else {
                b.add_cathodic_charge(q);
            }
            t += dt;
            let ph = b.ph();
            if (anodic && ph <= target_ph) || (!anodic && ph >= target_ph) {
                return t / 60.0;
            }
            assert!(t < 3.0e4, "порог pH не достигнут за 500 мин");
        }
    }

    #[test]
    fn anchor_8_tae_book_sours_at_41_and_51_minutes() {
        // №8 (M3-PLAN §6): TAE 1×, 150 мл, 108 мА. Независимое
        // подтверждение чисел — bookref.py: 41.0 и 51.4 мин.
        let ta = minutes_to_ph(&TAE_MEDIUM, 0.150, 0.108, 6.0, true);
        let tc = minutes_to_ph(&TAE_MEDIUM, 0.150, 0.108, 10.0, false);
        assert!(
            (ta - 41.0).abs() <= 4.1,
            "анод TAE до pH 6 за {ta:.1} мин против якоря 41 ± 10%"
        );
        assert!(
            (tc - 51.0).abs() <= 5.1,
            "катод TAE до pH 10 за {tc:.1} мин против якоря 51 ± 10%"
        );
    }

    #[test]
    fn anchor_9_tbe_lasts_4_5x_longer() {
        // №9: отношение времён устойчивее абсолютов (объём сокращается).
        // На меньшем объёме те же τ*: t = τ·F·V/I — отношение не зависит
        // от V и I. Независимо: 186.2/41.0 = 4.55 (bookref.py).
        let v = 0.015; // 1/10 объёма — в 10 раз быстрее, тот же ответ
        let ta_tae = minutes_to_ph(&TAE_MEDIUM, v, 0.108, 6.0, true);
        let ta_tbe = minutes_to_ph(&TBE_MEDIUM, v, 0.108, 6.0, true);
        let r = ta_tbe / ta_tae;
        assert!(
            (3.8..=5.2).contains(&r),
            "TBE/TAE = {r:.2} вне полосы 4.5 ± 15%"
        );
        // Катодная ветка TBE — тоже дольше TAE в разы (182 против 51).
        let tc_tae = minutes_to_ph(&TAE_MEDIUM, v, 0.108, 10.0, false);
        let tc_tbe = minutes_to_ph(&TBE_MEDIUM, v, 0.108, 10.0, false);
        assert!(
            tc_tbe / tc_tae > 2.5,
            "катод TBE/TAE = {:.2}",
            tc_tbe / tc_tae
        );
    }

    #[test]
    fn anchor_15_electroneutrality_of_titrated_book() {
        // №15: инвариант электронейтральности узла с титрантом — полная
        // форма после специации сходится к нулю на машинном уровне.
        let mut b = ReservoirBook::fresh(&TAE_MEDIUM, 0.150);
        b.add_anodic_charge(0.108 * 1200.0); // 20 минут тока
        let scale = b.node.tris_m.max(b.node.anion_m);
        let res = b.electroneutrality_residual_m().abs();
        assert!(
            res < 1e-12 * scale.max(1e-3),
            "остаток электронейтральности {res:.3e}"
        );
    }

    #[test]
    fn faraday_bookkeeping_is_exact() {
        // Бухгалтерия: τ·V·F == заряд бит-в-бит (одно деление на вводе).
        let mut b = ReservoirBook::fresh(&TAE_MEDIUM, 0.5);
        let q = 12.345_f64;
        b.add_anodic_charge(q);
        assert_eq!(b.charge_c, q);
        let tau = b.node.titrant_strong_acid_m;
        assert_eq!(tau, q / (FARADAY_C_MOL * 0.5));
        // Катодный ввод того же заряда возвращает титрант в ноль точно.
        b.add_cathodic_charge(q);
        assert_eq!(b.node.titrant_strong_acid_m, 0.0);
    }
}
