//! Кислотно-основное равновесие и проводимость буферов TAE и TBE при 25 °C.
//! Расчёт по предельным эквивалентным проводимостям даёт верхнюю оценку κ₂₅.
//! Рабочее значение учитывает экранировку по ионной силе.

use crate::fmath::DetMath;

/// Ионное произведение воды при 25 °C.
pub const KW: f64 = 1e-14;

/// pKa при 25 °C (литература: CRC Handbook).
pub const PK_TRIS: f64 = 8.06;
pub const PK_ACETIC: f64 = 4.76;
pub const PK_BORIC: f64 = 9.24;
/// Четыре ступени EDTA (H₄Y → Y⁴⁻).
pub const PK_EDTA: [f64; 4] = [2.0, 2.67, 6.16, 10.26];

/// Предельные эквивалентные проводимости при 25 °C, См·см²/моль-экв (CRC Handbook).
pub mod lambda {
    pub const TRIS_H: f64 = 29.5;
    pub const ACETATE: f64 = 40.9;
    pub const BORATE: f64 = 30.5;
    pub const H: f64 = 349.8;
    pub const OH: f64 = 198.6;
}

/// Приближённая средняя эквивалентная проводимость EDTA, См·см²/моль-экв.
/// Диапазон 40–70 меняет расчётную κ буфера TAE примерно на ±2,7%.
pub const EDTA_EQ_PLACEHOLDER: f64 = 55.0;

/// Предельные эквивалентные проводимости ионов пробы при 25 °C, См·см²/моль-экв
/// (Sawyer & McCarty, Bard & Faulkner).
pub mod lambda_strong {
    pub const NA: f64 = 50.11;
    pub const CL: f64 = 76.34;
    /// Эквивалентная (на заряд): ½Mg²⁺.
    pub const MG_EQ: f64 = 53.06;
}

/// Поправка к предельной проводимости при 25 °C по закону Онзагера.
/// Релаксационный член масштабируется с квадратом заряда, электрофоретический — с его модулем.
/// Экранировка по Bahga и Santiago (2010). Нижняя граница 0,5·λ° ограничивает поправку
/// для многозарядных ионов вне области предельного закона.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IonicStrengthModel {
    /// Предельные проводимости при бесконечном разбавлении.
    Ideal,
    /// Экранированная поправка с нижней границей проводимости.
    OnsagerScreened,
}

/// λ(I): поправка одной ионной формы. `z_abs` — |заряд| формы,
/// `i_molar` — ионная сила, моль/л.
pub fn lambda_at_ionic_strength(lambda0: f64, z_abs: f64, i_molar: f64) -> f64 {
    let s = i_molar.max(0.0).sqrt();
    let screen = s / (1.0 + 1.5 * s);
    let dl = (0.2289 * z_abs * z_abs * lambda0 + 30.16 * z_abs) * screen;
    (lambda0 - dl).max(0.5 * lambda0)
}

/// Доля депротонированной формы одноосновной кислоты.
#[inline]
pub fn frac_acid(ph: f64, pka: f64) -> f64 {
    1.0 / (1.0 + 10f64.dpowf(pka - ph))
}

/// Доля протонированной формы основания ([TrisH⁺] / c_Tris).
#[inline]
pub fn frac_base_protonated(ph: f64, pka: f64) -> f64 {
    1.0 / (1.0 + 10f64.dpowf(ph - pka))
}

/// Fractions of the five protonation states of EDTA and their mean signed charge.
pub fn edta_species(ph: f64) -> ([f64; 5], f64) {
    static KA: std::sync::LazyLock<[f64; 4]> =
        std::sync::LazyLock::new(|| PK_EDTA.map(|pka| 10.0_f64.dpowf(-pka)));
    let h = 10f64.dpowf(-ph);
    let mut terms = [1.0f64; 5];
    for (i, ka) in KA.iter().enumerate() {
        terms[i + 1] = terms[i] * ka / h;
    }
    let tot: f64 = terms.iter().sum();
    let mut fr = [0.0f64; 5];
    for (f, t) in fr.iter_mut().zip(terms.iter()) {
        *f = t / tot;
    }
    let charge: f64 = fr.iter().enumerate().map(|(i, f)| -(i as f64) * f).sum();
    (fr, charge)
}

/// Аналитические концентрации (моль/л) в буфере Tris, одноосновной кислоты и EDTA.
pub struct BufferSpec {
    pub name: &'static str,
    pub tris_m: f64,
    pub anion_m: f64,
    pub edta_m: f64,
    pub pka_anion: f64,
    /// Предельная проводимость аниона, См·см²/моль-экв.
    pub lambda_anion: f64,
}

/// Состав TAE 1× эквивалентен разведению прописи 50× (242 г Tris и 57,1 мл AcOH на литр).
pub const TAE_1X: BufferSpec = BufferSpec {
    name: "TAE 1x",
    tris_m: 0.040,
    anion_m: 0.020,
    edta_m: 0.001,
    pka_anion: PK_ACETIC,
    lambda_anion: lambda::ACETATE,
};

/// Состав TBE 1× соответствует прописи 5× (54 г Tris и 27,5 г борной кислоты на литр).
pub const TBE_1X: BufferSpec = BufferSpec {
    name: "TBE 1x",
    tris_m: 0.089,
    anion_m: 0.089,
    edta_m: 0.002,
    pka_anion: PK_BORIC,
    lambda_anion: lambda::BORATE,
};

/// Источник значения проводимости при 25 °C.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KappaSource {
    /// Расчёт равновесия (верхняя оценка: λ° без поправки на ионную силу).
    EquilibriumEstimate,
    Measured,
}

/// Результат равновесного расчёта буфера при 25 °C.
#[derive(Debug, Clone, Copy)]
pub struct BufferEquilibrium {
    pub spec_name: &'static str,
    pub ph: f64,
    /// [TrisH⁺], моль/л.
    pub tris_h_m: f64,
    /// [анион⁻], моль/л.
    pub anion_m: f64,
    /// Средний заряд EDTA по модулю, экв/моль.
    pub edta_mean_charge: f64,
    /// Экранированная κ₂₅, См/м, для цепи и тепловой модели.
    pub kappa_s_m: f64,
    /// Верхняя оценка κ₂₅ по предельным λ°, См/м.
    pub kappa_limiting_s_m: f64,
    /// Ионная сила, моль/л.
    pub ionic_strength_m: f64,
    pub source: KappaSource,
}

impl BufferSpec {
    /// Равновесный состав и проводимость рассчитываются той же [`crate::medium::NodeComposition`],
    /// что задаёт поле в дорожке.
    pub fn equilibrium(&self) -> BufferEquilibrium {
        use crate::medium::{IonicStrengthModel, NodeComposition, StrongIons};
        let node = NodeComposition {
            tris_m: self.tris_m,
            anion_m: self.anion_m,
            edta_m: self.edta_m,
            pka_anion: self.pka_anion,
            lambda_anion: self.lambda_anion,
            titrant_strong_acid_m: 0.0,
        };
        let eq = node.speciate();
        BufferEquilibrium {
            spec_name: self.name,
            ph: eq.ph,
            tris_h_m: eq.tris_h_m,
            anion_m: eq.anion_m,
            edta_mean_charge: eq.edta_mean_charge,
            kappa_s_m: node.conductivity_s_m(
                &eq,
                &StrongIons::ZERO,
                IonicStrengthModel::OnsagerScreened,
            ),
            kappa_limiting_s_m: node.conductivity_s_m(
                &eq,
                &StrongIons::ZERO,
                IonicStrengthModel::Ideal,
            ),
            ionic_strength_m: node.ionic_strength_m(&eq, &StrongIons::ZERO),
            source: KappaSource::EquilibriumEstimate,
        }
    }

    /// Заменяет расчётную проводимость результатом кондуктометрии, См/м.
    pub fn measured(&self, kappa_s_m: f64) -> BufferEquilibrium {
        let eq = self.equilibrium();
        BufferEquilibrium {
            kappa_s_m,
            source: KappaSource::Measured,
            ..eq
        }
    }
}

/// Проводимость TAE 1× при 25 °C, См/м, с поправкой на ионную силу.
pub const KAPPA_TAE_S_M: f64 = 0.1369;

/// Верхняя оценка для TAE 1× по предельным λ°, См/м (PHYSICS §13).
pub const KAPPA_TAE_LIMITING_S_M: f64 = 0.166;

/// Проводимость TBE 1× при 25 °C, См/м, с поправкой на ионную силу.
pub const KAPPA_TBE_S_M: f64 = 0.1087;

/// Верхняя оценка для TBE 1× по предельным λ°, См/м (PHYSICS §13).
pub const KAPPA_TBE_LIMITING_S_M: f64 = 0.143;

#[test]
fn edta_fractions_match_direct_dissociation_constants() {
    for hundredths in 0..=1400 {
        let ph = hundredths as f64 / 100.0;
        let h = 10.0_f64.dpowf(-ph);
        let mut weights = [1.0; 5];
        for (i, pka) in PK_EDTA.iter().enumerate() {
            weights[i + 1] = weights[i] * 10.0_f64.dpowf(-pka) / h;
        }
        let total: f64 = weights.iter().sum();
        let fractions = weights.map(|value| value / total);
        let mean_charge: f64 = fractions
            .iter()
            .enumerate()
            .map(|(i, f)| -(i as f64) * f)
            .sum();
        let (observed, charge) = edta_species(ph);
        assert_eq!(observed.map(f64::to_bits), fractions.map(f64::to_bits));
        assert_eq!(charge.to_bits(), mean_charge.to_bits());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_hardcoded_kappa_literals_left() {
        // Страж этапа Э1 (ревизия V12), брат
        // `no_direct_std_transcendentals_left_in_sim`.
        //
        // Переименования константы НЕДОСТАТОЧНО: копии κ жили литералами в
        // `examples/presets_failures.rs` и `tests/m4_chamber.rs`
        // (`const KAPPA_TAE: f64 = 0.166;`), поэтому камера продолжала
        // считать по предельной проводимости даже после того, как общая
        // константа стала рабочей. Дефект пойман не тестом, а сверкой
        // печатного тока пресета с предсказанием — то есть случайно.
        // Теперь ловится здесь.
        //
        // Правило: числа 0.166 / 0.143 (предельные κ TAE/TBE) в коде
        // движка, тестов и примеров пишутся ТОЛЬКО через
        // `KAPPA_*_LIMITING_S_M`, рабочие — через `KAPPA_*_S_M`.
        let root = env!("CARGO_MANIFEST_DIR");
        let mut offenders = Vec::new();
        let mut stack: Vec<std::path::PathBuf> = ["src", "tests", "examples"]
            .iter()
            .map(|d| std::path::PathBuf::from(root).join(d))
            .collect();
        while let Some(dir) = stack.pop() {
            let Ok(entries) = std::fs::read_dir(&dir) else {
                continue;
            };
            for entry in entries {
                let path = entry.expect("запись").path();
                if path.is_dir() {
                    stack.push(path);
                    continue;
                }
                if path.extension().and_then(|e| e.to_str()) != Some("rs") {
                    continue;
                }
                // В этом файле литералы — определения самих констант.
                if path.file_name().and_then(|f| f.to_str()) == Some("buffers.rs") {
                    continue;
                }
                let text = std::fs::read_to_string(&path).expect("читается");
                for (i, line) in text.lines().enumerate() {
                    // Комментарии — и целые строки, и хвосты — легитимны:
                    // в них число ЦИТИРУЕТСЯ, а не используется.
                    let code = line.split("//").next().unwrap_or("").trim();
                    for bad in ["0.166", "0.143"] {
                        if code.contains(bad) {
                            offenders.push(format!(
                                "{}:{}: {}",
                                path.file_name().unwrap().to_string_lossy(),
                                i + 1,
                                code.trim()
                            ));
                        }
                    }
                }
            }
        }
        assert!(
            offenders.is_empty(),
            "κ прописана литералом — заведётся вторая проводимость:
{}",
            offenders.join(
                "
"
            )
        );
    }

    #[test]
    fn tae_equilibrium_matches_review_13() {
        // Страж §13 PHYSICS.md: pH 7.93, [TrisH+] ≈ 23.0 мМ, [Ac−] = 20.0 мМ.
        // Проводимостей ДВЕ, и обе названы: предельная 0.166 (число ревью
        // §3.2 и таблицы §13, верхняя оценка) и рабочая 0.1369 с
        // экранировкой по ионной силе. В расчёты идёт вторая — Э1 ревизии
        // V12; первая держится здесь, чтобы переход был виден в диффе.
        let eq = TAE_1X.equilibrium();
        assert_eq!(eq.source, KappaSource::EquilibriumEstimate);
        assert!(
            (eq.ph - 7.93).abs() < 5e-3,
            "pH(TAE) = {:.4}, эталон 7.93",
            eq.ph
        );
        assert!(
            (eq.tris_h_m * 1e3 - 23.0).abs() < 0.15,
            "[TrisH+](TAE) = {:.2} мМ, эталон 23.0",
            eq.tris_h_m * 1e3
        );
        assert!(
            (eq.anion_m * 1e3 - 20.0).abs() < 0.05,
            "[Ac-](TAE) = {:.2} мМ, эталон 20.0",
            eq.anion_m * 1e3
        );
        assert!(
            (eq.kappa_limiting_s_m - KAPPA_TAE_LIMITING_S_M).abs() < 5e-4,
            "κ_предельная(TAE) = {:.4} См/м, ревью §3.2 даёт {}",
            eq.kappa_limiting_s_m,
            KAPPA_TAE_LIMITING_S_M
        );
        assert!(
            (eq.kappa_s_m - KAPPA_TAE_S_M).abs() < 5e-4,
            "κ_рабочая(TAE) = {:.4} См/м, документ {}",
            eq.kappa_s_m,
            KAPPA_TAE_S_M
        );
        // Направление поправки: экранировка обязана СНИЖАТЬ проводимость,
        // и §13 называет масштаб «−10…−15 %».
        let drop = 1.0 - eq.kappa_s_m / eq.kappa_limiting_s_m;
        assert!(
            (0.10..=0.25).contains(&drop),
            "экранировка сняла {:.1} % — вне ожидаемого коридора",
            drop * 100.0
        );
    }

    #[test]
    fn tbe_equilibrium_matches_review_13() {
        // Страж §13: pH 8.56, [TrisH+] ≈ 21.4 мМ, [борат−] ≈ 15.4 мМ,
        // κ ≈ 0.143 См/м. Борная кислота ионизована лишь на ~17% —
        // именно поэтому κ(TBE) < κ(TAE), а не наоборот.
        let eq = TBE_1X.equilibrium();
        assert!(
            (eq.ph - 8.56).abs() < 5e-3,
            "pH(TBE) = {:.4}, эталон 8.56",
            eq.ph
        );
        assert!(
            (eq.tris_h_m * 1e3 - 21.4).abs() < 0.15,
            "[TrisH+](TBE) = {:.2} мМ, эталон 21.4",
            eq.tris_h_m * 1e3
        );
        assert!(
            (eq.anion_m * 1e3 - 15.4).abs() < 0.2,
            "[борат-](TBE) = {:.2} мМ, эталон 15.4",
            eq.anion_m * 1e3
        );
        assert!(
            (eq.kappa_limiting_s_m - KAPPA_TBE_LIMITING_S_M).abs() < 5e-4,
            "kappa_limiting(TBE) = {:.4}, review {}",
            eq.kappa_limiting_s_m,
            KAPPA_TBE_LIMITING_S_M
        );
        assert!(
            (eq.kappa_s_m - KAPPA_TBE_S_M).abs() < 5e-4,
            "κ(TBE) = {:.4} См/м, документ {}",
            eq.kappa_s_m,
            KAPPA_TBE_S_M
        );
    }

    #[test]
    fn kappa_tae_exceeds_kappa_tbe() {
        // Направление сравнения, на которое опирается якорь №3:
        // наивный расчёт «вся соль диссоциирована» давал обратное.
        assert!(KAPPA_TAE_S_M > KAPPA_TBE_S_M);
        assert!(TAE_1X.equilibrium().kappa_s_m > TBE_1X.equilibrium().kappa_s_m);
    }

    #[test]
    fn measured_source_replaces_estimate() {
        // После кондуктометрии меняется одна строка и флаг: значение
        // подменяется, происхождение видно из типа.
        let m = TAE_1X.measured(0.171);
        assert_eq!(m.source, KappaSource::Measured);
        assert!((m.kappa_s_m - 0.171).abs() < 1e-12);
        // Состав при этом остаётся равновесным.
        assert!((m.ph - TAE_1X.equilibrium().ph).abs() < 1e-12);
    }

    #[test]
    fn named_constants_track_equilibrium() {
        // Константы, на которые опирается thermal.rs, не оторваны от
        // расчёта: расхождение > допуска печати документа ловится здесь.
        for (spec, c) in [(&TAE_1X, KAPPA_TAE_S_M), (&TBE_1X, KAPPA_TBE_S_M)] {
            let k = spec.equilibrium().kappa_s_m;
            assert!(
                (k - c).abs() <= 5e-4,
                "{}: расчёт {:.4} против константы {}",
                spec.name,
                k,
                c
            );
        }
    }
}
