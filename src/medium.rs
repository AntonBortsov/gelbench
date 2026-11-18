//! Кислотно-основное равновесие буфера при заданных концентрациях семейств.
//! Перенос меняет суммарную концентрацию каждого семейства, протолиз — доли его ионных форм.
//! Равновесный pH задаёт корень баланса заряда. Локальная электронейтральность оправдана
//! нанометровой длиной Дебая и временем диэлектрической релаксации порядка 7 нс.
//!
//! В перемешанном резервуаре титрант учитывает прошедший заряд. Его противоион —
//! анион буфера на аноде и TrisH⁺ на катоде. В узлах геля скалярный титрант равен нулю.

use crate::buffers::{
    edta_species, frac_acid, frac_base_protonated, lambda, KW, PK_ACETIC, PK_BORIC, PK_TRIS,
};
use crate::fmath::DetMath;

/// Аналитические концентрации кислотно-основных семейств, моль/л.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NodeComposition {
    pub tris_m: f64,
    pub anion_m: f64,
    pub edta_m: f64,
    /// pKa слабой кислоты буфера (4.76 ацетат / 9.24 борная).
    pub pka_anion: f64,
    /// λ° аниона слабой кислоты, См·см²/моль-экв.
    pub lambda_anion: f64,
    /// Эквиваленты сильной кислоты в резервуаре, моль/л. Отрицательное значение задаёт щёлочь.
    pub titrant_strong_acid_m: f64,
}

/// Равновесие может выйти за рабочий диапазон pH при истощении буфера.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpeciationError {
    /// Титрант вывел состав за буферную ёмкость: корень вне скобок [0; 14].
    OutOfBufferCapacity { titrant_strong_acid_m: f64 },
}

/// Ниже [`NodeComposition::MIN_BETA_FOR_PH`] шум концентраций делает оценку pH неопределённой.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PhStatus {
    Defined(f64),
    Undefined,
}

/// Равновесный состав узла при 25 °C (выход [`NodeComposition::try_speciate`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NodeEquilibrium {
    pub ph: f64,
    /// [H⁺], моль/л.
    pub h_m: f64,
    /// [OH⁻], моль/л.
    pub oh_m: f64,
    /// [TrisH⁺], моль/л.
    pub tris_h_m: f64,
    /// [анион⁻] слабой кислоты, моль/л.
    pub anion_m: f64,
    /// Средний заряд EDTA по модулю, экв/моль.
    pub edta_mean_charge: f64,
    /// σ₂₅, См/м (предельные λ°, без поправки на ионную силу). Вклад
    /// контриона титранта включается только в резервуарной конвенции.
    pub kappa_s_m: f64,
}

impl NodeComposition {
    /// Баланс монотонно убывает с pH, что позволяет искать корень бисекцией.
    fn charge_balance(&self, ph: f64) -> f64 {
        let h = 10f64.dpowf(-ph);
        let oh = KW / h;
        let trish = self.tris_m * frac_base_protonated(ph, PK_TRIS);
        let an = self.anion_m * frac_acid(ph, self.pka_anion);
        let (_, z_edta) = edta_species(ph);
        trish + h - oh - an - self.titrant_strong_acid_m + z_edta * self.edta_m
    }

    /// Равновесие в интервале pH от 0 до 14.
    pub fn try_speciate(&self) -> Result<NodeEquilibrium, SpeciationError> {
        self.speciate_bracketed(0.0, 14.0)
    }

    /// Начинает бисекцию в пределах трёх единиц от предыдущего pH.
    /// Если состав резко изменился, использует полный интервал.
    pub fn try_speciate_warm(&self, ph_prev: f64) -> Result<NodeEquilibrium, SpeciationError> {
        let lo = (ph_prev - 3.0).max(0.0);
        let hi = (ph_prev + 3.0).min(14.0);
        match crate::thermal::bisect(|ph| self.charge_balance(ph), lo, hi) {
            Some(ph) => Ok(self.equilibrium_from_ph(ph)),
            None => self.speciate_bracketed(0.0, 14.0),
        }
    }

    fn speciate_bracketed(&self, lo: f64, hi: f64) -> Result<NodeEquilibrium, SpeciationError> {
        match crate::thermal::bisect(|ph| self.charge_balance(ph), lo, hi) {
            Some(ph) => Ok(self.equilibrium_from_ph(ph)),
            None => Err(SpeciationError::OutOfBufferCapacity {
                titrant_strong_acid_m: self.titrant_strong_acid_m,
            }),
        }
    }

    /// Паникует при выходе за буферную ёмкость. Для меняющегося состава используйте [`Self::try_speciate`].
    pub fn speciate(&self) -> NodeEquilibrium {
        match self.try_speciate() {
            Ok(eq) => eq,
            Err(e) => panic!(
                "специация узла: титрант {:.4e} М вне буферной ёмкости ({e:?}); \
                 политика перетитрованного узла решается на интеграции",
                self.titrant_strong_acid_m
            ),
        }
    }

    fn equilibrium_from_ph(&self, ph: f64) -> NodeEquilibrium {
        let h = 10f64.dpowf(-ph);
        let oh = KW / h;
        let trish = self.tris_m * frac_base_protonated(ph, PK_TRIS);
        let an = self.anion_m * frac_acid(ph, self.pka_anion);
        let (_, z_edta) = edta_species(ph);
        // Противоион титранта в резервуаре приходит из буфера: ацетат/борат к аноду, TrisH⁺ к катоду.
        let counter = if self.titrant_strong_acid_m > 0.0 {
            self.lambda_anion * self.titrant_strong_acid_m
        } else {
            lambda::TRIS_H * (-self.titrant_strong_acid_m)
        };
        let kappa_s_cm = lambda::TRIS_H * trish
            + self.lambda_anion * an
            + lambda::H * h
            + lambda::OH * oh
            + crate::buffers::EDTA_EQ_PLACEHOLDER * z_edta.abs() * self.edta_m
            + counter;
        NodeEquilibrium {
            ph,
            h_m: h,
            oh_m: oh,
            tris_h_m: trish,
            anion_m: an,
            edta_mean_charge: z_edta.abs(),
            kappa_s_m: kappa_s_cm / 10.0,
        }
    }

    /// При малой ёмкости шум переноса сравним с вкладом диссоциации воды в баланс заряда.
    pub const MIN_BETA_FOR_PH: f64 = 1e-5; // М/pH

    /// Возвращает `Undefined` при малой ёмкости и ошибку при выходе pH за рабочий диапазон.
    pub fn ph_status(&self) -> Result<PhStatus, SpeciationError> {
        let eq = self.try_speciate()?;
        if self.buffer_capacity(eq.ph) < Self::MIN_BETA_FOR_PH {
            Ok(PhStatus::Undefined)
        } else {
            Ok(PhStatus::Defined(eq.ph))
        }
    }

    /// Сдвиг pH по невязке заряда, пересчитанной через буферную ёмкость в точке равновесия.
    pub fn delta_ph(&self, eq: &NodeEquilibrium, charge_residual_m: f64) -> f64 {
        charge_residual_m / self.buffer_capacity(eq.ph)
    }

    /// Остаток баланса по равновесным формам с учётом противоиона титранта, моль/л.
    pub fn electroneutrality_residual(&self) -> f64 {
        let eq = self.try_speciate().expect("узел в пределах ёмкости");
        let (_, z_edta) = edta_species(eq.ph);
        (eq.tris_h_m + eq.h_m - eq.oh_m - eq.anion_m - self.titrant_strong_acid_m
            + z_edta * self.edta_m)
            .abs()
    }

    /// Буферная ёмкость Ван Слайка, моль/(л·pH): минус производная концентрации добавленной кислоты по pH.
    /// Для EDTA вклад определяется дисперсией заряда ионных форм.
    /// Округлённое ln(10) сохраняет совместимость с эталонными значениями ёмкости.
    #[allow(clippy::approx_constant)]
    pub fn buffer_capacity(&self, ph: f64) -> f64 {
        let h = 10f64.dpowf(-ph);
        let f_tris = frac_base_protonated(ph, PK_TRIS);
        let f_an = frac_acid(ph, self.pka_anion);
        let (fr, z_edta) = edta_species(ph);
        let sum_i2: f64 = fr.iter().enumerate().map(|(i, f)| (i * i) as f64 * f).sum();
        let var_edta = sum_i2 - z_edta * z_edta;
        2.302585
            * (self.tris_m * f_tris * (1.0 - f_tris)
                + self.anion_m * f_an * (1.0 - f_an)
                + self.edta_m * var_edta
                + h
                + KW / h)
    }

    /// Эквиваленты сильной кислоты (моль/л) для достижения `target_ph`. Отрицательный результат означает щёлочь.
    pub fn titrate_to(&self, target_ph: f64) -> f64 {
        let f = |s: f64| {
            let shifted = NodeComposition {
                titrant_strong_acid_m: self.titrant_strong_acid_m + s,
                ..*self
            };
            match shifted.try_speciate() {
                Ok(eq) => eq.ph - target_ph,
                Err(_) => {
                    if s > 0.0 {
                        -target_ph
                    } else {
                        14.0 - target_ph
                    }
                }
            }
        };
        crate::thermal::bisect(f, -1.0, 1.0)
            .expect("титрование: pH монотонен по s, скобки [−1;+1] М накрывают корень")
    }
}

/// Рецептурный состав буфера.
#[derive(Debug, Clone, Copy)]
pub struct MediumSpec {
    pub name: &'static str,
    pub tris_m: f64,
    pub anion_m: f64,
    pub edta_m: f64,
    pub pka_anion: f64,
    pub lambda_anion: f64,
}

pub const TAE_MEDIUM: MediumSpec = MediumSpec {
    name: "TAE 1x",
    tris_m: 0.040,
    anion_m: 0.020,
    edta_m: 0.001,
    pka_anion: PK_ACETIC,
    lambda_anion: lambda::ACETATE,
};

pub const TBE_MEDIUM: MediumSpec = MediumSpec {
    name: "TBE 1x",
    tris_m: 0.089,
    anion_m: 0.089,
    edta_m: 0.002,
    pka_anion: PK_BORIC,
    lambda_anion: lambda::BORATE,
};

/// Загрузочный буфер, 10 мМ Tris и 1 мМ EDTA. Его ёмкость около 5,8·10⁻³ М/pH
/// против 4,6·10⁻⁷ М/pH у чистой воды.
pub const LOADING_BUFFER: MediumSpec = MediumSpec {
    name: "Загрузочный 10 мМ Tris + 1 мМ EDTA",
    tris_m: 0.010,
    anion_m: 0.0,
    edta_m: 0.001,
    pka_anion: PK_ACETIC,
    lambda_anion: lambda::ACETATE,
};

impl MediumSpec {
    /// Свежий (неистощённый) узел среды: нулевой титрант.
    pub fn fresh_node(&self) -> NodeComposition {
        NodeComposition {
            tris_m: self.tris_m,
            anion_m: self.anion_m,
            edta_m: self.edta_m,
            pka_anion: self.pka_anion,
            lambda_anion: self.lambda_anion,
            titrant_strong_acid_m: 0.0,
        }
    }
}

// В одномерной модели ток постоянен вдоль дорожки. Диффузионный член в замыкании
// поля подавляет ложные пики pH на границах зон (Thormann и Mosher, PMC9292373).

/// Ионная форма с зарядом, концентрацией (моль/л) и эквивалентной проводимостью (См·см²/моль-экв).
#[derive(Debug, Clone, Copy)]
pub struct IonForm {
    pub z: f64,
    pub c_m: f64,
    pub lambda: f64,
}

/// Газовая постоянная, Дж/(моль·К).
const R_GAS: f64 = 8.31446261815324;
/// Постоянная Фарадея, Кл/моль (CODATA).
pub const FARADAY: f64 = 96485.33212;

/// Диффузия одновалентного иона по Нернсту–Эйнштейну. Проводимость вводится в См·см²/моль, результат в м²/с.
#[inline]
pub fn diffusivity_m2_s(lambda_s_cm2_mol: f64, t_kelvin: f64) -> f64 {
    lambda_s_cm2_mol * 1.0e-4 * R_GAS * t_kelvin / (FARADAY * FARADAY)
}

/// Tracer diffusivity (m²/s) from equivalent conductance (S cm² per mole of unit charge).
/// At equal equivalent conductance, a divalent ion has half the monovalent diffusivity.
#[inline]
pub fn ion_diffusivity_m2_s(lambda_eq: f64, abs_charge: f64, temperature_k: f64) -> f64 {
    debug_assert!(abs_charge > 0.0);
    diffusivity_m2_s(lambda_eq, temperature_k) / abs_charge
}

impl NodeComposition {
    /// Равновесные формы для диффузионного тока, включая противоион титранта в резервуаре.
    pub fn ion_forms(&self) -> Vec<IonForm> {
        let eq = self.speciate();
        let (fr, _) = edta_species(eq.ph);
        let mut out = vec![
            IonForm {
                z: 1.0,
                c_m: eq.tris_h_m,
                lambda: lambda::TRIS_H,
            },
            IonForm {
                z: -1.0,
                c_m: eq.anion_m,
                lambda: self.lambda_anion,
            },
            IonForm {
                z: 1.0,
                c_m: eq.h_m,
                lambda: lambda::H,
            },
            IonForm {
                z: -1.0,
                c_m: eq.oh_m,
                lambda: lambda::OH,
            },
        ];
        for (i, f) in fr.iter().enumerate() {
            if *f > 0.0 && i > 0 {
                out.push(IonForm {
                    z: -(i as f64),
                    c_m: self.edta_m * f,
                    lambda: crate::buffers::EDTA_EQ_PLACEHOLDER,
                });
            }
        }
        if self.titrant_strong_acid_m > 0.0 {
            out.push(IonForm {
                z: -1.0,
                c_m: self.titrant_strong_acid_m,
                lambda: self.lambda_anion,
            });
        } else if self.titrant_strong_acid_m < 0.0 {
            out.push(IonForm {
                z: 1.0,
                c_m: -self.titrant_strong_acid_m,
                lambda: lambda::TRIS_H,
            });
        }
        out
    }

    /// Катионная и анионная доли локальной проводимости по предельным λ°.
    pub fn transport_numbers(&self) -> (f64, f64) {
        let mut s_plus = 0.0;
        let mut s_minus = 0.0;
        for f in self.ion_forms() {
            let w = f.lambda * f.z.abs() * f.c_m;
            if f.z > 0.0 {
                s_plus += w;
            } else if f.z < 0.0 {
                s_minus += w;
            }
        }
        let total = s_plus + s_minus;
        (s_plus / total, s_minus / total)
    }
}

/// Диффузионный ток одной формы, А/м². Концентрация задана в моль/л.
/// Градиент центральный внутри сетки и односторонний на границах.
pub fn diffusion_current_term(
    species_c_m: &[f64],
    z: f64,
    lambda_s_cm2_mol: f64,
    t_kelvin: f64,
    dx_m: f64,
) -> Vec<f64> {
    let n = species_c_m.len();
    if z == 0.0 {
        return vec![0.0; n];
    }
    let d = ion_diffusivity_m2_s(lambda_s_cm2_mol, z.abs(), t_kelvin);
    // моль/л → моль/м³ (×1000) в производной.
    let grad = |i: usize| -> f64 {
        let dc = if n >= 3 && i > 0 && i < n - 1 {
            (species_c_m[i + 1] - species_c_m[i - 1]) / (2.0 * dx_m)
        } else if i == 0 && n >= 2 {
            (species_c_m[1] - species_c_m[0]) / dx_m
        } else if n >= 2 {
            (species_c_m[n - 1] - species_c_m[n - 2]) / dx_m
        } else {
            0.0
        };
        dc * 1000.0
    };
    (0..n).map(|i| FARADAY * z * d * grad(i)).collect()
}

/// Проводимость (См/м) и диффузионный ток (А/м²) полностью диссоциированного
/// электролита 1:1 с одинаковым профилем концентраций катиона и аниона.
pub fn strong_electrolyte_terms(
    c_profile_m: &[f64],
    lambda_cation: f64,
    lambda_anion: f64,
    t_kelvin: f64,
    dx_m: f64,
) -> (Vec<f64>, Vec<f64>) {
    let sigma: Vec<f64> = c_profile_m
        .iter()
        .map(|&c| (lambda_cation + lambda_anion) * c / 10.0)
        .collect();
    let s_cat = diffusion_current_term(c_profile_m, 1.0, lambda_cation, t_kelvin, dx_m);
    let s_an = diffusion_current_term(c_profile_m, -1.0, lambda_anion, t_kelvin, dx_m);
    let s: Vec<f64> = s_cat.iter().zip(s_an.iter()).map(|(a, b)| a + b).collect();
    (sigma, s)
}

/// ЭДС жидкостного соединения при нулевом токе. Интеграл −S/σ соответствует знаку
/// потенциала Хендерсона для отношения концентраций c₁/c₂.
pub fn junction_potential_v(sigma_s_m: &[f64], diff_term_a_m2: &[f64], dx_m: f64) -> f64 {
    assert_eq!(sigma_s_m.len(), diff_term_a_m2.len());
    assert!(dx_m > 0.0);
    let n = sigma_s_m.len();
    if n < 2 {
        return 0.0;
    }
    let mut u = 0.0;
    for i in 0..n - 1 {
        let g0 = -diff_term_a_m2[i] / sigma_s_m[i];
        let g1 = -diff_term_a_m2[i + 1] / sigma_s_m[i + 1];
        u += 0.5 * (g0 + g1) * dx_m;
    }
    u
}

/// Поле в В/м. Координата x направлена к аноду, положительный ток течёт вдоль x.
/// Для вычисления поля проводимость должна быть положительной в каждом узле.
pub fn field_profile_e(j_a_m2: f64, sigma_s_m: &[f64], diff_term_a_m2: &[f64]) -> Vec<f64> {
    assert_eq!(
        sigma_s_m.len(),
        diff_term_a_m2.len(),
        "профили σ и диффузионного тока разной длины: {} против {}",
        sigma_s_m.len(),
        diff_term_a_m2.len()
    );
    sigma_s_m
        .iter()
        .zip(diff_term_a_m2.iter())
        .map(|(s, d)| {
            assert!(*s > 0.0, "σ(x) ≤ 0 в узле: обеднение среды, шаг недопустим");
            (j_a_m2 + d) / s
        })
        .collect()
}

/// Ионы пробы, не участвующие в протолизе, моль/л.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct StrongIons {
    pub na_m: f64,
    pub cl_m: f64,
    pub mg_m: f64,
}

impl StrongIons {
    pub const ZERO: StrongIons = StrongIons {
        na_m: 0.0,
        cl_m: 0.0,
        mg_m: 0.0,
    };
    /// Суммарный заряд, моль экв/л (входит в баланс заряда узла).
    pub fn charge_m(&self) -> f64 {
        self.na_m + 2.0 * self.mg_m - self.cl_m
    }
}

pub use crate::buffers::{lambda_at_ionic_strength, IonicStrengthModel};

impl NodeComposition {
    /// Баланс заряда с учётом сильных ионов пробы, моль экв/л.
    pub fn charge_balance_with(&self, ph: f64, s: &StrongIons) -> f64 {
        self.charge_balance(ph) + s.charge_m()
    }

    /// Равновесие с сильными ионами пробы, начиная с предыдущего pH.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    pub fn try_speciate_warm_with(
        &self,
        s: &StrongIons,
        ph_prev: f64,
    ) -> Result<NodeEquilibrium, SpeciationError> {
        // Производная баланса заряда для шага Ньютона равна минус буферной ёмкости Ван Слайка.
        let scale = (self.tris_m + self.anion_m + s.na_m + s.cl_m).max(1e-3);
        let mut ph = ph_prev.clamp(0.0, 14.0);
        for _ in 0..16 {
            let q = self.charge_balance_with(ph, s);
            if q.abs() <= 1e-15 * scale {
                return Ok(self.equilibrium_from_ph(ph));
            }
            let beta = self.buffer_capacity(ph);
            // Неположительная ёмкость и NaN требуют возврата к бисекции.
            if !(beta > 0.0) {
                break;
            }
            let step = (q / beta).clamp(-1.5, 1.5);
            let next = (ph + step).clamp(0.0, 14.0);
            if (next - ph).abs() < 1e-14 {
                ph = next;
                break;
            }
            ph = next;
        }
        if self.charge_balance_with(ph, s).abs() <= 1e-12 * scale {
            return Ok(self.equilibrium_from_ph(ph));
        }
        // Запасной путь — бисекция (узкая скобка, затем полная).
        let lo = (ph_prev - 3.0).max(0.0);
        let hi = (ph_prev + 3.0).min(14.0);
        match crate::thermal::bisect(|ph| self.charge_balance_with(ph, s), lo, hi) {
            Some(ph) => Ok(self.equilibrium_from_ph(ph)),
            None => match crate::thermal::bisect(|ph| self.charge_balance_with(ph, s), 0.0, 14.0) {
                Some(ph) => Ok(self.equilibrium_from_ph(ph)),
                None => Err(SpeciationError::OutOfBufferCapacity {
                    titrant_strong_acid_m: self.titrant_strong_acid_m,
                }),
            },
        }
    }

    /// Ионная сила, моль/л. Для EDTA используется квадрат среднего заряда:
    /// приближение пренебрегает дисперсией, когда при pH 7–9 доминирует одна форма.
    pub fn ionic_strength_m(&self, eq: &NodeEquilibrium, s: &StrongIons) -> f64 {
        let z = eq.edta_mean_charge;
        0.5 * (eq.tris_h_m
            + eq.h_m
            + eq.oh_m
            + eq.anion_m
            + s.na_m
            + s.cl_m
            + 4.0 * s.mg_m
            + z * z * self.edta_m)
    }

    /// Проводимость узла геля, См/м, с поправкой на ионную силу для всех семейств и сильных ионов.
    pub fn conductivity_s_m(
        &self,
        eq: &NodeEquilibrium,
        s: &StrongIons,
        model: IonicStrengthModel,
    ) -> f64 {
        use crate::buffers::{lambda_strong, EDTA_EQ_PLACEHOLDER};
        let i_m = self.ionic_strength_m(eq, s);
        let f = |l0: f64, z: f64| match model {
            IonicStrengthModel::Ideal => l0,
            IonicStrengthModel::OnsagerScreened => lambda_at_ionic_strength(l0, z, i_m),
        };
        let z_ed = eq.edta_mean_charge;
        let kappa_s_cm = f(lambda::TRIS_H, 1.0) * eq.tris_h_m
            + f(self.lambda_anion, 1.0) * eq.anion_m
            + f(lambda::H, 1.0) * eq.h_m
            + f(lambda::OH, 1.0) * eq.oh_m
            + f(EDTA_EQ_PLACEHOLDER, z_ed) * z_ed * self.edta_m
            + f(lambda_strong::NA, 1.0) * s.na_m
            + f(lambda_strong::CL, 1.0) * s.cl_m
            + f(lambda_strong::MG_EQ, 2.0) * 2.0 * s.mg_m;
        kappa_s_cm / 10.0
    }

    /// Плотность диффузионного тока узла S = F·Σ z·D·c, А/м
    /// (∂S/∂x — член поля E = (j + ∂ₓS)/σ). D — из λ той же модели.
    pub fn s_density_a_m(
        &self,
        eq: &NodeEquilibrium,
        s: &StrongIons,
        model: IonicStrengthModel,
        t_kelvin: f64,
    ) -> f64 {
        use crate::buffers::{lambda_strong, EDTA_EQ_PLACEHOLDER};
        let i_m = self.ionic_strength_m(eq, s);
        let f = |l0: f64, z: f64| match model {
            IonicStrengthModel::Ideal => l0,
            IonicStrengthModel::OnsagerScreened => lambda_at_ionic_strength(l0, z, i_m),
        };
        let d = |l: f64| diffusivity_m2_s(l, t_kelvin);
        let z_ed = eq.edta_mean_charge;
        FARADAY
            * 1000.0
            * (d(f(lambda::TRIS_H, 1.0)) * eq.tris_h_m + d(f(lambda::H, 1.0)) * eq.h_m
                - d(f(lambda::OH, 1.0)) * eq.oh_m
                - d(f(self.lambda_anion, 1.0)) * eq.anion_m
                - d(f(EDTA_EQ_PLACEHOLDER, z_ed)) * z_ed * self.edta_m
                + d(f(lambda_strong::NA, 1.0)) * s.na_m
                - d(f(lambda_strong::CL, 1.0)) * s.cl_m
                + ion_diffusivity_m2_s(f(lambda_strong::MG_EQ, 2.0), 2.0, t_kelvin) * 2.0 * s.mg_m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fresh_nodes_match_buffers_equilibrium() {
        // Связь вех: специация свежего узла обязана совпасть с равновесием
        // M2 (buffers.rs): pH(TAE)=7.93, pH(TBE)=8.56 — один решатель.
        let tae = TAE_MEDIUM.fresh_node().speciate();
        let tbe = TBE_MEDIUM.fresh_node().speciate();
        assert!((tae.ph - 7.93).abs() < 5e-3, "pH(TAE) = {:.4}", tae.ph);
        assert!((tbe.ph - 8.56).abs() < 5e-3, "pH(TBE) = {:.4}", tbe.ph);
        // Предельная κ узла — те же 0.166 / 0.143 См/м из §13 PHYSICS.md.
        assert!(
            (tae.kappa_s_m - crate::buffers::KAPPA_TAE_LIMITING_S_M).abs() < 5e-4,
            "σ_предельная(TAE) = {:.4}",
            tae.kappa_s_m
        );
        assert!(
            (tbe.kappa_s_m - crate::buffers::KAPPA_TBE_LIMITING_S_M).abs() < 5e-4,
            "σ_предельная(TBE) = {:.4}",
            tbe.kappa_s_m
        );
    }

    #[test]
    fn diffusion_current_dxd_term_is_measured_not_hidden() {
        // РЕВИЗИЯ V12, этап Э3 — с ИСПРАВЛЕННОЙ классификацией.
        //
        // Поле считается как E = (j + ∂ₓS)/σ. Физически диффузионный вклад
        // равен Σ F·z·D·∂ₓc, а обе реализации — и движок (`medium_loop`), и
        // независимый эталон `docs/m3-transport/gelref.py` (строка
        // `dS = np.gradient(S, h)`) — дифференцируют ПРОИЗВЕДЕНИЕ
        // S = Σ F·z·D·c. Разница есть лишний член Σ F·z·c·∂ₓD, живущий там,
        // где λ(I) зависит от координаты, то есть на солевых фронтах.
        //
        // Это НЕ расхождение движка с эталоном (они согласованы), а общее
        // приближение обеих реализаций. Поэтому здесь оно не «чинится»
        // односторонне — иначе рассыплется паритет ворот E, — а МЕРЯЕТСЯ:
        // тот же приём, что у `reservoir_model_uncertainty_is_measured_not_hidden`
        // (§23 — «ответ есть число, а не да/нет»).
        //
        // Модель фронта: две зоны, свежий TAE слева и он же плюс солевой
        // слиток справа, скат — одна ячейка (максимально резкий градиент,
        // то есть ВЕРХНЯЯ оценка обоих членов).
        let node = TAE_MEDIUM.fresh_node();
        let salt = StrongIons {
            na_m: 71.7e-3,
            cl_m: 68.37e-3,
            mg_m: 8.3e-3,
        };
        let dx = 0.1e-3_f64;
        let t = 298.15;
        // Перечень форм повторяет `s_density_a_m` пооперационно — тот же
        // порядок и те же λ, иначе сравнивались бы разные величины.
        let ions = |st: &StrongIons| -> (f64, Vec<(f64, f64, f64)>) {
            use crate::buffers::{lambda_strong, EDTA_EQ_PLACEHOLDER};
            let eq = node.try_speciate_warm_with(st, 8.0).expect("специация");
            let i_m = node.ionic_strength_m(&eq, st);
            let z_ed = eq.edta_mean_charge;
            // (заряд, концентрация моль/л, λ°)
            let forms = vec![
                (1.0, eq.tris_h_m, lambda::TRIS_H),
                (-1.0, eq.anion_m, node.lambda_anion),
                (1.0, eq.h_m, lambda::H),
                (-1.0, eq.oh_m, lambda::OH),
                (-z_ed, node.edta_m, EDTA_EQ_PLACEHOLDER),
                (1.0, st.na_m, lambda_strong::NA),
                (-1.0, st.cl_m, lambda_strong::CL),
                (2.0, st.mg_m, lambda_strong::MG_EQ),
            ];
            (i_m, forms)
        };
        let (i_left, left) = ions(&StrongIons::ZERO);
        let (i_right, right) = ions(&salt);
        // Обе суммы считаются по ОДНОМУ перечню форм: физическая
        // Σ z·D·Δc и паразитная Σ z·c·ΔD (в тех же единицах, А/м²).
        let d_of = |lam: f64, z: f64, i_m: f64| {
            let screened = lambda_at_ionic_strength(lam, z.abs(), i_m);
            if z == 2.0 {
                ion_diffusivity_m2_s(screened, 2.0, t)
            } else {
                // Weak-acid EDTA uses a single effective family diffusivity.
                diffusivity_m2_s(screened, t)
            }
        };
        let (mut physical, mut spurious) = (0.0_f64, 0.0_f64);
        for (&(z_l, c_l, lam), &(z_r, c_r, _)) in left.iter().zip(right.iter()) {
            let d_l = d_of(lam, z_l, i_left);
            let d_r = d_of(lam, z_r, i_right);
            let z_mid = 0.5 * (z_l + z_r);
            let c_mid = 0.5 * (c_l + c_r) * 1000.0; // моль/м³
            let d_mid = 0.5 * (d_l + d_r);
            physical += FARADAY * z_mid * d_mid * (c_r - c_l) * 1000.0 / dx;
            spurious += FARADAY * z_mid * c_mid * (d_r - d_l) / dx;
        }
        let ratio = spurious.abs() / physical.abs();
        println!(
            "Э3: паразитный член ∂ₓD составляет {:.2} % физического              (I: {:.1} → {:.1} мМ)",
            ratio * 100.0,
            i_left * 1e3,
            i_right * 1e3
        );
        // Записанная величина модельной неопределённости. Окно широкое
        // намеренно: это МЕРА, а не калибровка; сужать её при уточнении
        // модели — нужно, ослаблять — нельзя.
        assert!(
            ratio < 0.25,
            "паразитный член вырос до {:.1} % — приближение ∂ₓ(ΣzDc) больше              не второго порядка, форму придётся менять В ОБЕИХ реализациях",
            ratio * 100.0
        );
    }

    #[test]
    fn one_conductivity_for_one_buffer() {
        // ГЛАВНЫЙ СТРАЖ ЭТАПА Э1 (ревизия V12). До него движок держал два
        // ответа на вопрос «сколько проводит TAE»: 0.166 в цепи, тепле и
        // камере против 0.1369 в петле среды, которая двигает ДНК. Расхождение
        // ×1.176 шло прямо в ток прибора и в джоулев нагрев.
        //
        // Здесь сверяются ДВА НЕЗАВИСИМЫХ ПУТИ одного числа: сводка буфера
        // (`buffers::BufferSpec::equilibrium`) и поузловая σ петли среды
        // (`NodeComposition::conductivity_s_m`). Допуск машинный, а не
        // «инженерный»: обе стороны обязаны считать по одной формуле
        // ([`crate::buffers::lambda_at_ionic_strength`]), и любое расхождение
        // означает, что формул снова стало две.
        for (spec, medium, name) in [
            (&crate::buffers::TAE_1X, &TAE_MEDIUM, "TAE"),
            (&crate::buffers::TBE_1X, &TBE_MEDIUM, "TBE"),
        ] {
            let node = medium.fresh_node();
            let eq = node.speciate();
            let sigma_loop =
                node.conductivity_s_m(&eq, &StrongIons::ZERO, IonicStrengthModel::OnsagerScreened);
            let kappa_book = spec.equilibrium().kappa_s_m;
            assert!(
                (sigma_loop - kappa_book).abs() < 1e-12,
                "{name}: петля среды даёт {sigma_loop:.12}, сводка буфера {kappa_book:.12}                  — формула проводимости снова раздвоилась"
            );
            // И тот же путь при ВЫКЛЮЧЕННОЙ экранировке обязан дать предельное
            // значение — то, что раньше уходило в цепь.
            let sigma_ideal =
                node.conductivity_s_m(&eq, &StrongIons::ZERO, IonicStrengthModel::Ideal);
            let kappa_lim = spec.equilibrium().kappa_limiting_s_m;
            assert!(
                (sigma_ideal - kappa_lim).abs() < 1e-12,
                "{name}: предельные пути разошлись, {sigma_ideal:.12} против {kappa_lim:.12}"
            );
        }
    }

    #[test]
    fn titration_table_matches_plan() {
        // Титрование до заданного pH (M3-PLAN §6, проверено расчётом):
        // TAE: до 7.0 — 14.0 мМ H⁺; до 6.0 — 18.3; до 4.6 — 29.8; до 10 — 23.0 мМ OH⁻.
        // TBE: до 6.0 — 83.4; до 10 — 81.6 мМ OH⁻.
        let tae = TAE_MEDIUM.fresh_node();
        let tbe = TBE_MEDIUM.fresh_node();
        for (node, target, expect_mm, label) in [
            (&tae, 7.0, 14.0, "TAE→7.0"),
            (&tae, 6.0, 18.3, "TAE→6.0"),
            (&tae, 4.6, 29.8, "TAE→4.6"),
            (&tae, 10.0, -23.0, "TAE→10"),
            (&tbe, 6.0, 83.4, "TBE→6.0"),
            (&tbe, 10.0, -81.6, "TBE→10"),
        ] {
            let got = node.titrate_to(target) * 1e3;
            assert!(
                (got - expect_mm).abs() < 0.15,
                "{label}: {got:.2} мМ против эталона {expect_mm} мМ"
            );
        }
    }

    #[test]
    fn buffer_capacity_at_8_3_matches_plan() {
        // β(8.3), мМ/pH: TAE 21.4 (Tris 21.4, ацетат 0.01, EDTA 0.04);
        // TBE 66.5 (Tris 47.5 + борат 18.9). Ацетат при pH 8.3 НЕ буферит —
        // вся ёмкость TAE это Tris; отсюда трёхкратная разница с TBE.
        let tae = TAE_MEDIUM.fresh_node();
        let tbe = TBE_MEDIUM.fresh_node();
        let b_tae = tae.buffer_capacity(8.3) * 1e3;
        let b_tbe = tbe.buffer_capacity(8.3) * 1e3;
        assert!((b_tae - 21.4).abs() < 0.5, "β(TAE, 8.3) = {b_tae:.2} мМ/pH");
        assert!((b_tbe - 66.5).abs() < 0.8, "β(TBE, 8.3) = {b_tbe:.2} мМ/pH");
        assert!(b_tbe / b_tae > 2.5, "трёхкратная разница ёмкостей потеряна");
        // Компонентный вклад Tris в β(TAE): 2.303·40·f(1−f), f=f(8.3; 8.06).
        let f = frac_base_protonated(8.3, PK_TRIS);
        #[allow(clippy::approx_constant)] // тот же ln10-литерал эталона
        let tris_contrib = 2.302585 * 0.040 * f * (1.0 - f) * 1e3;
        assert!(
            (tris_contrib - 21.4).abs() < 0.15,
            "вклад Tris в β(TAE) = {tris_contrib:.2}, эталон 21.4"
        );
    }

    #[test]
    fn electroneutrality_holds_to_1e12_everywhere() {
        // Якорь №15 (groundwork): полная сумма зарядов решённого узла — до
        // 1e-12 относительных. Настоящая проверка электронейтральности
        // (впрыск электрода + миграция = 0) заводится после интеграции
        // переноса — в изоляции узла она была бы проверкой конвенции,
        // а не физики.
        for node in [
            TAE_MEDIUM.fresh_node(),
            TBE_MEDIUM.fresh_node(),
            NodeComposition {
                titrant_strong_acid_m: 0.015,
                ..TAE_MEDIUM.fresh_node()
            },
            NodeComposition {
                titrant_strong_acid_m: -0.02,
                ..TBE_MEDIUM.fresh_node()
            },
        ] {
            let eq = node.speciate();
            let scale = eq
                .tris_h_m
                .max(eq.anion_m)
                .max(eq.h_m)
                .max(eq.oh_m)
                .max(1e-9);
            let res = node.electroneutrality_residual();
            assert!(
                res / scale < 1e-12,
                "узел Tris {:.1} мМ: остаток {:.3e} против масштаба {:.3e}",
                node.tris_m * 1e3,
                res,
                scale
            );
        }
    }

    #[test]
    fn diffusivity_of_proton_matches_literature() {
        // Эйнштейн–Нернст: D(H⁺) = λ·R·T/F² = 9.31e-9 м²/с — литературные
        // 9.3e-9. Проверка пересчёта единиц λ(См·см²/моль) → D(м²/с).
        let d = diffusivity_m2_s(lambda::H, 298.15);
        assert!((d - 9.31e-9).abs() / 9.31e-9 < 0.01, "D(H+) = {d:.3e} м²/с");
    }

    #[test]
    fn two_zone_conductivity_splits_field_as_inverse_sigma() {
        // Двухзонная среда без градиентов состава (S = 0): E_i = j/σ_i,
        // отношение полей обратно отношению проводимостей — последовательные
        // сопротивления. База для якоря №14 (соль пробы): солёная зона имеет
        // МЕНЬШЕЕ поле при том же токе.
        let n = 40;
        let sigma_hi = 0.30; // См/м (солёная зона)
        let sigma_lo = 0.15;
        let mut sigma = vec![sigma_lo; n];
        for s in sigma.iter_mut().take(n / 2) {
            *s = sigma_hi;
        }
        let zero = vec![0.0; n];
        let e = field_profile_e(100.0, &sigma, &zero); // j = 100 А/м²
        assert!((e[5] - 100.0 / sigma_hi).abs() < 1e-9);
        assert!((e[35] - 100.0 / sigma_lo).abs() < 1e-9);
        assert!((e[5] / e[35] - sigma_lo / sigma_hi).abs() < 1e-9);
    }

    #[test]
    fn diffusion_term_matches_analytic_gaussian_derivative() {
        // S(x) = F·z·D·∂c/∂x: численная производная против аналитики
        // для гаусса c(x) = c₀·exp(−x²/2s²), вдали от краёв сетки.
        let n = 201;
        let dx = 1e-5; // м
        let x0 = 100.0_f64 * dx;
        // Ширина 20 ячеек: центральная разница имеет ошибку O((dx/s)²·x²/s²),
        // при |x| ≤ 2s это < 0.2% — допуск ниже с запасом.
        let width = 20.0 * dx;
        let c0 = 0.05; // моль/л
        let c: Vec<f64> = (0..n)
            .map(|i| {
                let d = i as f64 * dx - x0;
                c0 * (-d * d / (2.0 * width * width)).dexp()
            })
            .collect();
        let z = -1.0;
        let lam = lambda::ACETATE;
        let t = 298.15;
        let s = diffusion_current_term(&c, z, lam, t, dx);
        let d = diffusivity_m2_s(lam, t);
        for i in [60, 80, 100, 120, 140] {
            let x = i as f64 * dx - x0;
            let dc_dx = c[i] * (-x / (width * width)) * 1000.0; // моль/л→моль/м³
            let expect = FARADAY * z * d * dc_dx;
            assert!(
                (s[i] - expect).abs() / expect.abs().max(1e-30) < 5e-3,
                "i={i}: S={:.6e} против аналитики {:.6e}",
                s[i],
                expect
            );
        }
    }

    #[test]
    fn extreme_titrant_stays_solvable_and_loud() {
        // При полном балансе сильный ион сам доминирует на краях pH:
        // титрант +30 мМ → pH ≈ 4.58 (переход БФБ, сверка со скриптом),
        // −25 мМ → pH > 11; даже ±100 мМ решается (H⁺/OH⁻ титранта
        // заполняют скобки) — OutOfBufferCapacity через титрант практически
        // недостижим и остаётся защитным путём для вырожденных входов.
        let tae = TAE_MEDIUM.fresh_node();
        let acid = NodeComposition {
            titrant_strong_acid_m: 0.030,
            ..tae
        }
        .speciate()
        .ph;
        assert!((acid - 4.58).abs() < 0.05, "pH при +30 мМ = {acid}");
        let base = NodeComposition {
            titrant_strong_acid_m: -0.025,
            ..tae
        }
        .speciate()
        .ph;
        assert!(base > 11.0, "pH при −25 мМ = {base}");

        let very_acid = NodeComposition {
            titrant_strong_acid_m: 0.100,
            ..tae
        }
        .try_speciate();
        assert!(
            very_acid.as_ref().expect("решается").ph < 2.0,
            "pH при +100 мМ = {:.3}",
            very_acid.unwrap().ph
        );
        let very_base = NodeComposition {
            titrant_strong_acid_m: -0.100,
            ..tae
        }
        .try_speciate();
        assert!(
            very_base.as_ref().expect("решается").ph > 12.0,
            "pH при −100 мМ = {:.3}",
            very_base.unwrap().ph
        );
    }

    #[test]
    fn degenerate_pure_water_node_speciates_neutral() {
        // Вырожденный узел (все семейства нулевые): баланс h − oh = 0 → pH 7,
        // κ > 0 (вода проводит чуть-чуть) — специация не делит на ноль.
        let water = NodeComposition {
            tris_m: 0.0,
            anion_m: 0.0,
            edta_m: 0.0,
            pka_anion: 4.76,
            lambda_anion: 40.9,
            titrant_strong_acid_m: 0.0,
        };
        let eq = water.speciate();
        assert!((eq.ph - 7.0).abs() < 1e-6, "pH воды = {:.6}", eq.ph);
        assert!(eq.kappa_s_m > 0.0);
    }

    #[test]
    fn warm_start_identity_with_cold_start() {
        // Тёплый старт — УСКОРЕНИЕ, а не другая физика: результат обязан
        // совпадать с холодным до 1e-10 по всем полям равновесия. Проверяем
        // три сценария догадки: точное ph_prev, у края узкой скобки (±2.9)
        // и за скобкой (±5 — корень не найдётся, сработает запасной путь).
        let nodes: Vec<NodeComposition> = vec![
            TAE_MEDIUM.fresh_node(),
            TBE_MEDIUM.fresh_node(),
            LOADING_BUFFER.fresh_node(),
            NodeComposition {
                titrant_strong_acid_m: 0.015,
                ..TAE_MEDIUM.fresh_node()
            },
            NodeComposition {
                titrant_strong_acid_m: -0.02,
                ..TBE_MEDIUM.fresh_node()
            },
            NodeComposition {
                tris_m: 0.0,
                anion_m: 0.0,
                edta_m: 0.0,
                pka_anion: 4.76,
                lambda_anion: 40.9,
                titrant_strong_acid_m: 0.0,
            },
        ];
        for node in &nodes {
            let cold = node.try_speciate().expect("холодный старт решается");
            for guess in [
                cold.ph,
                cold.ph + 2.9,
                cold.ph - 2.9,
                cold.ph + 5.0,
                cold.ph - 5.0,
            ] {
                let warm = node
                    .try_speciate_warm(guess)
                    .expect("тёплый старт обязан дойти до того же корня");
                let diff = (cold.ph - warm.ph).abs()
                    + (cold.h_m - warm.h_m).abs() / cold.h_m.max(1e-14)
                    + (cold.oh_m - warm.oh_m).abs() / cold.oh_m.max(1e-14)
                    + (cold.tris_h_m - warm.tris_h_m).abs() / cold.tris_h_m.max(1e-14)
                    + (cold.anion_m - warm.anion_m).abs() / cold.anion_m.max(1e-14)
                    + (cold.edta_mean_charge - warm.edta_mean_charge).abs()
                    + (cold.kappa_s_m - warm.kappa_s_m).abs();
                assert!(
                    diff < 1e-10,
                    "тёплый старт от {:.3} разошёлся с холодным: Δ = {diff:.3e}",
                    guess
                );
            }
        }
    }

    #[test]
    fn newton_speciation_identity_with_bisect() {
        // Ньютон в try_speciate_warm_with — ускорение, не другая физика:
        // pH обязан совпадать с бисекцией полного баланса до 1e-10 на
        // представительном наборе узлов и догадок.
        let salt = StrongIons {
            na_m: 0.05,
            cl_m: 0.0683,
            mg_m: 0.0083,
        };
        let cases = [
            (TAE_MEDIUM.fresh_node(), StrongIons::ZERO),
            (TBE_MEDIUM.fresh_node(), StrongIons::ZERO),
            (LOADING_BUFFER.fresh_node(), salt),
            (TAE_MEDIUM.fresh_node(), salt),
        ];
        for (node, s) in cases {
            let exact = crate::thermal::bisect(|ph| node.charge_balance_with(ph, &s), 0.0, 14.0)
                .expect("бисекция решается");
            for guess in [exact, exact + 2.5, exact - 2.5, 7.0] {
                let got = node.try_speciate_warm_with(&s, guess).unwrap().ph;
                assert!(
                    (got - exact).abs() < 1e-10,
                    "Ньютон от {guess:.2}: {got:.12} против {exact:.12}"
                );
            }
        }
    }

    #[test]
    fn strong_ions_zero_is_identity_and_neutral_salt_keeps_ph() {
        // StrongIons::ZERO обязан не менять ни pH, ни σ (Ideal);
        // электронейтральная соль (Na = Cl) не сдвигает pH вовсе.
        for spec in [TAE_MEDIUM, TBE_MEDIUM, LOADING_BUFFER] {
            let node = spec.fresh_node();
            let cold = node.try_speciate().unwrap();
            let with0 = node
                .try_speciate_warm_with(&StrongIons::ZERO, cold.ph)
                .unwrap();
            assert!((cold.ph - with0.ph).abs() < 1e-12);
            let sigma_id =
                node.conductivity_s_m(&with0, &StrongIons::ZERO, IonicStrengthModel::Ideal);
            assert!(
                (sigma_id - cold.kappa_s_m).abs() < 1e-12,
                "σ(Ideal, ZERO) = {sigma_id} против kappa {}",
                cold.kappa_s_m
            );
            let salt = StrongIons {
                na_m: 0.05,
                cl_m: 0.05,
                mg_m: 0.0,
            };
            let eq_salt = node.try_speciate_warm_with(&salt, cold.ph).unwrap();
            assert!(
                (eq_salt.ph - cold.ph).abs() < 1e-9,
                "нейтральная соль сдвинула pH: {} → {}",
                cold.ph,
                eq_salt.ph
            );
        }
    }

    #[test]
    fn lambda_ionic_strength_scale_and_floor() {
        // Масштаб поправки (разведка R2 §3): при I = 100 мМ падение
        // одновалентных 10–25%; f(I=0) = 1; монотонность; пол 0.5·λ⁰
        // для многозарядных (EDTA z̄≈3).
        use crate::buffers::lambda_strong;
        assert_eq!(lambda_at_ionic_strength(50.11, 1.0, 0.0), 50.11);
        let na = lambda_at_ionic_strength(lambda_strong::NA, 1.0, 0.1);
        let cl = lambda_at_ionic_strength(lambda_strong::CL, 1.0, 0.1);
        let drop_na = 1.0 - na / lambda_strong::NA;
        let drop_cl = 1.0 - cl / lambda_strong::CL;
        assert!((drop_na - 0.178).abs() < 0.01, "Na: {drop_na:.3}");
        assert!((drop_cl - 0.134).abs() < 0.01, "Cl: {drop_cl:.3}");
        let mut prev = f64::INFINITY;
        for i in [0.0, 0.005, 0.02, 0.05, 0.1, 0.2] {
            let v = lambda_at_ionic_strength(76.34, 1.0, i);
            assert!(v < prev + 1e-12, "немонотонность при I = {i}");
            prev = v;
        }
        // Пол: EDTA z=3 при 100 мМ упирается в 0.5·λ⁰.
        assert_eq!(lambda_at_ionic_strength(55.0, 3.0, 0.1), 27.5);
    }

    #[test]
    fn sigma_with_ionic_correction_matches_review_numbers() {
        // Числа сверены независимым скриптом ревью (этап 1): I(TAE) = 25.95 мМ,
        // σ_ons(TAE) = 0.1369 См/м (−17.5% от предельной); I(TBE) = 27.48 мМ,
        // σ_ons(TBE) = 0.1087 (−24.1%); TAE + 50 мМ NaCl: σ_ons = 0.6710.
        let tae = TAE_MEDIUM.fresh_node();
        let tbe = TBE_MEDIUM.fresh_node();
        let eq_tae = tae.try_speciate().unwrap();
        let eq_tbe = tbe.try_speciate().unwrap();
        let i_tae = tae.ionic_strength_m(&eq_tae, &StrongIons::ZERO) * 1e3;
        let i_tbe = tbe.ionic_strength_m(&eq_tbe, &StrongIons::ZERO) * 1e3;
        assert!((i_tae - 25.95).abs() < 0.15, "I(TAE) = {i_tae:.2} мМ");
        assert!((i_tbe - 27.48).abs() < 0.15, "I(TBE) = {i_tbe:.2} мМ");
        let s_tae = tae.conductivity_s_m(
            &eq_tae,
            &StrongIons::ZERO,
            IonicStrengthModel::OnsagerScreened,
        );
        let s_tbe = tbe.conductivity_s_m(
            &eq_tbe,
            &StrongIons::ZERO,
            IonicStrengthModel::OnsagerScreened,
        );
        assert!((s_tae - 0.1369).abs() < 0.001, "σ_ons(TAE) = {s_tae:.4}");
        assert!((s_tbe - 0.1087).abs() < 0.001, "σ_ons(TBE) = {s_tbe:.4}");
        let salt = StrongIons {
            na_m: 0.05,
            cl_m: 0.05,
            mg_m: 0.0,
        };
        let eq_s = tae.try_speciate_warm_with(&salt, eq_tae.ph).unwrap();
        let s_salt = tae.conductivity_s_m(&eq_s, &salt, IonicStrengthModel::OnsagerScreened);
        assert!(
            (s_salt - 0.6710).abs() < 0.005,
            "σ_ons(TAE+NaCl) = {s_salt:.4}"
        );
    }

    #[test]
    fn beta_table_matches_frozen_and_ph_review() {
        // Таблица β по Ван Слайку из FROZEN-AND-PH.md §2 — числа стража
        // δpH = δρ/β. β оценивается в ТОЧКАХ таблицы на свежих составах.
        let cases: [(&str, NodeComposition, f64, f64); 4] = [
            ("гель TAE 1×", TAE_MEDIUM.fresh_node(), 8.0, 2.30e-2),
            ("гель TBE 1×", TBE_MEDIUM.fresh_node(), 8.5, 6.68e-2),
            ("буфер загрузки", LOADING_BUFFER.fresh_node(), 8.0, 5.78e-3),
            (
                "вода",
                NodeComposition {
                    tris_m: 0.0,
                    anion_m: 0.0,
                    edta_m: 0.0,
                    pka_anion: 4.76,
                    lambda_anion: 40.9,
                    titrant_strong_acid_m: 0.0,
                },
                7.0,
                4.61e-7,
            ),
        ];
        for (label, node, ph, beta_ref) in cases {
            let beta = node.buffer_capacity(ph);
            let rel = (beta - beta_ref).abs() / beta_ref;
            assert!(
                rel < 0.01,
                "{label}: β({ph}) = {beta:.3e} против табличных {beta_ref:.3e} (отн. {rel:.2e})"
            );
        }
        // Пять порядков между забуференной зоной и водой — смысл метрики.
        let b_tae = TAE_MEDIUM.fresh_node().buffer_capacity(8.0);
        let b_water = cases[3].1.buffer_capacity(7.0);
        assert!(
            b_tae / b_water > 1e4,
            "контраст ёмкостей потерян: {b_tae:.2e} / {b_water:.2e}"
        );
    }

    #[test]
    fn ph_status_policy_undefined_below_min_beta() {
        // Политика «где β ниже порога — pH числом не печатать»: вода и соль
        // на воде → Undefined; гель и буфер загрузки → Defined. Перетитрованный
        // узел пробрасывает ошибку специации, а не выдаёт шумное число.
        let water = NodeComposition {
            tris_m: 0.0,
            anion_m: 0.0,
            edta_m: 0.0,
            pka_anion: 4.76,
            lambda_anion: 40.9,
            titrant_strong_acid_m: 0.0,
        };
        assert_eq!(water.ph_status().unwrap(), PhStatus::Undefined);
        // След Trиса (1 мкМ): узел формально «забуферен», но β ≈ 5.8e−7
        // ниже порога — Undefined. Соль на воде (инертный NaCl) в системе
        // семейств не выражается вовсе: сильный электролит без протолиза —
        // отдельная семья, не титрант.
        let trace_tris = NodeComposition {
            tris_m: 1e-6,
            ..water
        };
        assert_eq!(trace_tris.ph_status().unwrap(), PhStatus::Undefined);
        for node in [
            TAE_MEDIUM.fresh_node(),
            TBE_MEDIUM.fresh_node(),
            LOADING_BUFFER.fresh_node(),
        ] {
            match node.ph_status().expect("буферные узлы решаются") {
                PhStatus::Defined(ph) => assert!((ph - node.speciate().ph).abs() < 1e-12),
                PhStatus::Undefined => panic!("забуференный узел потерял определённость"),
            }
        }
        // δpH = δρ/β: 1 % невязки заряда геля TAE (0.046 М → 4.6e−4 М) стоит
        // ~0.02 pH по таблице ревью; та же невязка в воде — сотни pH.
        let tae = TAE_MEDIUM.fresh_node();
        let eq = tae.speciate();
        let d_tae = tae.delta_ph(&eq, 4.6e-4);
        assert!((d_tae - 0.020).abs() < 2e-3, "δpH(TAE) = {d_tae:.4}");
        let d_water = water.delta_ph(&water.speciate(), 4.6e-4);
        assert!(
            d_water > 100.0,
            "вода обязана показать катастрофу, а не шум"
        );
    }
}
