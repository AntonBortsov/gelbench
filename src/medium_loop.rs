//! Проводимость, поле и перенос ионных семейств в одномерных дорожках геля.
//! Специация задаёт локальные σ и диффузионный член S. Ток каждой дорожки постоянен вдоль x
//! и определяется падением напряжения на геле. Из поля на гранях вычисляются потоки
//! Шарфеттера–Гуммеля. Для ДНК поле фиксируется на начало макрошага.
//!
//! Граничные концентрации соответствуют свежему буферу. Дрейф учитывает локальную
//! ионную силу, диффузия использует предельные проводимости λ°.

use crate::buffers::{frac_acid, frac_base_protonated, lambda, lambda_strong, PK_TRIS};
use crate::fmath::DetMath;
use crate::medium::{
    diffusivity_m2_s, lambda_at_ionic_strength, IonicStrengthModel, MediumSpec, NodeComposition,
    StrongIons,
};
use crate::reservoir::ReservoirBook;
use crate::transport::sg_step_reservoir;

/// RT/F при 298.15 К, В — диффузионная длина фронта δ = (RT/F)/E.
const RT_F_V: f64 = 0.025_692_579;
/// На калибровке Кольрауша этот порог даёт ошибку плато до 0,01%.
pub const DEFAULT_REJECT_THRESHOLD: f64 = 0.05;

const F_TRIS: usize = 0;
const F_ANION: usize = 1;
const F_EDTA: usize = 2;
const F_NA: usize = 3;
const F_CL: usize = 4;
const F_MG: usize = 5;
const N_FAM: usize = 6;

/// Поузловые концентрации семейств и pH предыдущего шага.
#[derive(Debug, Clone)]
pub struct LaneMedium {
    /// `fam[k][i]` — семейство k в узле i, моль/л.
    pub fam: [Vec<f64>; N_FAM],
    /// Тёплый старт специации (pH прошлого шага).
    ph: Vec<f64>,
}

/// Состав лунки в моль/л. В центре он заменяет фоновый буфер.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WellLoad {
    pub tris_m: f64,
    pub anion_m: f64,
    pub edta_m: f64,
    pub na_m: f64,
    pub cl_m: f64,
    pub mg_m: f64,
}

impl WellLoad {
    /// Concentration and charge-balance limits apply to user-entered free ions.
    /// A neutral MgCl₂ addition contributes one Mg²⁺ and two Cl⁻ per formula unit.
    pub fn with_ions(mut self, ions: StrongIons) -> Result<Self, &'static str> {
        if [ions.na_m, ions.cl_m, ions.mg_m]
            .iter()
            .any(|c| !(0.0..=SALT_MAX_M).contains(c))
        {
            return Err("Na, Cl and free Mg concentrations must be finite and within 0 to 1 mol/L");
        }
        if ions.charge_m().abs() > SALT_IMBALANCE_MAX_M {
            return Err("well charge imbalance |Na + 2 Mg − Cl| exceeds 0.02 mol/L");
        }
        self.na_m += ions.na_m;
        self.cl_m += ions.cl_m;
        self.mg_m += ions.mg_m;
        Ok(self)
    }

    /// Лунка с составом фонового геля.
    pub fn background(spec: &MediumSpec) -> Self {
        Self {
            tris_m: spec.tris_m,
            anion_m: spec.anion_m,
            edta_m: spec.edta_m,
            na_m: 0.0,
            cl_m: 0.0,
            mg_m: 0.0,
        }
    }
}

/// Concentration limit (mol/L) and charge imbalance limit (mol charge/L).
pub const SALT_MAX_M: f64 = 1.0;
pub const SALT_IMBALANCE_MAX_M: f64 = 0.02;

/// Диагностика петли за прогон (счётчики — накопительные).
#[derive(Debug, Clone, Default)]
pub struct MediumDiag {
    /// Плотность тока по дорожкам на последнем шаге, А/м².
    pub j_lanes_a_m2: Vec<f64>,
    pub substeps: u64,
    pub rejects: u64,
    /// Масса, впрыснутая клэмпом отрицательных концентраций (моль/л·ячеек).
    pub clamped_mass: f64,
    /// Максимальный δpH = δρ/β перенесённых семейств за последний шаг.
    pub delta_ph_max: f64,
    /// Число ячеек на диффузионную длину фронта. При значении меньше единицы
    /// сетка занижает концентрационный контраст при стэкинге.
    pub front_resolution_cells: f64,
    pub dt_current_s: f64,
    /// Максимальная встроенная оценка ошибки за последний макрошаг.
    pub est_last: f64,
    /// pH книг (NaN — книги не включены).
    pub book_ph_anode: f64,
    pub book_ph_cathode: f64,
    /// Заряд, прошедший через камеру с включения книг, Кл.
    pub book_charge_c: f64,
}

/// При ошибке шага [`crate::SimState`] останавливает прогон и сохраняет [`crate::SimFault`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MediumFault {
    /// У баланса заряда нет корня в интервале pH от 0 до 14.
    Speciation { lane: usize, cell: usize },
    /// Подшаг отклонён `halvings` раз подряд.
    StepRejected { halvings: u32 },
}

impl std::fmt::Display for MediumFault {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Speciation { lane, cell } => write!(
                f,
                "среда: узел {cell} дорожки {lane} вне буферной ёмкости — \
                 баланс заряда не имеет корня в pH [0; 14]"
            ),
            Self::StepRejected { halvings } => write!(
                f,
                "среда: отбраковщик ужимал подшаг {halvings} раз подряд — \
                 ошибка структурная, шаг не считается"
            ),
        }
    }
}

impl std::error::Error for MediumFault {}

#[derive(Debug, Clone)]
pub struct MediumLoop {
    spec: MediumSpec,
    ionic: IonicStrengthModel,
    lanes: Vec<LaneMedium>,
    n: usize,
    dx_m: f64,
    /// Падение напряжения на геле, В.
    u_gel_v: f64,
    /// Номинальное поле U/L, В/м — знаменатель φ.
    e_nominal_v_m: f64,
    t_kelvin: f64,
    reject_threshold: f64,
    dt_cur: f64,
    /// φ по дорожкам с начала последнего макрошага.
    phi: Vec<Vec<f64>>,
    /// При свежих граничных концентрациях однородная дорожка остаётся в равновесии.
    /// Шаг переноса нужен только после загрузки, нарушающей однородность.
    touched: Vec<bool>,
    /// Ток однородной дорожки, А/м².
    j_fresh_a_m2: f64,
    /// Интеграл dx/σ по дорожке на начало макрошага, м²/См. Определяет сопротивление для внешней цепи.
    r_prime_lane: Vec<f64>,
    front_res_fresh: f64,
    books: Option<LoopBooks>,
    pub diag: MediumDiag,
}

/// Накопление электродного заряда у катода (x = 0) и анода (x = L).
/// Резервуары получают заряд только от принятых подшагов. Изменение их pH
/// не меняет граничные концентрации в геле.
#[derive(Debug, Clone)]
pub struct LoopBooks {
    pub anode: ReservoirBook,
    pub cathode: ReservoirBook,
    /// Площадь поперечного сечения одной дорожки, м².
    pub lane_area_m2: f64,
}

struct LaneProps {
    sigma: Vec<f64>,
    s_term: Vec<f64>,
    /// Скорости семейств НА ГРАНЯХ, м/с, при E = 1 В/м (масштабируются полем).
    u_face: [Vec<f64>; N_FAM],
    delta_ph_max: f64,
}

impl MediumLoop {
    /// Буфер дорожек и фон для [`WellLoad::background`].
    pub fn spec(&self) -> &MediumSpec {
        &self.spec
    }

    pub fn new(
        spec: MediumSpec,
        lanes: usize,
        n: usize,
        dx_m: f64,
        u_gel_v: f64,
        t_kelvin: f64,
        ionic: IonicStrengthModel,
    ) -> Self {
        assert!(n >= 8 && lanes >= 1 && dx_m > 0.0 && u_gel_v > 0.0);
        let fresh = || LaneMedium {
            fam: [
                vec![spec.tris_m; n],
                vec![spec.anion_m; n],
                vec![spec.edta_m; n],
                vec![0.0; n],
                vec![0.0; n],
                vec![0.0; n],
            ],
            ph: vec![spec.fresh_node().speciate().ph; n],
        };
        let e_nom = u_gel_v / (n as f64 * dx_m);
        let mut me = Self {
            spec,
            ionic,
            lanes: (0..lanes).map(|_| fresh()).collect(),
            n,
            dx_m,
            u_gel_v,
            e_nominal_v_m: e_nom,
            t_kelvin,
            reject_threshold: DEFAULT_REJECT_THRESHOLD,
            dt_cur: f64::INFINITY,
            phi: (0..lanes).map(|_| vec![1.0; n]).collect(),
            touched: vec![false; lanes],
            j_fresh_a_m2: 0.0,
            r_prime_lane: vec![f64::NAN; lanes],
            front_res_fresh: f64::INFINITY,
            books: None,
            diag: MediumDiag {
                j_lanes_a_m2: vec![0.0; lanes],
                dt_current_s: 0.0,
                front_resolution_cells: f64::INFINITY,
                book_ph_anode: f64::NAN,
                book_ph_cathode: f64::NAN,
                ..Default::default()
            },
        };
        let p0 = me
            .lane_props(0)
            .expect("свежий буфер рецепта специируется всегда");
        let (e0, j0, rp0) = me.field_of(&p0);
        me.j_fresh_a_m2 = j0;
        for r in me.r_prime_lane.iter_mut() {
            *r = rp0;
        }
        let e_max = e0.iter().cloned().fold(0.0, f64::max);
        me.front_res_fresh = RT_F_V / (e_max * dx_m);
        for j in me.diag.j_lanes_a_m2.iter_mut() {
            *j = j0;
        }
        me
    }

    /// Заменяет состав внутри лунки с косинусным переходом к фону у её краёв.
    /// `half_m` и `ramp_m` — полуширина лунки и ширина перехода в метрах.
    pub fn load_well(
        &mut self,
        lane: usize,
        center_m: f64,
        half_m: f64,
        ramp_m: f64,
        load: &WellLoad,
    ) {
        self.try_load_well(lane, center_m, half_m, ramp_m, load)
            .expect("load_well: состав лунки вне буферной ёмкости — используйте try_load_well")
    }

    /// При [`MediumFault::Speciation`] сохраняет прежний состав, pH и активность дорожки.
    pub fn try_load_well(
        &mut self,
        lane: usize,
        center_m: f64,
        half_m: f64,
        ramp_m: f64,
        load: &WellLoad,
    ) -> Result<(), MediumFault> {
        let n = self.n;
        let backup = self.lanes[lane].clone();
        let bg = [
            self.spec.tris_m,
            self.spec.anion_m,
            self.spec.edta_m,
            0.0,
            0.0,
            0.0,
        ];
        let tgt = [
            load.tris_m,
            load.anion_m,
            load.edta_m,
            load.na_m,
            load.cl_m,
            load.mg_m,
        ];
        let x0 = center_m - half_m;
        let x1 = center_m + half_m;
        let lm = &mut self.lanes[lane];
        for i in 0..n {
            let xc = (i as f64 + 0.5) * self.dx_m;
            // Скаты лежат внутри границ лунки. Их ширина задана в метрах и не зависит от сетки.
            let w = if xc <= x0 || xc >= x1 {
                0.0
            } else if xc < x0 + ramp_m {
                0.5 * (1.0 - (std::f64::consts::PI * (xc - x0) / ramp_m).dcos())
            } else if xc > x1 - ramp_m {
                0.5 * (1.0 + (std::f64::consts::PI * (xc - (x1 - ramp_m)) / ramp_m).dcos())
            } else {
                1.0
            };
            for k in 0..N_FAM {
                lm.fam[k][i] = bg[k] * (1.0 - w) + tgt[k] * w;
            }
        }
        for i in 0..n {
            let (node, s) = Self::node_of(&self.spec, lm, i);
            if let Ok(eq) = node.try_speciate_warm_with(&s, lm.ph[i]) {
                lm.ph[i] = eq.ph;
            }
        }
        // Пересчитываем сопротивление до первого шага, чтобы камера сразу учитывала состав пробы.
        let p = match self.lane_props(lane) {
            Ok(p) => p,
            Err(fault) => {
                self.lanes[lane] = backup;
                return Err(fault);
            }
        };
        let (_, _, r_prime) = self.field_of(&p);
        self.r_prime_lane[lane] = r_prime;
        self.touched[lane] = true;
        Ok(())
    }

    pub fn ionic_model(&self) -> IonicStrengthModel {
        self.ionic
    }

    /// Проводимость свежего буфера при 25 °C, См/м, для мостов и шунта камеры.
    pub fn fresh_kappa_s_m(&self) -> f64 {
        let node = self.spec.fresh_node();
        let eq = node.speciate();
        node.conductivity_s_m(&eq, &crate::medium::StrongIons::ZERO, self.ionic)
    }

    pub fn set_reject_threshold(&mut self, thr: f64) {
        assert!(thr > 0.0);
        self.reject_threshold = thr;
    }

    /// Два резервуара свежего буфера по `volume_l` литров. Заряд электролиза равен сумме j·A·dt по дорожкам.
    pub fn enable_books(&mut self, volume_l: f64, lane_area_m2: f64) {
        assert!(lane_area_m2 > 0.0);
        self.books = Some(LoopBooks {
            anode: ReservoirBook::fresh(&self.spec, volume_l),
            cathode: ReservoirBook::fresh(&self.spec, volume_l),
            lane_area_m2,
        });
        let b = self.books.as_mut().expect("только что создано");
        self.diag.book_ph_anode = b.anode.ph();
        self.diag.book_ph_cathode = b.cathode.ph();
        self.diag.book_charge_c = 0.0;
    }

    pub fn books(&self) -> Option<&LoopBooks> {
        self.books.as_ref()
    }

    /// Удельное сопротивление геля для внешней цепи, Ом·м.
    /// Участки одной дорожки соединены последовательно, поэтому 1/σ усредняется вдоль x.
    /// Дорожки равного сечения включены параллельно: их удельные сопротивления усредняются гармонически.
    /// Множитель [`crate::circuit::KAPPA_GEL_OVER_BUFFER`] поправляет проводимость геля
    /// только для цепи. Перенос ионов использует проводимость раствора.
    pub fn gel_resistivity_ohm_m(&self) -> f64 {
        let domain_m = self.n as f64 * self.dx_m;
        let mut sum_inv_rho = 0.0;
        for &r_prime in &self.r_prime_lane {
            assert!(
                r_prime.is_finite() && r_prime > 0.0,
                "кэш Σdx/σ не заполнен: {r_prime}"
            );
            let rho_lane = r_prime / domain_m; // ⟨1/σ⟩ дорожки, Ом·м
            sum_inv_rho += 1.0 / rho_lane;
        }
        let rho_eff = self.r_prime_lane.len() as f64 / sum_inv_rho;
        rho_eff / crate::circuit::KAPPA_GEL_OVER_BUFFER
    }

    /// Update temperature-dependent circuit resistance immediately. The next
    /// transport step uses this temperature for conductivity, drift and diffusion.
    pub fn set_temperature_k(&mut self, t_kelvin: f64) {
        assert!(t_kelvin.is_finite() && t_kelvin > 0.0);
        if self.t_kelvin == t_kelvin {
            return;
        }
        // At fixed composition all lane conductivities scale by the same factor, including inactive lanes.
        let ratio =
            crate::circuit::kappa_at(1.0, t_kelvin) / crate::circuit::kappa_at(1.0, self.t_kelvin);
        assert!(ratio > 0.0 && ratio.is_finite());
        self.j_fresh_a_m2 *= ratio;
        for resistance in &mut self.r_prime_lane {
            *resistance /= ratio;
        }
        self.t_kelvin = t_kelvin;
    }

    pub fn temperature_k(&self) -> f64 {
        self.t_kelvin
    }

    pub fn set_u_gel(&mut self, u_gel_v: f64) {
        assert!(u_gel_v >= 0.0, "U на пролёте не может быть отрицательным");
        if u_gel_v == self.u_gel_v {
            return;
        }
        if self.u_gel_v > 0.0 && u_gel_v > 0.0 {
            let k = u_gel_v / self.u_gel_v;
            self.e_nominal_v_m *= k;
            self.j_fresh_a_m2 *= k;
            self.front_res_fresh /= k;
        } else {
            // A zero current loses the prefactor required to scale it back up.
            self.e_nominal_v_m = u_gel_v / (self.n as f64 * self.dx_m);
            let kappa = crate::circuit::kappa_at(self.fresh_kappa_s_m(), self.t_kelvin);
            self.j_fresh_a_m2 = kappa * self.e_nominal_v_m;
            self.front_res_fresh = RT_F_V / (self.e_nominal_v_m * self.dx_m);
        }
        self.u_gel_v = u_gel_v;
    }

    /// Поле pH дорожки (последняя специация узлов; тёплый кэш петли).
    pub fn ph_profile(&self, lane: usize) -> &[f64] {
        &self.lanes[lane].ph
    }

    /// φ дорожки (с начала последнего макрошага).
    pub fn phi_lane(&self, lane: usize) -> &[f64] {
        &self.phi[lane]
    }

    pub fn lanes_count(&self) -> usize {
        self.lanes.len()
    }

    /// Профиль выбранного ионного семейства в моль/л.
    pub fn family(&self, lane: usize, k: usize) -> &[f64] {
        &self.lanes[lane].fam[k]
    }

    fn node_of(spec: &MediumSpec, lm: &LaneMedium, i: usize) -> (NodeComposition, StrongIons) {
        (
            NodeComposition {
                tris_m: lm.fam[F_TRIS][i],
                anion_m: lm.fam[F_ANION][i],
                edta_m: lm.fam[F_EDTA][i],
                pka_anion: spec.pka_anion,
                lambda_anion: spec.lambda_anion,
                titrant_strong_acid_m: 0.0,
            },
            StrongIons {
                na_m: lm.fam[F_NA][i],
                cl_m: lm.fam[F_CL][i],
                mg_m: lm.fam[F_MG][i],
            },
        )
    }

    /// На начало подшага: σ, S и ионные скорости на гранях при единичном поле.
    fn lane_props(&mut self, lane: usize) -> Result<LaneProps, MediumFault> {
        let n = self.n;
        let mut sigma = vec![0.0; n];
        let mut s_term = vec![0.0; n];
        let mut u_node: [Vec<f64>; N_FAM] = Default::default();
        for u in u_node.iter_mut() {
            *u = vec![0.0; n];
        }
        let mut dph_max = 0.0_f64;
        let t_scale = crate::circuit::kappa_at(1.0, self.t_kelvin);
        assert!(t_scale > 0.0, "κ(T) обнулилась: T вне области модели");
        for i in 0..n {
            let (node, s) = Self::node_of(&self.spec, &self.lanes[lane], i);
            let ph_prev = self.lanes[lane].ph[i];
            // Оцениваем сдвиг pH по невязке перенесённого состава при прежнем pH.
            let q0 = node.charge_balance_with(ph_prev, &s);
            let beta = node.buffer_capacity(ph_prev).max(1e-30);
            dph_max = dph_max.max((q0 / beta).abs());
            let eq = node
                .try_speciate_warm_with(&s, ph_prev)
                .map_err(|_| MediumFault::Speciation { lane, cell: i })?;
            self.lanes[lane].ph[i] = eq.ph;
            // Один и тот же температурный множитель нужен для σ и S в уравнении поля.
            sigma[i] = node.conductivity_s_m(&eq, &s, self.ionic) * t_scale;
            s_term[i] = node.s_density_a_m(&eq, &s, self.ionic, self.t_kelvin) * t_scale;
            let i_m = node.ionic_strength_m(&eq, &s);
            let lam_of = |l0: f64, z: f64| match self.ionic {
                IonicStrengthModel::Ideal => l0,
                IonicStrengthModel::OnsagerScreened => lambda_at_ionic_strength(l0, z, i_m),
            };
            let u_of = |l: f64| l * t_scale * 1.0e-4 / crate::medium::FARADAY;
            let a_tr = frac_base_protonated(eq.ph, PK_TRIS);
            let a_an = frac_acid(eq.ph, self.spec.pka_anion);
            u_node[F_TRIS][i] = -u_of(lam_of(lambda::TRIS_H, 1.0)) * a_tr;
            u_node[F_ANION][i] = u_of(lam_of(self.spec.lambda_anion, 1.0)) * a_an;
            u_node[F_EDTA][i] = u_of(lam_of(
                crate::buffers::EDTA_EQ_PLACEHOLDER,
                eq.edta_mean_charge,
            ));
            u_node[F_NA][i] = -u_of(lam_of(lambda_strong::NA, 1.0));
            u_node[F_CL][i] = u_of(lam_of(lambda_strong::CL, 1.0));
            u_node[F_MG][i] = -u_of(lam_of(lambda_strong::MG_EQ, 2.0));
        }
        let mut u_face: [Vec<f64>; N_FAM] = Default::default();
        for (k, uf) in u_face.iter_mut().enumerate() {
            let mut v = vec![0.0; n + 1];
            v[0] = u_node[k][0];
            v[n] = u_node[k][n - 1];
            for i in 1..n {
                v[i] = 0.5 * (u_node[k][i - 1] + u_node[k][i]);
            }
            *uf = v;
        }
        Ok(LaneProps {
            sigma,
            s_term,
            u_face,
            delta_ph_max: dph_max,
        })
    }

    /// Ток задаётся интегралом E·dx = U. Возвращает узловое поле, плотность тока и интеграл dx/σ.
    fn field_of(&self, p: &LaneProps) -> (Vec<f64>, f64, f64) {
        let n = self.n;
        let dx = self.dx_m;
        let mut ds = vec![0.0; n];
        for i in 0..n {
            ds[i] = if i == 0 {
                (p.s_term[1] - p.s_term[0]) / dx
            } else if i == n - 1 {
                (p.s_term[n - 1] - p.s_term[n - 2]) / dx
            } else {
                (p.s_term[i + 1] - p.s_term[i - 1]) / (2.0 * dx)
            };
        }
        let mut r_prime = 0.0; // Σ dx/σ
        let mut u_ds = 0.0; // Σ (∂ₓS)/σ·dx
        for i in 0..n {
            assert!(
                p.sigma[i] > 0.0,
                "σ(x) ≤ 0: обеднение среды, шаг недопустим"
            );
            r_prime += dx / p.sigma[i];
            u_ds += ds[i] / p.sigma[i] * dx;
        }
        let j = (self.u_gel_v - u_ds) / r_prime;
        let e: Vec<f64> = (0..n).map(|i| (j + ds[i]) / p.sigma[i]).collect();
        (e, j, r_prime)
    }

    /// Поток Шарфеттера–Гуммеля учитывает дрейф и диффузию на одной грани.
    /// Добавленная при отсечении отрицательных концентраций масса должна быть нулевой.
    fn move_lane(
        &self,
        lane: usize,
        p: &LaneProps,
        e: &[f64],
        dt: f64,
    ) -> ([Vec<f64>; N_FAM], f64) {
        let n = self.n;
        let mut ef = vec![0.0; n + 1];
        ef[0] = e[0];
        ef[n] = e[n - 1];
        for i in 1..n {
            ef[i] = 0.5 * (e[i - 1] + e[i]);
        }
        let res = [
            self.spec.tris_m,
            self.spec.anion_m,
            self.spec.edta_m,
            0.0,
            0.0,
            0.0,
        ];
        let t_scale = crate::circuit::kappa_at(1.0, self.t_kelvin);
        let d0 = [
            diffusivity_m2_s(lambda::TRIS_H, self.t_kelvin),
            diffusivity_m2_s(self.spec.lambda_anion, self.t_kelvin),
            diffusivity_m2_s(crate::buffers::EDTA_EQ_PLACEHOLDER, self.t_kelvin),
            diffusivity_m2_s(lambda_strong::NA, self.t_kelvin),
            diffusivity_m2_s(lambda_strong::CL, self.t_kelvin),
            crate::medium::ion_diffusivity_m2_s(lambda_strong::MG_EQ, 2.0, self.t_kelvin),
        ]
        .map(|d| d * t_scale);
        let mut out: [Vec<f64>; N_FAM] = Default::default();
        let mut clamped = 0.0;
        for k in 0..N_FAM {
            let vf: Vec<f64> = (0..=n).map(|i| p.u_face[k][i] * ef[i]).collect();
            let mut dif = sg_step_reservoir(
                &self.lanes[lane].fam[k],
                self.dx_m,
                &vf,
                d0[k],
                dt,
                res[k],
                res[k],
            );
            for v in dif.iter_mut() {
                if *v < 0.0 {
                    clamped += -*v;
                    *v = 0.0;
                }
            }
            out[k] = dif;
        }
        (out, clamped)
    }

    /// Адаптивный перенос среды. Нормированное поле для ДНК фиксируется на начало макрошага.
    /// При ошибке часть подшагов уже может быть принята. Такой прогон останавливается.
    pub fn step(&mut self, dt_macro: f64) -> Result<(), MediumFault> {
        assert!(dt_macro > 0.0);
        let lanes = self.lanes.len();
        self.dt_cur = self.dt_cur.min(dt_macro);
        let mut remaining = dt_macro;
        let mut first_substep_of_macro = true;
        let mut est_macro = 0.0_f64;
        let mut halvings_in_row = 0u32;
        let scale = self.spec.tris_m.max(self.spec.anion_m).max(0.05);
        let mut dph_step = 0.0_f64;
        let mut front_res = f64::INFINITY;

        let active: Vec<usize> = (0..lanes).filter(|&l| self.touched[l]).collect();
        for lane in 0..lanes {
            if !self.touched[lane] {
                self.diag.j_lanes_a_m2[lane] = self.j_fresh_a_m2;
            }
        }
        if active.len() < lanes {
            front_res = front_res.min(self.front_res_fresh);
        }

        while remaining > 1e-12 * dt_macro {
            let dt = self.dt_cur.min(remaining);
            let mut props = Vec::with_capacity(active.len());
            let mut fields = Vec::with_capacity(active.len());
            for &lane in &active {
                let p = self.lane_props(lane)?;
                dph_step = dph_step.max(p.delta_ph_max);
                let (e, j, r_prime) = self.field_of(&p);
                self.diag.j_lanes_a_m2[lane] = j;
                if first_substep_of_macro {
                    self.r_prime_lane[lane] = r_prime;
                }
                let e_max = e.iter().cloned().fold(0.0, f64::max);
                front_res = front_res.min(RT_F_V / (e_max * self.dx_m));
                props.push(p);
                fields.push(e);
            }
            if first_substep_of_macro {
                // Без внешнего поля ДНК не дрейфует. Диффузионный потенциал всё ещё действует на ионы.
                // При ненулевом поле деление сохраняет округление эталонного расчёта.
                let nom = self.e_nominal_v_m;
                for (a, &lane) in active.iter().enumerate() {
                    for i in 0..self.n {
                        self.phi[lane][i] = if nom > 0.0 { fields[a][i] / nom } else { 0.0 };
                    }
                }
            }
            let mut pred: Vec<[Vec<f64>; N_FAM]> = Vec::with_capacity(active.len());
            let mut clamp_p = 0.0;
            for (a, &lane) in active.iter().enumerate() {
                let (fam, cl) = self.move_lane(lane, &props[a], &fields[a], dt);
                clamp_p += cl;
                pred.push(fam);
            }
            // Поле корректора усредняется по исходному и предсказанному составу.
            // Специация предиктора временно перезаписывает pH дорожки.
            let mut est = 0.0_f64;
            let mut corr: Vec<[Vec<f64>; N_FAM]> = Vec::with_capacity(active.len());
            let mut clamp_q = 0.0;
            for (a, &lane) in active.iter().enumerate() {
                let ph_saved = self.lanes[lane].ph.clone();
                std::mem::swap(&mut self.lanes[lane].fam, &mut pred[a]);
                let props_pred = self.lane_props(lane);
                // Восстанавливаем и семейства, и pH до возможного выхода по ошибке.
                std::mem::swap(&mut self.lanes[lane].fam, &mut pred[a]);
                self.lanes[lane].ph = ph_saved;
                let p1 = props_pred?;
                let (e1, _, _) = self.field_of(&p1);
                let e_mid: Vec<f64> = fields[a]
                    .iter()
                    .zip(e1.iter())
                    .map(|(x, y)| 0.5 * (x + y))
                    .collect();
                let (fam_q, cl) = self.move_lane(lane, &props[a], &e_mid, dt);
                clamp_q += cl;
                for k in 0..N_FAM {
                    for i in 0..self.n {
                        est = est.max((fam_q[k][i] - pred[a][k][i]).abs() / scale);
                    }
                }
                corr.push(fam_q);
            }
            est_macro = est_macro.max(est);
            if est > self.reject_threshold {
                self.diag.rejects += 1;
                self.dt_cur = dt * 0.5;
                halvings_in_row += 1;
                if halvings_in_row > 60 {
                    return Err(MediumFault::StepRejected {
                        halvings: halvings_in_row,
                    });
                }
                continue;
            }
            halvings_in_row = 0;
            self.diag.clamped_mass += clamp_p + clamp_q;
            // Резервуары получают заряд только от принятого подшага.
            if let Some(b) = self.books.as_mut() {
                let j_sum: f64 = self.diag.j_lanes_a_m2.iter().sum();
                let q = j_sum * b.lane_area_m2 * dt;
                b.anode.add_anodic_charge(q);
                b.cathode.add_cathodic_charge(q);
                self.diag.book_charge_c += q;
            }
            for (a, &lane) in active.iter().enumerate() {
                self.lanes[lane].fam = std::mem::take(&mut corr[a]);
            }
            self.diag.substeps += 1;
            remaining -= dt;
            first_substep_of_macro = false;
            if est < 0.25 * self.reject_threshold && self.dt_cur < dt_macro {
                self.dt_cur = (self.dt_cur * 2.0).min(dt_macro);
            }
        }
        self.diag.dt_current_s = self.dt_cur;
        self.diag.est_last = est_macro;
        self.diag.delta_ph_max = dph_step;
        self.diag.front_resolution_cells = front_res;
        if let Some(b) = self.books.as_mut() {
            self.diag.book_ph_anode = b.anode.ph();
            self.diag.book_ph_cathode = b.cathode.ph();
        }
        Ok(())
    }
}

pub mod family {
    pub const TRIS: usize = super::F_TRIS;
    pub const ANION: usize = super::F_ANION;
    pub const EDTA: usize = super::F_EDTA;
    pub const NA: usize = super::F_NA;
    pub const CL: usize = super::F_CL;
    pub const MG: usize = super::F_MG;
}
