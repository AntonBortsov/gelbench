//! Analytic limits, grid convergence, lateral smile and thermal runaway tests.

/// Путь полосы, на котором считались пороги улыбки эталона REVIEW-M2:
/// фольклорные 2.98 см прежнего якоря №1. Здесь это ВХОД эталона —
/// как и предельная κ (см. `review_4a_config`), — поэтому литерал
/// остаётся литералом. Продуктовые числа при ДЕЙСТВУЮЩЕМ якоре
/// (1.42 см из модели V11) живут в `smile_and_red_line_at_working_physics`.
const REVIEW_PATH_5CM_CM: f64 = 2.98;

use super::*;
use crate::fmath::DetMath;

#[test]
fn thomas_parabolic_limit_h_to_zero() {
    // Предел h→0 (сток исчезает): профиль между пришпиленными краями —
    // парабола с максимумом q_v·W²/(8·k_g).
    let cfg = ThermalConfig::default();
    let s_tiny = 1e-12;
    let t = cfg.steady_profile_explicit_sink_pinned(s_tiny);
    let n = t.len();
    let q_v = cfg.q_v();
    let expected_max = q_v * cfg.width_m.powi(2) / (8.0 * K_GEL);
    let got_max = t[n / 2] - t[0];
    assert!(
        (got_max - expected_max).abs() / expected_max < 1e-4,
        "парабола: {got_max:.6} против {expected_max:.6}"
    );
}

#[test]
fn thomas_matches_cosh_solution_dirichlet_finite_h() {
    // Основной тест ревью: конечный сток, пришпиленные края — точное
    // решение ΔT(y) = (q_v·t_g/2h)·[1 − cosh(y/ℓ)/cosh(W/2ℓ)],
    // ℓ = sqrt(k_g·t_g/(2h)).
    let cfg = ThermalConfig::default(); // вода, κ=0.166, E=5, Робин по умолчанию
    let h_water = 200.0_f64; // типичный порядок для воды
    let s = 2.0 * h_water / cfg.t_gel_m;
    let t = cfg.steady_profile_explicit_sink_pinned(s);
    let n = t.len();
    let l_ell = (K_GEL * cfg.t_gel_m / (2.0 * h_water)).sqrt();
    let q_v = cfg.q_v();
    let prefactor = q_v * cfg.t_gel_m / (2.0 * h_water);
    let w_half = cfg.width_m / 2.0;
    let mut worst = 0.0f64;
    for i in 0..n {
        let y = -w_half + (i as f64) * (cfg.width_m / (n - 1) as f64);
        let exact = prefactor * (1.0 - (y / l_ell).dcosh() / (w_half / l_ell).dcosh());
        let dev = (t[i] - exact).abs() / prefactor.max(1e-9);
        worst = worst.max(dev);
    }
    // Допуск 1e-4 ограничен точностью квадратуры второго порядка
    // на сетке из 801 точки (у краёв), а не алгеброй решения:
    // вне погранслоя расхождение ≤ 2e-5.
    assert!(
        worst < 1e-4,
        "расхождение с cosh-решением {worst:.3e} (доля от масштаба перегрева)"
    );
}

#[test]
#[allow(clippy::field_reassign_with_default)] // cfg мутируется дважды
fn robin_edges_give_hotter_center_than_dirichlet() {
    // Робин (край конвектирует) держит центр горячее, чем пришпиленный
    // край, при прочих равных.
    let mut cfg = ThermalConfig::default();
    cfg.edge = EdgeKind::Pinned;
    let (t_pin, _) = cfg.steady_profile().expect("профиль Дирихле");
    cfg.edge = EdgeKind::Convective;
    let (t_conv, _) = cfg.steady_profile().expect("профиль Робина");
    let mid_pin = t_pin[t_pin.len() / 2];
    let mid_conv = t_conv[t_conv.len() / 2];
    assert!(
        mid_conv > mid_pin,
        "Робин {mid_conv:.4} должен быть горячее Дирихле {mid_pin:.4}"
    );
}

#[test]
fn smile_threshold_open_cassette_in_review_window() {
    // Ревью §3.3: открытая кассета (воздух) — порог 5.35 [5.04; 5.63] В/см.
    let cfg = ThermalConfig {
        medium: AIR,
        kappa_s_m: crate::buffers::KAPPA_TAE_LIMITING_S_M,
        cooling: CoolingModel::FreeConvection { c_coef: 0.59 },
        edge: EdgeKind::Convective,
        ..ThermalConfig::default()
    };
    let e_star = cfg
        .smile_threshold_e(REVIEW_PATH_5CM_CM)
        .expect("тепловой решатель")
        .expect("порог открытой кассеты должен существовать");
    assert!(
        (5.04..=5.63).contains(&e_star),
        "U* = {e_star:.3} вне окна ревью [5.04; 5.63]"
    );
}

// ---- Эталон §4а ревью (docs/REVIEW-M2.md): конфигурация и утилиты ----

/// Полный набор входных данных §4а без умолчаний: вода, κ=0.166 См/м,
/// E=5 В/см, W=0.100 м, t_g=0.006 м, L_c=0.08 м, C=0.59, n=6401.
/// Вход эталона REVIEW-M2 §4а. κ здесь — ПРЕДЕЛЬНАЯ 0.166 намеренно:
/// именно её брал независимый python ревью, и она тут ВХОД задачи, а
/// не рабочее значение движка. Если подставить сюда рабочую κ (Э1
/// ревизии V12), страж перестанет проверять решатель и начнёт
/// проверять выбор константы — это разные вещи. Продуктовые числа
/// при рабочей κ живут в `smile_and_red_line_at_working_physics`.
fn review_4a_config(edge: EdgeKind) -> ThermalConfig {
    ThermalConfig {
        medium: WATER,
        kappa_s_m: crate::buffers::KAPPA_TAE_LIMITING_S_M,
        e_v_per_cm: 5.0,
        h_scale: 1.0,
        // Независимый эталон ревью проверяет прежнюю корреляцию как
        // численный решатель; продуктовая Э7 проверяется отдельно.
        cooling: CoolingModel::FreeConvection { c_coef: 0.59 },
        width_m: 0.10,
        t_gel_m: 0.006,
        n_points: 6401,
        l_c_m: 0.08,
        kappa_feedback: false,
        edge,
    }
}

/// ΔT в точке |y| (мм) линейной интерполяцией профиля.
fn dt_at(cfg: &ThermalConfig, t: &[f64], y_mm: f64) -> f64 {
    let dy = cfg.width_m / (t.len() - 1) as f64;
    let pos = (y_mm * 1e-3 + cfg.width_m / 2.0) / dy;
    let i = (pos.floor() as usize).min(t.len() - 2);
    let frac = pos - i as f64;
    t[i] * (1.0 - frac) + t[i + 1] * frac
}

/// Перепад центр — крайняя дорожка (1 см от торца) — метрика A ревью.
/// Делегирует единому определению; тесты метрики A против таблицы §14
/// сторожат именно его.
fn drop_lane_a(cfg: &ThermalConfig, t: &[f64]) -> f64 {
    cfg.drop_center_to_lane(t)
}

#[test]
fn profile_matches_review_reference_4a() {
    // Эталонный профиль ΔT(y) в девяти точках, оба типа края.
    // Сверяется профиль, а не порог: порог — производная величина.
    const Y_MM: [f64; 9] = [0.0, 10.0, 20.0, 30.0, 40.0, 45.0, 48.0, 49.0, 50.0];
    const REF_PINNED: [f64; 9] = [
        0.61461, 0.61461, 0.61458, 0.61386, 0.59314, 0.49975, 0.30040, 0.17516, 0.00000,
    ];
    const REF_CONV: [f64; 9] = [
        0.61461, 0.61461, 0.61459, 0.61423, 0.60383, 0.55692, 0.45680, 0.39389, 0.30592,
    ];
    for (edge, reference) in [
        (EdgeKind::Pinned, &REF_PINNED),
        (EdgeKind::Convective, &REF_CONV),
    ] {
        let cfg = review_4a_config(edge);
        let (t, h) = cfg.steady_profile().expect("эталонный профиль");
        assert!(
            (h - 202.569).abs() < 5e-3,
            "{edge:?}: самосогласованный h = {h:.4}, эталон 202.569"
        );
        for (&y, &r) in Y_MM.iter().zip(reference.iter()) {
            let v = dt_at(&cfg, &t, y);
            assert!(
                (v - r).abs() < 1e-4,
                "{edge:?} |y|={y} мм: движок {v:.5}, эталон {r:.5}"
            );
        }
    }
}

#[test]
fn invariant_edge_never_colder_than_ambient() {
    // Инвариант 1: конвекция не может охладить край ниже среды —
    // min ΔT(y) ≥ 0. Нарушение означает ошибку в краевой строке
    // матрицы, а не в физике.
    for med in [AIR, WATER] {
        for edge in [EdgeKind::Pinned, EdgeKind::Convective] {
            for e in [3.0, 5.0, 10.0, 50.0] {
                let cfg = ThermalConfig {
                    medium: med,
                    e_v_per_cm: e,
                    edge,
                    ..ThermalConfig::default()
                };
                let (t, _) = cfg.steady_profile().expect("профиль инварианта");
                let mn = t.iter().cloned().fold(f64::INFINITY, f64::min);
                assert!(
                    mn >= 0.0,
                    "{med:?} {edge:?} E={e}: min ΔT = {mn:.3e} < 0 — краевая строка охлаждает ниже среды"
                );
            }
        }
    }
}

#[test]
fn invariant_robin_drop_le_dirichlet() {
    // Инвариант 2: пришпиленный к T_amb край — максимально возможное
    // охлаждение торца, поэтому перепад центр–дорожка (метрика A)
    // у Робина не больше, чем у Дирихле, при прочих равных.
    for med in [AIR, WATER] {
        for e in [3.0, 5.0, 10.0] {
            let base = ThermalConfig {
                medium: med,
                e_v_per_cm: e,
                ..ThermalConfig::default()
            };
            let drop = |edge: EdgeKind| {
                let cfg = ThermalConfig { edge, ..base };
                let (t, _) = cfg.steady_profile().expect("профиль сетки");
                drop_lane_a(&cfg, &t)
            };
            let d_pin = drop(EdgeKind::Pinned);
            let d_conv = drop(EdgeKind::Convective);
            assert!(
                d_conv <= d_pin,
                "{med:?} E={e}: перепад Робин {d_conv:.5} > Дирихле {d_pin:.5} — Робин охлаждает сильнее предела"
            );
        }
    }
}

#[test]
fn invariant_center_independent_of_edge_kind() {
    // Инвариант 3: при ℓ ≪ W кривизна живёт в ~3ℓ у торцов, центр —
    // на плато, поэтому ΔT(0) не зависит от типа края (эталон §4а:
    // 0.61461 в обоих случаях).
    let pin = review_4a_config(EdgeKind::Pinned);
    let conv = review_4a_config(EdgeKind::Convective);
    let (tp, _) = pin.steady_profile().expect("профиль Дирихле");
    let (tc, _) = conv.steady_profile().expect("профиль Робина");
    let mid_p = tp[tp.len() / 2];
    let mid_c = tc[tc.len() / 2];
    assert!(
        (mid_p - mid_c).abs() < 1e-4 && (mid_p - 0.61461).abs() < 1e-4,
        "ΔT(0): Дирихле {mid_p:.5}, Робин {mid_c:.5}, эталон 0.61461"
    );
}

#[test]
fn invariant_grid_convergence_one_sided_first_order() {
    // Инвариант 4: сходимость по сетке односторонняя и первого порядка —
    // перепад центр–торец (Робин) монотонно убывает к пределу:
    // 0.3584 (n=51) → 0.3115 (n=801) → 0.3087 (n=6401). Грубая сетка
    // занижает температуру края и завышает улыбку. Требование: не менее
    // 20 узлов на тепловую длину ℓ.
    let expected = [0.3584, 0.3115, 0.3087];
    let mut prev = f64::INFINITY;
    for (i, &n) in [51usize, 801, 6401].iter().enumerate() {
        let cfg = ThermalConfig {
            n_points: n,
            ..review_4a_config(EdgeKind::Convective)
        };
        let (t, h) = cfg.steady_profile().expect("эталонный профиль");
        let drop = t[t.len() / 2] - t[t.len() - 1];
        assert!(
            drop < prev,
            "n={n}: перепад {drop:.4} не убывает (пред. {prev:.4}) — сходимость обязана быть односторонней"
        );
        assert!(
            (drop - expected[i]).abs() < 2e-3,
            "n={n}: перепад {drop:.4}, эталон ревью {}",
            expected[i]
        );
        // Требование ≥20 узлов на ℓ — для рабочих сеток; грубая n=51
        // включена в последовательность именно чтобы показать, что
        // недостаток разрешения (1.5 узла на ℓ) занижает температуру
        // края и завышает улыбку на ~16%.
        let ell = (K_GEL * cfg.t_gel_m / (2.0 * h)).sqrt();
        let nodes_per_ell = ell / (cfg.width_m / (n - 1) as f64);
        if n >= 801 {
            assert!(
                nodes_per_ell >= 20.0,
                "n={n}: только {nodes_per_ell:.1} узлов на ℓ (нужно ≥20)"
            );
        } else {
            assert!(
                nodes_per_ell < 20.0,
                "n={n}: ожидалась грубая сетка (<20 узлов на ℓ), а их {nodes_per_ell:.1}"
            );
        }
        prev = drop;
    }
}

#[test]
fn smile_threshold_submerged_metric_a() {
    // Ревью §3.3 (ревизия 2) / §4а, метрика A — центр против крайней
    // дорожки в 1 см от торца: погружение в буфер даёт U* = 84.0 В/см
    // (Робин) и 62.3 В/см (Дирихле). История: числа «6–7» первой
    // редакции ревью считались по метрике B (центр против самого
    // торца), а 30.6 движка были дефектом краевой строки Робина —
    // смешение объёмной формы k/dy² и потоковой k/dy.
    let base = ThermalConfig {
        medium: WATER,
        kappa_s_m: crate::buffers::KAPPA_TAE_LIMITING_S_M,
        // Историческая таблица §4а остаётся тестом старого решателя.
        cooling: CoolingModel::FreeConvection { c_coef: 0.59 },
        ..ThermalConfig::default()
    };
    let pinned = ThermalConfig {
        edge: EdgeKind::Pinned,
        ..base
    };
    let e_pin = pinned
        .smile_threshold_e(REVIEW_PATH_5CM_CM)
        .expect("тепловой решатель Дирихле")
        .expect("порог Дирихле должен существовать");
    assert!(
        (e_pin - 62.3).abs() < 0.5,
        "U*(Дирихле) = {e_pin:.2}, эталон 62.3 В/см"
    );
    let conv = ThermalConfig {
        edge: EdgeKind::Convective,
        ..base
    };
    let e_conv = conv
        .smile_threshold_e(REVIEW_PATH_5CM_CM)
        .expect("тепловой решатель Робина")
        .expect("порог Робин должен существовать");
    assert!(
        (e_conv - 84.0).abs() < 0.5,
        "U*(Робин) = {e_conv:.2}, эталон 84.0 В/см"
    );
    assert!(
        e_conv > e_pin,
        "Робин обязан держать улыбку до большего поля, чем Дирихле"
    );
}

#[test]
fn smile_thresholds_metric_b_review_table() {
    // Страж таблицы §14 по метрике B (центр против самого торца),
    // ревью §3.3 ревизия 2: погружение 9.37 (Робин, сошедшееся) /
    // 7.36 (Дирихле); открытая кассета 4.21 / 2.06. Именно расхождение
    // метрик A/B стоило целого раунда ревью — без стража эта таблица
    // разъедется первой. Заодно кассета-Дирихле по метрике A: 2.59.
    //
    // Сетка n=6401 ОБЯЗАТЕЛЬНА: пороги метрики B зависят от n в третьей
    // значащей цифре — робиновский край сходится односторонне и первого
    // порядка (9.2963 @401 → 9.3339 @801 → 9.3529 @1601 → 9.3673 @6401
    // → 9.3697 @12801). На n=801 значение 9.33 было занижено; допуск
    // ±0.15 это прятал. Достижимое согласие с таблицей ~0.005,
    // допуск ±0.02 — запас 4×.
    let mk = |medium: Medium, edge: EdgeKind| ThermalConfig {
        medium,
        kappa_s_m: crate::buffers::KAPPA_TAE_LIMITING_S_M,
        // Историческая таблица ревью, не продуктовая модель Э7.
        cooling: CoolingModel::FreeConvection { c_coef: 0.59 },
        edge,
        n_points: 6401,
        ..ThermalConfig::default()
    };
    let cases = [
        (WATER, EdgeKind::Convective, 9.37),
        (WATER, EdgeKind::Pinned, 7.36),
        (AIR, EdgeKind::Convective, 4.21),
        (AIR, EdgeKind::Pinned, 2.06),
    ];
    for (med, edge, expected) in cases {
        let cfg = mk(med, edge);
        let e_b = cfg
            .smile_threshold_e_metric_b(REVIEW_PATH_5CM_CM)
            .expect("тепловой решатель B")
            .unwrap_or_else(|| panic!("{med:?} {edge:?}: порог B не найден"));
        assert!(
            (e_b - expected).abs() < 0.02,
            "{med:?} {edge:?}: U*(B) = {e_b:.4}, эталон {expected}"
        );
    }
    // Кассета-Дирихле, метрика A.
    let e_a = mk(AIR, EdgeKind::Pinned)
        .smile_threshold_e(REVIEW_PATH_5CM_CM)
        .expect("тепловой решатель A")
        .expect("порог кассеты A не найден");
    assert!(
        (e_a - 2.59).abs() < 0.15,
        "кассета Дирихле: U*(A) = {e_a:.2}, эталон 2.59"
    );
    // Метрика B всегда даёт меньший порог: торец холоднее дорожки
    // в 1 см от него, перепад больше.
    for (med, edge, _) in cases {
        let cfg = mk(med, edge);
        let ea = cfg
            .smile_threshold_e(REVIEW_PATH_5CM_CM)
            .expect("тепловой решатель A")
            .expect("порог A");
        let eb = cfg
            .smile_threshold_e_metric_b(REVIEW_PATH_5CM_CM)
            .expect("тепловой решатель B")
            .expect("порог B");
        assert!(eb < ea, "{med:?} {edge:?}: B {eb:.2} ≥ A {ea:.2}");
    }
}

// ---- Нестационарное тепло (закрытие M2) ----

/// Конфигурация для нестационарных тестов: умеренная сетка, TAE.
fn transient_cfg(medium: Medium) -> ThermalConfig {
    let cooling = if medium.k == AIR.k {
        CoolingModel::FreeConvection { c_coef: 0.59 }
    } else {
        ThermalConfig::default().cooling
    };
    ThermalConfig {
        medium,
        cooling,
        kappa_s_m: crate::buffers::KAPPA_TAE_LIMITING_S_M,
        edge: EdgeKind::Convective,
        width_m: 0.10,
        t_gel_m: 0.006,
        n_points: 161,
        l_c_m: 0.08,
        ..ThermalConfig::default()
    }
}

/// Классификатор режима нагрева при полной связке κ(T), h(ΔT):
/// true — разгон (температура ушла за 500 K), false — плато.
/// Горизонт — параметр: около порога рост медленный.
fn is_runaway(cfg: &ThermalConfig, e_v_per_cm: f64, horizon_s: f64) -> bool {
    const KAPPA25: f64 = crate::buffers::KAPPA_TAE_LIMITING_S_M;
    let dt = 10.0;
    let n = cfg.n_points.max(3);
    let mut t = vec![0.0; n];
    let steps = (horizon_s / dt) as usize;
    for _ in 0..steps {
        let d_t_mean = t.iter().sum::<f64>() / n as f64;
        let kappa = crate::circuit::kappa_at(KAPPA25, 298.0 + d_t_mean);
        let q_v = kappa * (e_v_per_cm * 100.0).powi(2);
        let h = cfg.mean_surface_h(d_t_mean) * cfg.h_scale;
        t = cfg.transient_step(&t, dt, h, q_v);
        if t.iter().cloned().fold(f64::MIN, f64::max) > 500.0 {
            return true;
        }
    }
    false
}

/// Бисекция порога разгона по классификатору; порог находит движок.
fn find_runaway_threshold(cfg: &ThermalConfig, lo: f64, hi: f64, horizon_s: f64) -> f64 {
    let (mut lo, mut hi) = (lo, hi);
    for _ in 0..12 {
        let mid = 0.5 * (lo + hi);
        if is_runaway(cfg, mid, horizon_s) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn transient_converges_to_accepted_steady() {
    // Сильнейшая регрессия тепловой части: нестационар при t → ∞ обязан
    // совпасть с принятым steady_profile до 1e-9. Неявный Эйлер имеет ТОТ
    // ЖЕ дискретный неподвижный узел, что и стационарная прогонка
    // (члены ρc_p/dt сокращаются), так что расхождение — чисто
    // экспоненциальный хвост медленной моды: e^{-30} ≈ 1e-13.
    for edge in [EdgeKind::Pinned, EdgeKind::Convective] {
        let mut cfg = review_4a_config(edge);
        cfg.n_points = 321;
        let tau = cfg.heating_time_constant_s();
        let (final_t, _, _) = cfg.transient_run(30.0 * tau, tau / 5.0);
        let (steady_t, _) = cfg.steady_profile().expect("стационарный профиль");
        let worst = final_t
            .iter()
            .zip(steady_t.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        assert!(
            worst < 1e-9,
            "{edge:?}: max |ΔT_нестац − ΔT_стационар| = {worst:.2e} K"
        );
    }
}

#[test]
fn heating_time_constants_reflect_e7_tray_resistance() {
    // Э7: поддон снижает эквивалентный h погружения с прежних ~203 до
    // ~133 Вт/(м²·К), поэтому τ растёт с исторических 62 до ~95 с.
    // Воздушная кассета этой фазой не меняется.
    let water = ThermalConfig {
        medium: WATER,
        kappa_s_m: crate::buffers::KAPPA_TAE_LIMITING_S_M,
        ..ThermalConfig::default()
    };
    let air = ThermalConfig {
        medium: AIR,
        kappa_s_m: crate::buffers::KAPPA_TAE_LIMITING_S_M,
        cooling: CoolingModel::FreeConvection { c_coef: 0.59 },
        ..ThermalConfig::default()
    };
    let tau_w = water.heating_time_constant_s();
    let tau_a = air.heating_time_constant_s();
    println!("τ(вода) = {tau_w:.1} с; τ(воздух) = {tau_a:.0} с");
    assert!(
        (tau_w - 94.7).abs() < 1.0,
        "τ(вода) = {tau_w:.2} с, Э7 ожидает ≈94.7 с"
    );
    assert!(
        (tau_a - 2074.0).abs() < 20.0,
        "τ(воздух) = {tau_a:.1} с, эталон ≈2074 с"
    );
    assert!(
        tau_a > 0.5 * 45.0 * 60.0,
        "τ(воздух) обязана быть сравнима с 45-минутным прогоном"
    );

    let legacy_water = ThermalConfig {
        cooling: CoolingModel::FreeConvection { c_coef: 0.59 },
        ..water
    };
    let legacy_tau = legacy_water.heating_time_constant_s();
    assert!(
        (legacy_tau - 62.0).abs() < 2.0,
        "исторический решатель перестал воспроизводить τ≈62 с: {legacy_tau:.2}"
    );
}

#[test]
fn e7_immersed_layer_is_conductive_at_three_mm_and_tray_is_explicit() {
    let cfg = ThermalConfig::default();
    let (h_top, h_bottom, ra, nu) = cfg.surface_heat_transfer(0.615);
    assert!((ra - 234.78).abs() < 0.1, "Ra={ra:.3}");
    assert_eq!(nu, 1.0, "Ra<1708: слой не должен конвектировать");
    assert!((h_top - WATER.k / 0.003).abs() < 1e-12);
    assert!((h_bottom - 0.20 / 0.003).abs() < 1e-12);
    assert!((cfg.mean_surface_h(0.615) - 133.0).abs() < 1e-12);
}

#[test]
fn e7_buffer_depth_changes_cooling_and_crosses_rayleigh_onset() {
    let cfg_for = |depth_m| ThermalConfig {
        cooling: CoolingModel::ImmersedHorizontal {
            buffer_depth_m: depth_m,
            tray_thickness_m: 0.003,
            tray_k_w_m_k: 0.20,
        },
        ..ThermalConfig::default()
    };
    let (_, _, ra2, nu2) = cfg_for(0.002).surface_heat_transfer(0.615);
    let (h3, _, ra3, nu3) = cfg_for(0.003).surface_heat_transfer(0.615);
    let (h5, _, ra5, nu5) = cfg_for(0.005).surface_heat_transfer(0.615);
    let (h10, _, ra10, nu10) = cfg_for(0.010).surface_heat_transfer(0.615);
    assert!(ra2 < ra3 && ra3 < ra5 && ra5 < 1708.0 && ra10 > 1708.0);
    assert_eq!((nu2, nu3, nu5), (1.0, 1.0, 1.0));
    assert!(nu10 > 1.0, "10 мм: Ra={ra10:.0}, Nu={nu10:.3}");
    assert!(h3 > h10 && h10 > h5, "h3={h3:.1}, h5={h5:.1}, h10={h10:.1}");
}

#[test]
fn air_never_reaches_steady_during_45min() {
    // В открытой кассете гель греется всю дорогу: средний за прогон
    // перегрев меньше половины стационарного, финальный — между
    // половиной и стационаром. Стационарная модель там систематически
    // завышает нагрев (и улыбку).
    let air = transient_cfg(AIR);
    let (final45, mean45, _) = air.transient_run(45.0 * 60.0, 5.0);
    let (steady, _) = air.steady_profile().expect("стационар воздуха");
    let mid = steady.len() / 2;
    assert!(
        mean45[mid] < 0.5 * steady[mid],
        "средний ΔT {:.4} ≥ 0.5·стационарного {:.4}",
        mean45[mid],
        steady[mid]
    );
    assert!(
        final45[mid] > 0.5 * steady[mid] && final45[mid] < steady[mid],
        "финальный ΔT {:.4} вне ({:.4}; {:.4})",
        final45[mid],
        0.5 * steady[mid],
        steady[mid]
    );
}

#[test]
fn runaway_threshold_engine_found_and_ordered_after_smile() {
    // Порог разгона движок находит САМ (бисекция по классификатору
    // уход/плато при полной связке κ(T), h(ΔT)); числа ревью
    // (≈10.04 В/см кассета, ≈89.85 погружение) не зашиты — сверяются
    // независимо. Инвариант вехи: улыбка наступает раньше разгона.
    let air = transient_cfg(AIR);
    let crit_air = find_runaway_threshold(&air, 6.0, 30.0, 6.0 * 3600.0);
    println!("E_разгона(кассета, U-режим) = {crit_air:.2} В/см");
    assert!(
        (5.35..20.0).contains(&crit_air),
        "E_разгона(кассета) = {crit_air:.2} вне окна порядка (5.35; 20)"
    );
    let smile_air = air
        .smile_threshold_e(REVIEW_PATH_5CM_CM)
        .expect("тепловой решатель воздуха")
        .expect("порог улыбки кассеты");
    assert!(
        smile_air < crit_air,
        "инвариант нарушен: улыбка {smile_air:.2} ≥ разгона {crit_air:.2}"
    );

    let water = transient_cfg(WATER);
    let crit_water = find_runaway_threshold(&water, 60.0, 150.0, 2.0 * 3600.0);
    println!("E_разгона(погружение, U-режим) = {crit_water:.2} В/см");
    let smile_water = water
        .smile_threshold_e(REVIEW_PATH_5CM_CM)
        .expect("тепловой решатель воды")
        .expect("порог улыбки погружения");
    assert!(
        smile_water < crit_water && crit_water < 150.0,
        "погружение: улыбка {smile_water:.2}, разгон {crit_water:.2}"
    );
}

#[test]
fn current_mode_plateaus_above_voltage_runaway_field() {
    // Якорь №6b как процесс: в I-режиме E = I/(κS) падает с T, значит
    // q_v = I²/(S²κ(T)) УБЫВАЕТ — плато обязано быть даже при поле,
    // заведомо выше U-порога разгона кассеты.
    let cfg = transient_cfg(AIR);
    let kappa25 = crate::buffers::KAPPA_TAE_LIMITING_S_M;
    let e0 = 15.0_f64; // выше U-порога (~10 В/см)
    let q0 = kappa25 * (e0 * 100.0).powi(2);
    let dt = 10.0;
    let n = cfg.n_points.max(3);
    let mut t = vec![0.0; n];
    for _ in 0..1800 {
        // 5 часов
        let d_t_mean = t.iter().sum::<f64>() / n as f64;
        let kappa = crate::circuit::kappa_at(kappa25, 298.0 + d_t_mean);
        let q_v = q0 * kappa25 / kappa;
        let h = cfg.mean_surface_h(d_t_mean) * cfg.h_scale;
        t = cfg.transient_step(&t, dt, h, q_v);
    }
    let max_t = t.iter().cloned().fold(f64::MIN, f64::max);
    println!("I-режим E₀={e0}: плато max ΔT = {max_t:.1} K");
    assert!(
        max_t < 400.0,
        "I-режим ушёл в разгон: max ΔT = {max_t:.1} K"
    );
}

/// Плато max ΔT (K) при полной связке κ(T), h(ΔT) после горизонта;
/// уход классифицируется как 500+.
fn plateau_max_d_t(cfg: &ThermalConfig, e_v_per_cm: f64, horizon_s: f64) -> f64 {
    // АУДИТ УРОВНЯ 2. κ₂₅ берётся ИЗ КОНФИГА. Прежде помощник жёстко
    // держал `KAPPA_TAE_LIMITING_S_M` и МОЛЧА игнорировал
    // `cfg.kappa_s_m`, который вызывающий заботливо выставлял. Из-за
    // этого продуктовая красная линия считалась ПРЕДЕЛЬНОЙ κ = 0.166,
    // хотя и комментарий рядом с вызовом, и §18 заявляли рабочую
    // 0.1369 — тот самый класс «одна величина, два расчёта», ради
    // которого делалась вся ревизия V12 (§35). Сканер литералов его
    // не поймал: в коде стояла именованная константа, а не число.
    //
    // Эталонный страж `gel_integrity_red_line_bracketed` подаёт
    // предельную κ САМ, через `transient_cfg`, поэтому его вход не
    // меняется и REVIEW-M2 остаётся воспроизведённым (та же
    // дисциплина, что с `FreeConvection { 0.59 }` в §43.4).
    let kappa_25 = cfg.kappa_s_m;
    let dt = 10.0;
    let n = cfg.n_points.max(3);
    let mut t = vec![0.0; n];
    let steps = (horizon_s / dt) as usize;
    for _ in 0..steps {
        let d_t_mean = t.iter().sum::<f64>() / n as f64;
        let kappa = crate::circuit::kappa_at(kappa_25, 298.0 + d_t_mean);
        let q_v = kappa * (e_v_per_cm * 100.0).powi(2);
        let h = cfg.mean_surface_h(d_t_mean) * cfg.h_scale;
        t = cfg.transient_step(&t, dt, h, q_v);
        if t.iter().cloned().fold(f64::MIN, f64::max) > 500.0 {
            return 500.0;
        }
    }
    t.iter().cloned().fold(f64::MIN, f64::max)
}

#[test]
fn gel_integrity_red_line_bracketed() {
    // Красная линия интерфейса — порог ЦЕЛОСТНОСТИ геля, не онсет
    // разгона: кассета при ~6 В/см уже подходит к деградации (66 °C),
    // а онсет (~11.3 В/см) наружу не выходит. Страж: движок бракетует
    // поле достижения ΔT = 65 K (плавление агарозы ~88 °C минус среда)
    // в окно [5; 8.5] В/см. Само поле печатается — сверяется независимо.
    let cfg = transient_cfg(AIR);
    let limit = 65.0_f64;
    let cross = |e: f64| plateau_max_d_t(&cfg, e, 6.0 * 3600.0) > limit;
    let (mut lo, mut hi) = (4.0_f64, 10.0_f64);
    assert!(!cross(lo), "при {lo} В/см не должно плавиться");
    assert!(cross(hi), "при {hi} В/см обязано плавиться");
    for _ in 0..12 {
        let mid = 0.5 * (lo + hi);
        if cross(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let e_red = 0.5 * (lo + hi);
    println!("красная линия целостности геля = {e_red:.2} В/см");
    assert!(
        (5.0..=8.5).contains(&e_red),
        "красная линия {e_red:.2} вне окна [5; 8.5] В/см"
    );
}

#[test]
fn smile_shift_factorises_by_input() {
    // Порог улыбки уехал от эталонных 5.35 к рабочим 7.45 по ДВУМ
    // причинам сразу (κ и путь полосы, §35). Страж проверяет, что
    // вклады МУЛЬТИПЛИКАТИВНО РАЗДЕЛИМЫ: иначе «разложение» в §14 было
    // бы пересказом, а не измерением. Заодно фиксирует базовую точку
    // на входах эталона — 5.35 ревью.
    let cfg = |kappa: f64| ThermalConfig {
        medium: AIR,
        kappa_s_m: kappa,
        cooling: CoolingModel::FreeConvection { c_coef: 0.59 },
        edge: EdgeKind::Convective,
        ..ThermalConfig::default()
    };
    let lim = crate::buffers::KAPPA_TAE_LIMITING_S_M;
    let work = crate::buffers::KAPPA_TAE_S_M;
    let path_work = reference_band_path_cm();
    let e = |k: f64, p: f64| {
        cfg(k)
            .smile_threshold_e(p)
            .expect("тепловой решатель")
            .expect("порог")
    };
    let base = e(lim, REVIEW_PATH_5CM_CM);
    assert!(
        (base - 5.35).abs() < 0.05,
        "базовая точка {base:.4} разошлась с эталоном ревью 5.35"
    );
    let only_kappa = e(work, REVIEW_PATH_5CM_CM) / base;
    let only_path = e(lim, path_work) / base;
    let both = e(work, path_work) / base;
    assert!(
        (only_kappa * only_path / both - 1.0).abs() < 1e-3,
        "вклады не разделяются: {only_kappa:.4}·{only_path:.4} = {:.4} против {both:.4}",
        only_kappa * only_path
    );
    // Направления обязаны быть вверх у обоих: меньше κ — меньше тепла;
    // короче путь — меньше накопленного расхождения фронтов.
    assert!(only_kappa > 1.0 && only_path > 1.0);
}

#[test]
fn open_cassette_uses_horizontal_correlations_e7() {
    // Доработка Э7. §12 приводит ТРИ корреляции McAdams, но в код была
    // заведена только вертикальная 0.59, и ею считались оба лица
    // горизонтального геля. Горячая плита отдаёт вверх вдвое лучше,
    // чем вниз (восходящий факел против запертого под плитой слоя),
    // поэтому одно число на два лица физически неверно.
    let cfg = ThermalConfig {
        medium: AIR,
        cooling: CoolingModel::OpenHorizontal {
            c_up: NU_HORIZONTAL_UP,
            c_down: NU_HORIZONTAL_DOWN,
            tray_thickness_m: DEFAULT_TRAY_THICKNESS_M,
            tray_k_w_m_k: DEFAULT_TRAY_K_W_M_K,
        },
        ..ThermalConfig::default()
    };
    let (up, down, ra, nu) = cfg.surface_heat_transfer(15.0);
    assert!(ra.is_nan() && nu.is_nan(), "у открытой плиты слоя нет");
    assert!(
        up > down,
        "вверх обязано отводить лучше: {up:.2} против {down:.2}"
    );
    // До поддона отношение задано коэффициентами 0.54/0.27. Сам
    // поддон обязан дополнительно уменьшить нижний h.
    let bare_down = h_natural_convection(&AIR, 15.0, cfg.l_c_m, NU_HORIZONTAL_DOWN);
    let ratio = up / bare_down;
    assert!(
        (ratio - NU_HORIZONTAL_UP / NU_HORIZONTAL_DOWN).abs() < 1e-9,
        "отношение лиц {ratio:.4} не равно 0.54/0.27"
    );
    assert!(down < bare_down, "поддон не уменьшил нижний h");
    // Горизонталь охлаждает ХУЖЕ вертикали — иначе доработка ничего
    // не меняла бы и старое число 0.59 было бы консервативным.
    let vertical = ThermalConfig {
        cooling: CoolingModel::FreeConvection { c_coef: 0.59 },
        ..cfg
    };
    let bare_horizontal = ThermalConfig {
        cooling: CoolingModel::OpenHorizontal {
            c_up: NU_HORIZONTAL_UP,
            c_down: NU_HORIZONTAL_DOWN,
            tray_thickness_m: 0.0,
            tray_k_w_m_k: DEFAULT_TRAY_K_W_M_K,
        },
        ..cfg
    };
    let (t_h, h_h) = cfg.steady_profile().expect("профиль с поддоном");
    let (t_bare, h_bare) = bare_horizontal
        .steady_profile()
        .expect("профиль голой плиты");
    let (t_v, h_v) = vertical.steady_profile().expect("вертикальный профиль");
    assert!(
        h_h < h_v,
        "горизонталь охлаждает лучше вертикали: {h_h:.2} против {h_v:.2}"
    );
    assert!(
        h_h < h_bare,
        "поддон не уменьшил средний h: {h_h:.2} против {h_bare:.2}"
    );
    let (dt_h, dt_v) = (
        t_h.iter().cloned().fold(f64::MIN, f64::max),
        t_v.iter().cloned().fold(f64::MIN, f64::max),
    );
    let dt_bare = t_bare.iter().cloned().fold(f64::MIN, f64::max);
    assert!(
        dt_h > dt_v,
        "гель не стал горячее: {dt_h:.2} против {dt_v:.2}"
    );
    assert!(
        dt_h > dt_bare,
        "поддон не сделал гель горячее: {dt_h:.2} против {dt_bare:.2}"
    );
    println!(
        "B cassette: bare h {h_bare:.9} -> tray {h_h:.9} W/(m2 K), max dT {dt_bare:.9} -> {dt_h:.9} K; vertical reference h {h_v:.9}, dT {dt_v:.9}"
    );
}

/// Коридор аудита уровня 2: у положительной κ(T)-связки есть предел
/// существования стационара, и решатель проходит его МОЛЧА.
#[test]
fn kappa_feedback_runaway_is_recorded_not_silent() {
    let path = reference_band_path_cm();
    let air = ThermalConfig {
        medium: AIR,
        kappa_s_m: crate::buffers::KAPPA_TAE_S_M,
        kappa_feedback: true,
        cooling: CoolingModel::OpenHorizontal {
            c_up: NU_HORIZONTAL_UP,
            c_down: NU_HORIZONTAL_DOWN,
            tray_thickness_m: DEFAULT_TRAY_THICKNESS_M,
            tray_k_w_m_k: DEFAULT_TRAY_K_W_M_K,
        },
        ..ThermalConfig::default()
    };
    // Параметр петли: k·A, где A — средний перегрев БЕЗ связки.
    // Стационар существует только при k·A < 1.
    let loop_gain = |e: f64| {
        let (t, _) = ThermalConfig {
            e_v_per_cm: e,
            kappa_feedback: false,
            ..air
        }
        .steady_profile()
        .expect("профиль без κ(T)-связки");
        crate::circuit::K_KAPPA_PER_C * t.iter().sum::<f64>() / t.len() as f64
    };
    let threshold = air
        .smile_threshold_e(path)
        .expect("тепловой решатель")
        .expect("порог");
    let gain_at_threshold = loop_gain(threshold);
    assert!(
        gain_at_threshold < 1.0,
        "опубликованный порог {threshold:.4} В/см стоит в РАЗГОНЕ              (k·A = {gain_at_threshold:.3}) — число нефизично"
    );
    // Граница существования стационара — и она внутри интервала
    // бисекции [0.3; 100], которым порог и ищется.
    let e_runaway =
        bisect(|e| loop_gain(e) - 1.0, 0.3, 100.0).expect("граница разгона обязана бракетиться");
    assert!(
        threshold < e_runaway,
        "порог {threshold:.3} перелез границу разгона {e_runaway:.3} В/см"
    );
    let valid_below = ThermalConfig {
        e_v_per_cm: e_runaway * 0.99,
        ..air
    }
    .spread_mm(path);
    let invalid_above = ThermalConfig {
        e_v_per_cm: e_runaway * 1.01,
        ..air
    }
    .spread_mm(path);
    assert!(valid_below.is_ok(), "профиль ниже границы обязан сходиться");
    assert!(
        matches!(
            invalid_above,
            Err(SteadyProfileError::KappaFeedbackRunaway { .. })
        ),
        "профиль выше границы обязан сообщить κ(T)-разгон: {invalid_above:?}"
    );
    println!(
        "коридор разгона: порог {threshold:.4} В/см при k·A = {gain_at_threshold:.3};              граница стационара {e_runaway:.4} В/см; выше границы возвращается ошибка"
    );
}

#[test]
fn steady_profile_reports_kappa_nonconvergence_level2() {
    let runaway = ThermalConfig {
        medium: AIR,
        e_v_per_cm: 20.0,
        kappa_s_m: crate::buffers::KAPPA_TAE_S_M,
        kappa_feedback: true,
        cooling: CoolingModel::OpenHorizontal {
            c_up: NU_HORIZONTAL_UP,
            c_down: NU_HORIZONTAL_DOWN,
            tray_thickness_m: DEFAULT_TRAY_THICKNESS_M,
            tray_k_w_m_k: DEFAULT_TRAY_K_W_M_K,
        },
        ..ThermalConfig::default()
    };
    assert!(
        runaway.steady_profile().is_err(),
        "за границей κ(T)-разгона решатель обязан вернуть ошибку, а не профиль"
    );
}

#[test]
fn open_cassette_bottom_includes_tray_series_resistance_level2_b() {
    let cfg = ThermalConfig {
        medium: AIR,
        cooling: CoolingModel::OpenHorizontal {
            c_up: NU_HORIZONTAL_UP,
            c_down: NU_HORIZONTAL_DOWN,
            tray_thickness_m: DEFAULT_TRAY_THICKNESS_M,
            tray_k_w_m_k: DEFAULT_TRAY_K_W_M_K,
        },
        ..ThermalConfig::default()
    };
    let delta_t = 15.0;
    let (_, actual_bottom, _, _) = cfg.surface_heat_transfer(delta_t);
    let bare_air = h_natural_convection(&AIR, delta_t, cfg.l_c_m, NU_HORIZONTAL_DOWN);
    let tray = DEFAULT_TRAY_K_W_M_K / DEFAULT_TRAY_THICKNESS_M;
    let expected_bottom = 1.0 / (1.0 / tray + 1.0 / bare_air);
    println!(
        "B guard: h_down bare {bare_air:.9}, tray series expected {expected_bottom:.9}, actual {actual_bottom:.9}"
    );
    assert!(actual_bottom < bare_air, "поддон не уменьшил h_низ");
    assert!((actual_bottom - expected_bottom).abs() < 1e-12);
}

#[test]
fn tray_series_resistance_has_exact_bare_plate_limits_level2_b() {
    let bare_h = 2.571725554_f64;
    let config = |tray_thickness_m, tray_k_w_m_k| ThermalConfig {
        medium: AIR,
        cooling: CoolingModel::OpenHorizontal {
            c_up: NU_HORIZONTAL_UP,
            c_down: NU_HORIZONTAL_DOWN,
            tray_thickness_m,
            tray_k_w_m_k,
        },
        ..ThermalConfig::default()
    };
    let bare_down = |cfg: &ThermalConfig| cfg.surface_heat_transfer(15.0).1;
    assert_eq!(
        bare_down(&config(0.0, DEFAULT_TRAY_K_W_M_K)).to_bits(),
        h_natural_convection(&AIR, 15.0, 0.08, NU_HORIZONTAL_DOWN).to_bits(),
        "нулевая толщина обязана быть бит-в-бит старой плитой"
    );
    assert_eq!(
        bare_down(&config(DEFAULT_TRAY_THICKNESS_M, f64::INFINITY)).to_bits(),
        h_natural_convection(&AIR, 15.0, 0.08, NU_HORIZONTAL_DOWN).to_bits(),
        "идеально проводящий поддон обязан быть бит-в-бит старой плитой"
    );
    assert!(bare_down(&config(DEFAULT_TRAY_THICKNESS_M, DEFAULT_TRAY_K_W_M_K)) < bare_h);
}

#[test]
fn smile_and_red_line_at_working_physics() {
    // ПРОДУКТОВЫЕ числа §14/§18 при ДЕЙСТВУЮЩЕЙ физике движка, в
    // отличие от соседних стражей, воспроизводящих эталон REVIEW-M2 на
    // его собственных входах. Действующая физика — это два числа,
    // которые после M2 успели уехать, а тепловая задача этого не
    // заметила:
    //   * κ рабочая 0.1369 вместо предельной 0.166 (Э1 ревизии V12);
    //   * путь опорной полосы 1.42 см из модели V11 вместо
    //     фольклорных 2.98 см отменённого якоря №1 (PHYSICS §7).
    // Оба сдвига поднимают порог: спред ∝ κ·путь·E³, значит
    // E* ∝ (κ·путь)^(−1/3).
    let path = reference_band_path_cm();
    assert!(
        (1.2..=1.7).contains(&path),
        "опорный путь {path:.3} см разошёлся с якорем №1 V11 (1.47 ± 15 %)"
    );
    let base = ThermalConfig {
        medium: WATER,
        kappa_s_m: crate::buffers::KAPPA_TAE_S_M,
        // Продуктовая таблица обязана считать тот же feedback, что
        // связанный шаг: κ получает среднюю температуру геля.
        kappa_feedback: true,
        ..ThermalConfig::default()
    };
    // Геометрия — РОВНО эталонная (l_c по умолчанию 0.08 м). После
    // Э7 и уровня 2 продуктовая кассета включает горизонтальные
    // корреляции, поддон снизу и κ(T) по средней температуре геля.
    let air = ThermalConfig {
        medium: AIR,
        cooling: CoolingModel::OpenHorizontal {
            c_up: NU_HORIZONTAL_UP,
            c_down: NU_HORIZONTAL_DOWN,
            tray_thickness_m: DEFAULT_TRAY_THICKNESS_M,
            tray_k_w_m_k: DEFAULT_TRAY_K_W_M_K,
        },
        ..base
    };
    let cassette_robin = air
        .smile_threshold_e(path)
        .expect("тепловой решатель кассеты")
        .expect("порог кассеты");
    let cassette_pinned = ThermalConfig {
        edge: EdgeKind::Pinned,
        ..air
    }
    .smile_threshold_e(path)
    .expect("тепловой решатель кассеты, Дирихле")
    .expect("порог кассеты, Дирихле");
    println!(
        "рабочие пороги улыбки (метрика A): кассета Робин {cassette_robin:.2},              Дирихле {cassette_pinned:.2} В/см; опорный путь {path:.3} см"
    );
    // Порядок обязателен: Дирихле охлаждает торец сильнее → улыбка позже.
    assert!(
        cassette_pinned < cassette_robin,
        "Дирихле {cassette_pinned:.2} обязан быть ниже Робина {cassette_robin:.2}"
    );
    // АУДИТ УРОВНЯ 2: прежнее «пороги обязаны быть ВЫШЕ эталонных
    // 5.35/2.59» после включения κ(T)-связки держится на 0.03 В/см и
    // проверять перестало что-либо — два сдвига разного знака сошлись
    // почти в точку эталона. Утверждение снято, вместо него — узкий
    // страж опубликованной таблицы и проверка ОБЛАСТИ ОПРЕДЕЛЕНИЯ.
    assert!(
        // V14 (C1, §53): опорный путь 1.42 → 1.67 см ⇒ пороги ∝
        // путь^(−1/3) опустились с 5.38 / 3.12 до 5.19 / 2.95 В/см.
        (cassette_robin - 5.19).abs() < 0.02 && (cassette_pinned - 2.95).abs() < 0.02,
        "таблица §46 разошлась с расчётом: {cassette_robin:.4} / {cassette_pinned:.4} В/см"
    );
    // Область определения по ТЕМПЕРАТУРЕ. Главный урок §18 —
    // «проверка области определения обязана предшествовать проверке
    // точности» — до сих пор применялся только к геометрии (влезает
    // ли полоса в гель). Температуру не проверял никто, а открытая
    // кассета на своём же пороге улыбки уже РАСПЛАВЛЕНА: движок
    // классифицирует это состояние как MeltRisk.
    let max_overheat = |cfg: &ThermalConfig, e: f64| -> f64 {
        let (t, _) = ThermalConfig {
            e_v_per_cm: e,
            ..*cfg
        }
        .steady_profile()
        .expect("профиль области определения");
        t.iter().cloned().fold(f64::MIN, f64::max)
    };
    let melt_overheat_k = crate::chamber::MELT_RISK_K - crate::mobility::T_REF_K;
    let e_melt = bisect(|e| max_overheat(&air, e) - melt_overheat_k, 0.5, 8.0)
        .expect("порог плавления кассеты обязан бракетиться");
    let overheat_at_smile = max_overheat(&air, cassette_robin);
    println!(
        "область определения кассеты: плавление при {e_melt:.4} В/см,              перегрев на пороге улыбки {overheat_at_smile:.2} K (предел {melt_overheat_k:.1} K)"
    );
    assert!(
        e_melt < cassette_robin,
        "кассета перестала плавиться раньше улыбки ({e_melt:.3} против              {cassette_robin:.3} В/см) — область определения изменилась,              перевыпусти §46 и сними эту запись"
    );
    assert!(
        overheat_at_smile > melt_overheat_k,
        "перегрев на пороге улыбки {overheat_at_smile:.2} K опустился              под предел плавления — то же самое"
    );
    // Область определения (§14): путь при пороговом поле помещается в гель.
    let travel = air.band_path_cm(cassette_robin, path);
    assert!(
        travel <= 10.0,
        "порог {cassette_robin:.2} В/см гонит полосу на {travel:.1} см — вне области определения"
    );
    // Красная линия целостности при рабочей κ.
    // Красная линия — тоже ПРОДУКТОВОЕ число (§18 выносит её в
    // интерфейс), поэтому и она считается горизонтальной плитой.
    let cfg_red = ThermalConfig {
        kappa_s_m: crate::buffers::KAPPA_TAE_S_M,
        cooling: CoolingModel::OpenHorizontal {
            c_up: NU_HORIZONTAL_UP,
            c_down: NU_HORIZONTAL_DOWN,
            tray_thickness_m: DEFAULT_TRAY_THICKNESS_M,
            tray_k_w_m_k: DEFAULT_TRAY_K_W_M_K,
        },
        ..transient_cfg(AIR)
    };
    let cross = |e: f64| plateau_max_d_t(&cfg_red, e, 6.0 * 3600.0) > 65.0;
    let (mut lo, mut hi) = (4.0_f64, 12.0_f64);
    assert!(
        !cross(lo) && cross(hi),
        "красная линия не бракетится в [4; 12]"
    );
    for _ in 0..12 {
        let mid = 0.5 * (lo + hi);
        if cross(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let e_red = 0.5 * (lo + hi);
    println!("рабочая красная линия целостности = {e_red:.2} В/см");
    // κ упала ⇒ при том же поле греет слабее ⇒ линия обязана подняться.
    assert!(
        (5.5..=9.5).contains(&e_red),
        "рабочая красная линия {e_red:.2} вне окна [5.5; 9.5] В/см"
    );
    assert!(
        (e_red - 6.2334).abs() < 0.02,
        "красная линия таблицы §47 разошлась с расчётом: {e_red:.4} В/см"
    );
}

#[test]
fn smile_threshold_requires_band_inside_gel() {
    // Страж области определения (§14): порог улыбки осмыслен, только
    // если путь полосы при пороговом поле помещается в гель.
    // Кассета: путь при 5.35 В/см — 3.2 см < 10 см, порог физический.
    // Погружение: путь при «пороге» 84 В/см — 50 см, в 5 раз длиннее
    // геля ⟹ 84.0/62.3 — граница модели, не рабочие величины.
    let air = transient_cfg(AIR);
    let path_cassette = air.band_path_cm(5.35, REVIEW_PATH_5CM_CM);
    assert!(
        path_cassette < 10.0,
        "кассета: путь {path_cassette:.1} см при пороговом поле не помещается"
    );
    let water = transient_cfg(WATER);
    let path_submerged = water.band_path_cm(84.0, REVIEW_PATH_5CM_CM);
    assert!(
        path_submerged > 10.0,
        "погружение: путь {path_submerged:.1} см обязан выходить за гель"
    );
}
