//! Calibration, limiting behavior and conformation regression tests.

use super::*;
use crate::{fmath::DetMath, mobility::T_REF_K};

#[test]
fn holmes1990_table2_absolute_curve() {
    // ГЛАВНОЕ КАЛИБРОВОЧНОЕ ВОРОТО (V11, 2026-08-27). Holmes &
    // Stellwagen 1990, таблица 2, столбец «LE, TAE, без EtBr»,
    // 1 % агароза, µ₂₀ экстраполированы к E = 0; нормировка — на
    // «intrinsic mobility» 2.7e-4 той же работы (docs/datasets/
    // holmes1990.md, раздел «НОРМИРОВКА»). Полоса ±8 % — это
    // измерительная точность таблицы плюс округление констант;
    // ослаблять её нельзя, сужать при улучшении модели — нужно.
    for (bp, tgt) in [
        (510.0, 1.85 / 2.7),
        (1000.0, 1.52 / 2.7),
        (2000.0, 1.07 / 2.7),
        (4100.0, 0.62 / 2.7),
        (8100.0, 0.29 / 2.7),
        (12200.0, 0.16 / 2.7),
    ] {
        let m = MobilityV2::mobility_ratio(bp, 1.0, 0.0);
        let d = m / tgt;
        assert!(
            (0.92..=1.08).contains(&d),
            "Holmes {bp:.0} бп: модель {m:.4} против {tgt:.4} (×{d:.3})"
        );
    }
}

#[test]
fn holmes1990_fig1c_field_ratios() {
    // Fig. 1c той же работы: отношение µ(4.6 В/см)/µ(0.5 В/см)
    // в 1 % LE/TAE растёт с массой — 1 кб почти не ускоряется,
    // 12.2 кб ускоряется втрое. Полоса ±25 % (числа сняты с графика).
    for (bp, tgt) in [
        (1000.0, 1.03),
        (4100.0, 1.45),
        (8100.0, 2.20),
        (12200.0, 3.30),
    ] {
        let r = MobilityV2::mobility_ratio(bp, 1.0, 4.6) / MobilityV2::mobility_ratio(bp, 1.0, 0.5);
        assert!(
            (0.75 * tgt..=1.25 * tgt).contains(&r),
            "{bp:.0} бп: µ(4.6)/µ(0.5) = {r:.2}, у Holmes {tgt}"
        );
    }
}

#[test]
fn plateau_and_field_ratio_of_lambda_dna() {
    // Плато 48.5 кб при 1 %/5 В/см — ИЗМЕРЕННОЕ окно (V14, §53.2;
    // §40.5 требовал заменить «порядок 3» Viovy измерением): Rill 2002
    // µl(5 В/см, 1 %) = 0.67e-4 в TBE/20 °C → в TAE/25 °C
    // ×(3.68/4.5)·1.10 ≈ 0.60e-4; Akerman 1996 (по Viovy) 1.1e-4 в TBE
    // → ≈ 1.0e-4. Окно [0.55e-4; 1.10e-4]; V14 даёт 0.75e-4.
    // Отношение 8/1 В/см в 2–3 (R1 §2: сублинейность v(E)).
    let mu = MobilityV2::mobility_cm2_v_s(48502.0, 1.0, 5.0, T_REF_K);
    assert!(
        (0.55e-4..=1.10e-4).contains(&mu),
        "плато {mu:.3e} вне измеренного окна [0.55e-4; 1.10e-4]"
    );
    let r = MobilityV2::mobility_ratio(48502.0, 1.0, 8.0)
        / MobilityV2::mobility_ratio(48502.0, 1.0, 1.0);
    assert!((2.0..=3.0).contains(&r), "отношение 8/1 В/см = {r:.2}");
}

#[test]
fn no_band_inversion_at_1pct_matches_rill_v_q13() {
    // ПЕРЕПИСАНО В V14 (C1, §53.2). Ворота V11 требовали СЛАБОЙ
    // инверсии 15↔23 кб при 1 В/см в 1 % геле — числа, которого в
    // данных нет: Rill 2002 при 1.00 В/см на 1.0 % даёт монотонный
    // фит до 48.5 кб (χ² 0.9997, отклонения лишь > 25 кб), а Heller
    // H6 видит инверсию в 2 % гелях при низком поле. Поэтому при 1 %
    // лестница обязана быть МОНОТОННОЙ при 1 и 5 В/см на всём
    // диапазоне 100 бп – 48.5 кб; отношение 23/15 кб записывается.
    for e in [1.0, 5.0] {
        let mut prev = f64::INFINITY;
        let mut bp = 100.0f64;
        while bp <= 48_500.0 {
            let m = MobilityV2::mobility_ratio(bp, 1.0, e);
            // Допуск 0.1 %: выше ~25 кб при 1 В/см кривая — плато с
            // рябью в сотые процента (хвост окна трэппинга), а не
            // инверсия; сама статья видит отклонения там же (> 25 кб).
            assert!(
                m <= prev * (1.0 + 1e-3),
                "инверсия при 1 %, {e} В/см у {bp:.0} бп: {m:.4} > {prev:.4} — \
                 Rill 2002 такого не видит (V-Q13)"
            );
            prev = m;
            bp *= 1.1;
        }
    }
    let inv = |c: f64, e: f64| {
        MobilityV2::mobility_ratio(23_000.0, c, e) / MobilityV2::mobility_ratio(15_000.0, c, e)
    };
    let (r1, r5) = (inv(1.0, 1.0), inv(1.0, 5.0));
    assert!(
        (0.85..1.0).contains(&r1) && (0.85..1.0).contains(&r5),
        "запись V14: µ(23)/µ(15) при 1 % = {r1:.3} (1 В/см), {r5:.3} (5 В/см)"
    );
    // 23 кб НИКОГДА не обгоняет 1 кб — на сетке C × E.
    for c in [0.7, 1.0, 1.45, 2.0] {
        for e in [0.5, 1.0, 2.0, 5.0] {
            let r = MobilityV2::mobility_ratio(23_000.0, c, e)
                / MobilityV2::mobility_ratio(1000.0, c, e);
            assert!(r < 1.0, "23 кб обогнала 1 кб при C={c}, E={e}");
        }
    }
    // Точка инверсии модели в плотном геле при слабом поле (условия
    // Heller H6) — ЗАПИСЬ, не ворота: V-Q13 остаётся открытым, пока
    // центр окна трэппинга не едет с ε (VALIDATION V-Q13).
    let r2 = inv(2.0, 0.5);
    assert!(
        r2.is_finite() && (0.5..2.0).contains(&r2),
        "µ(23)/µ(15) при 2 %, 0.5 В/см = {r2:.3}"
    );
}

#[test]
fn mu_linear_in_e_never_quadratic_and_n_star_scales() {
    // §4.7 + Heller H3/H1: µ монотонна по E, локальный показатель
    // dlnµ/dlnE ≤ 1 ВЕЗДЕ (не только в насыщении), N* ∝ E^−0.65…−0.80.
    for &n in &[500.0, 5000.0, 12000.0, 23000.0, 80000.0] {
        let mut prev = MobilityV2::mobility_ratio(n, 1.0, 0.5);
        let mut e = 0.5_f64;
        while e < 10.0 {
            let e2 = e * 1.15;
            let cur = MobilityV2::mobility_ratio(n, 1.0, e2);
            assert!(cur >= prev - 1e-12, "µ(E) немонотонна при N={n}");
            let slope = (cur / prev).dln() / (e2 / e).dln();
            // Heller 1994, ГЛАВНЫЙ вывод (H3): «could not find any
            // regime with a dependence of mobility on field strength
            // that is quadratic or even larger than linear» — то есть
            // показатель ≤ 1 ВЕЗДЕ, а не только в насыщении. До V11
            // модель доходила до 3.1 (разрывный трэппинг).
            assert!(slope <= 1.0 + 1e-9, "µ~E^{slope:.2} > E¹ при N={n}, E={e}");
            prev = cur;
            e = e2;
        }
    }
    // Heller 1994 (H1): N* ∝ E^−q с ИЗМЕРЕННЫМ q = 0.65…0.80 — это
    // BRF, а не BRM (который дал бы q = 2). До V11 модель держала
    // ровно q = 1, то есть была круче всех измеренных наклонов.
    let r = MobilityV2::n_star_blobs(1.0, 1.0) / MobilityV2::n_star_blobs(4.0, 1.0);
    let q = r.dln() / 4.0_f64.dln();
    assert!(
        (0.63..=0.83).contains(&q),
        "показатель N* ∝ E^−{q:.3} вне измеренного окна Heller 0.65–0.80"
    );
}

#[test]
fn entropic_trapping_is_bounded_windowed_and_field_quenched() {
    // ПЕРЕПИСАНО В V11 (закрывает V-Q8). Требования литературы:
    // (1) Heller H5 — эффект в ДЕСЯТКИ ПРОЦЕНТОВ, µ конечна всегда;
    // (2) Rill 2002 — аномалия гаснет выше ~2.5 В/см, ГЛАДКО;
    // (3) Heller H5 — локализован по размеру (5–20 кб), а не растёт
    //     монотонно с N: 48.5 кб при 2 % его уже не чувствуют;
    // (4) Holmes — наклон круче −1 виден уже при 1 %, значит порога
    //     по концентрации нет.
    let f = |n: f64, c: f64, e: f64| MobilityV2::trapping_factor(MobilityV2::blobs(n, c), c, e);
    // (1) ограниченность снизу — по построению 1 − A_TRAP.
    for n in [1e3, 1e4, 1e5, 1e6] {
        for c in [0.5, 1.0, 2.0, 3.0] {
            for e in [0.0, 0.13, 1.0, 5.0, 20.0] {
                let v = f(n, c, e);
                assert!(
                    (1.0 - A_TRAP..=1.0).contains(&v),
                    "множитель {v:.3} вне [1−A_TRAP; 1] при N={n}, C={c}, E={e}"
                );
            }
        }
    }
    // (2) гладкое гашение полем и монотонность роста µ(E) у порога.
    assert!(f(12000.0, 1.0, 0.5) < f(12000.0, 1.0, 2.5));
    assert!(f(12000.0, 1.0, 2.5) < f(12000.0, 1.0, 10.0));
    let at = MobilityV2::mobility_ratio(12000.0, 1.4, E_TRAP_OFF_V_CM);
    let above = MobilityV2::mobility_ratio(12000.0, 1.4, E_TRAP_OFF_V_CM + 0.01);
    assert!(above >= at, "рост µ(E) сломан у поля гашения трэппинга");
    // (3) окно по размеру: 12 кб при 1 % внутри, 48.5 кб при 2 % — вне.
    assert!(f(12000.0, 1.0, 0.5) < 0.85, "окно не ловит 12 кб при 1 %");
    assert!(
        f(48502.0, 2.0, 0.13) > 0.98,
        "48.5 кб при 2 % не должна быть в окне"
    );
    // (4) при 1 % трэппинг ЕСТЬ (порога по C нет), и µ конечна.
    assert!(MobilityV2::mobility_ratio(48502.0, 2.0, 0.13) > 1e-3);
}

#[test]
fn crossover_matches_holmes_onset_criterion() {
    // ПЕРЕПИСАНО В V11. Прежде центр кроссовера сравнивался с
    // параметром γ фита Ван Винкла (8–12 кб) — но γ у него описывает
    // ВСЮ кривую 1/µ(L), а не точку смены механизма, и калибровка на
    // него уводила модель от абсолютов Holmes.
    // Критерий Holmes (§3.1) прямой и измеренный: полевая зависимость
    // µ начинается при 1.6 кб в 1 % TAE — там, где Rg сравнивается с
    // радиусом поры. Это и есть Огстон→рептация.
    let lx = MobilityV2::crossover_bp(1.0);
    assert!(
        (600.0..=2600.0).contains(&lx),
        "L× = {lx:.0} п.о. далеко от начала полевой зависимости              Holmes (1.6 кб)"
    );
    // Прямая проверка самого критерия: ниже L× поле µ не трогает,
    // выше — трогает.
    let ratio = |bp: f64| {
        MobilityV2::mobility_ratio(bp, 1.0, 4.6) / MobilityV2::mobility_ratio(bp, 1.0, 0.5)
    };
    assert!(ratio(0.5 * lx) < 1.02, "поле действует ниже кроссовера");
    assert!(ratio(4.0 * lx) > 1.15, "поле не действует выше кроссовера");
    // a(C) растёт при разбавлении → L× тоже.
    assert!(MobilityV2::crossover_bp(0.7) > lx);
}

#[test]
fn ladder_order_below_crossover_is_monotone() {
    // Порядок лестницы ≤5 кб при 1%/5 В/см строго убывает —
    // видимая лестница не перепутана.
    let ns = [
        100.0, 200.0, 300.0, 500.0, 800.0, 1200.0, 2000.0, 3000.0, 5000.0,
    ];
    let mut prev = f64::INFINITY;
    for &n in &ns {
        let m = MobilityV2::mobility_ratio(n, 1.0, 5.0);
        assert!(m < prev, "порядок нарушен на {n} бп");
        prev = m;
    }
}

#[test]
fn pluen_diffusion_branches_and_field_dispersion() {
    // Ветки Pluen: внутри Rouse-ветки D(1к)/D(5к) = 5^0.52 ≈ 2.31;
    // сшивка на 6 кб непрерывна; рептационная ветка круче
    // (D(10к)/D(40к) = 4^1.55 ≈ 8.6); полевая добавка растёт с E и
    // гаснет при E → 0.
    let d = |n: f64, e: f64| MobilityV2::dna_diffusivity_cm2_s(n, 1.0, e);
    let r_rouse = d(1000.0, 0.0) / d(5000.0, 0.0);
    assert!((r_rouse - 5.0_f64.dpowf(0.52)).abs() < 1e-9);
    let below = d(5999.9, 0.0);
    let above = d(6000.1, 0.0);
    assert!((below / above - 1.0).abs() < 1e-3, "излом разорван");
    let r_rept = d(10000.0, 0.0) / d(40000.0, 0.0);
    assert!((r_rept - 4.0_f64.dpowf(1.55)).abs() < 1e-6);
    assert!(d(23000.0, 8.0) > d(23000.0, 2.0));
    // «Низкие поля/короткие фрагменты — D почти не зависит от E»
    // (Тинланд). R2: ниже одного блоба включение s(nb) = 1 − exp(−nb²)
    // оставляет 500 бп при 1 % (0.15 блоба) не более ×1.15 даже при
    // 8 В/см — конвенция, данных короче 5.7 кб нет (PHYSICS §54.3).
    let short_gain = d(500.0, 8.0) / d(500.0, 0.0);
    assert!(
        (1.0..1.15).contains(&short_gain),
        "короткий уширяется полем: {short_gain:.4}"
    );
    assert!(MobilityV2::blobs(2000.0, 1.0) < 1.0, "2 кб внутри блоба");
    assert!(MobilityV2::blobs(5000.0, 1.0) > 1.0, "5 кб длиннее блоба");
    let long_gain = d(48502.0, 8.0) / d(48502.0, 0.0);
    assert!(
        long_gain > 3.0,
        "крупный не уширяется полем: {long_gain:.2}"
    );
    assert!(
        (d(23000.0, 0.0) / (D_1KB_CM2_S * 6.0_f64.dpowf(-0.52) * (23.0_f64 / 6.0).dpowf(-1.55))
            - 1.0)
            .abs()
            < 1e-9
    );
}

#[test]
fn conformer_order_i_iii_ii_at_1pct() {
    // Классический порядок 0.7–1% без EtBr (R1 §4, BitesizeBio):
    // суперскрученная (I) > линейная (III) > релаксированная (II).
    let bp = 5400.0; // pSM1-масштаб
    let sc =
        MobilityV2::form_mobility_ratio(bp, DnaForm::Supercoiled { sigma0: -0.06 }, 1.0, 5.0, 0.0);
    let lin = MobilityV2::form_mobility_ratio(bp, DnaForm::Linear, 1.0, 5.0, 0.0);
    let rel = MobilityV2::form_mobility_ratio(bp, DnaForm::Relaxed, 1.0, 5.0, 0.0);
    assert!(sc > lin * 1.03, "I не быстрее III: {sc:.3} vs {lin:.3}");
    assert!(lin > rel * 1.05, "III не быстрее II: {lin:.3} vs {rel:.3}");
}

#[test]
fn e6_zero_etbr_is_bitwise_old_geometry() {
    for bp in [500.0, 5000.0, 23_130.0] {
        let plain = MobilityV2::mobility_ratio(bp, 1.0, 5.0);
        let through_form = MobilityV2::form_mobility_ratio(bp, DnaForm::Linear, 1.0, 5.0, 0.0);
        assert_eq!(plain.to_bits(), through_form.to_bits(), "{bp} п.о.");
    }
}

#[test]
fn e6_intercalation_never_accelerates_dna() {
    // ФИЗИЧЕСКИЙ ЗАПРЕТ, найденный аудитом Э6. Геометрический маршрут
    // ν → удлинение контура двигает молекулу вдоль немонотонной µ(N),
    // и у резкого кроссовера Огстон↔BRF (ширина 0.128 декады)
    // перебрасывает её в БОЛЕЕ БЫСТРУЮ ветвь. До пола сетка
    // (N, C, E, доза) давала 1061 отрицательную точку из 8160 (13 %),
    // худшая −7.99 % при 2388 бп / 0.5 % / 8 В/см: краситель разгонял
    // ДНК почти на восемь процентов.
    //
    // Удлинение контура, рост Rg и добавленная масса — все механизмы
    // тормозят, поэтому знак здесь не обсуждается.
    let mut worst = 0.0f64;
    let mut checked = 0usize;
    for &c in &[0.5f64, 0.7, 1.0, 1.5, 2.0] {
        for &e in &[0.5f64, 0.83, 1.0, 2.0, 5.0, 8.0] {
            for &ug in &[0.25f64, 0.5, 1.0, 2.0] {
                let free = crate::stain::ug_ml_to_molar(ug);
                let mut bp = 100.0f64;
                while bp <= 200_000.0 {
                    for form in [DnaForm::Linear, DnaForm::Relaxed] {
                        let f = intercalation_factor(bp, form, c, e, free);
                        worst = worst.min(1.0 - f);
                        checked += 1;
                    }
                    bp *= 1.2;
                }
            }
        }
    }
    assert!(checked > 3_000, "сетка стража усохла: {checked} точек");
    assert!(
        worst >= -1e-12,
        "интеркаляция УСКОРЯЕТ ДНК: худшее торможение {:.4} %",
        worst * 100.0
    );

    // Пол не должен трогать путь без красителя — там всё побитово.
    for bp in [500.0, 5000.0, 23_130.0] {
        for form in [DnaForm::Linear, DnaForm::Relaxed] {
            let a = MobilityV2::form_mobility_ratio(bp, form, 1.0, 5.0, 0.0);
            let b = MobilityV2::form_mobility_ratio(bp, form, 1.0, 5.0, 0.0);
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    // Суперспираль полу НЕ подчиняется: там ν меняет ТОПОЛОГИЮ, и
    // немонотонность — якорь №13, а не артефакт сшивки. Проверяем
    // сам провал: у точки релаксации (σ_eff → 0) плектонема
    // распускается в кольцо и тормозится, а за ней снова компактится
    // положительными витками и ВОССТАНАВЛИВАЕТСЯ.
    let sc = DnaForm::Supercoiled { sigma0: -0.06 };
    let mu_sc = |ug: f64| {
        MobilityV2::form_mobility_ratio(5000.0, sc, 1.0, 5.0, crate::stain::ug_ml_to_molar(ug))
    };
    let (clean, dip, recovered) = (mu_sc(0.0), mu_sc(0.32), mu_sc(5.0));
    assert!(
        dip < clean,
        "провал дай-титрования исчез: {dip:.4} против {clean:.4}"
    );
    assert!(
        recovered > dip,
        "восстановление за точкой инверсии задавлено полом:              {recovered:.4} против провала {dip:.4}"
    );
}

/// Локальный помощник стража: множитель EtBr для формы.
fn intercalation_factor(bp: f64, form: DnaForm, c: f64, e: f64, free: f64) -> f64 {
    MobilityV2::etbr_geometry_factor(bp, form, c, e, free)
}

#[test]
fn e6_linear_etbr_is_geometric_and_exposes_upper_plateau_v_q19() {
    let factor = |bp: f64, ug_ml: f64| {
        MobilityV2::etbr_geometry_factor(
            bp,
            DnaForm::Linear,
            1.0,
            0.83,
            crate::stain::ug_ml_to_molar(ug_ml),
        )
    };
    let sizes = [2027.0, 4361.0, 9416.0, 23_130.0];
    let slow: Vec<f64> = sizes.iter().map(|&bp| 1.0 - factor(bp, 0.5)).collect();
    println!("Э6 slowdown λ/HindIII @0.5 мкг/мл (геометрия): {slow:?}");
    assert!(slow.iter().all(|&v| v > 0.0 && v < 0.4));
    assert!(slow[0] < slow[1] && slow[1] < slow[2]);
    // Геометрия одна: 23 кб почти не реагирует — верхняя ветвь сидит
    // на N-независимом ориентационном плато (было V-Q19). Записано как
    // структурная черта геометрической части: 9.4 кб ≈ 9 %, 23 кб ≈ 0.4 %.
    assert!(
        slow[3] < 0.02 && slow[2] > 0.08,
        "геометрический коридор уехал: {slow:?}; перепроверь форму верхней ветви"
    );
    // C2 (PHYSICS §55): полный множитель = геометрия × трение по
    // изотерме; трение от N не зависит, поэтому 23 кб тормозится на
    // ~18 %, а 9.4 кб — на ~25 % (Sigmon Table 1 держится, V-Q10 и V-Q19
    // закрыты). Записанный коридор полного торможения при 0.5 мкг/мл.
    let total: Vec<f64> = sizes
        .iter()
        .map(|&bp| {
            1.0 - MobilityV2::etbr_total_factor(
                bp,
                DnaForm::Linear,
                1.0,
                0.83,
                crate::stain::ug_ml_to_molar(0.5),
            )
        })
        .collect();
    println!("C2 slowdown λ/HindIII @0.5 мкг/мл (полный): {total:?}");
    assert!(
        (0.15..=0.22).contains(&total[3]) && (0.22..=0.28).contains(&total[2]),
        "коридор C2 уехал: {total:?}"
    );

    let almost_sat = factor(5000.0, 500.0);
    let saturated = factor(5000.0, 5000.0);
    assert!(
        (almost_sat - saturated).abs() < 0.01,
        "изотерма не дала насыщения: {almost_sat:.4} против {saturated:.4}"
    );
}

#[test]
fn anchor_sc_lin_reversal_by_gel_concentration() {
    // Ворота этапа 3 (MASTERPLAN, R1 §4 по Mickel 1977): существует
    // концентрация геля C*, где sc комигрирует с линейной той же
    // длины; ниже — sc быстрее, выше — медленнее. Для 15 кб C* в
    // рабочем окне [1.2; 2.4]; для 25 кб C* НИЖЕ (реверсия при MW,
    // убывающем с ростом C) — эмерджентно.
    let diff = |bp: f64, c: f64| {
        MobilityV2::form_mobility_ratio(bp, DnaForm::Supercoiled { sigma0: -0.06 }, c, 5.0, 0.0)
            - MobilityV2::form_mobility_ratio(bp, DnaForm::Linear, c, 5.0, 0.0)
    };
    let cstar = |bp: f64| -> f64 {
        let (mut lo, mut hi) = (0.5_f64, 3.0_f64);
        assert!(diff(bp, lo) > 0.0, "sc не быстрее при 0.5%");
        assert!(diff(bp, hi) < 0.0, "sc не медленнее при 3%");
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if diff(bp, mid) > 0.0 {
                lo = mid
            } else {
                hi = mid
            }
        }
        0.5 * (lo + hi)
    };
    let c15 = cstar(15000.0);
    let c25 = cstar(25000.0);
    assert!((1.2..=2.4).contains(&c15), "C*(15 кб) = {c15:.2}");
    assert!(c25 < c15, "C* не убывает с MW: {c25:.2} vs {c15:.2}");
}

#[test]
fn etbr_titration_dip_of_supercoil() {
    // Дай-титрование (R1 §5): µ_sc(c_EtBr) падает до минимума в точке
    // полной релаксации (σ_eff = 0) и растёт при позитивной
    // суперспирализации; точка минимума согласована с окном №13.
    let bp = 5400.0;
    let mu_at = |ug_ml: f64| {
        MobilityV2::form_mobility_ratio(
            bp,
            DnaForm::Supercoiled { sigma0: -0.06 },
            1.0,
            5.0,
            crate::stain::ug_ml_to_molar(ug_ml),
        )
    };
    let c_relax = crate::stain::inversion_free_ug_ml(-0.06);
    let dip = mu_at(c_relax);
    assert!(dip < mu_at(0.0) * 0.75, "провала при релаксации нет");
    assert!(
        dip < mu_at(3.0 * c_relax) * 0.999,
        "восстановления после нет"
    );
    // Минимум лежит у точки релаксации: сетка грубой проверки.
    let mut min_c = 0.0;
    let mut min_mu = f64::INFINITY;
    for i in 0..60 {
        let c = 0.02 * i as f64;
        let m = mu_at(c);
        if m < min_mu {
            min_mu = m;
            min_c = c;
        }
    }
    assert!(
        (min_c - c_relax).abs() < 0.1,
        "минимум при {min_c:.2} мкг/мл, релаксация при {c_relax:.2}"
    );
}

#[test]
fn form_ii_effectively_stalls_at_high_mw() {
    // Mickel 1977: форма II выше порогового MW «вообще не входит в
    // гель» — эмерджентно из огстон-ветки кольца с насаживанием.
    let lin = MobilityV2::form_mobility_ratio(30000.0, DnaForm::Linear, 1.3, 5.0, 0.0);
    let rel = MobilityV2::form_mobility_ratio(30000.0, DnaForm::Relaxed, 1.3, 5.0, 0.0);
    // Порог записан ЧЕРЕЗ ВИДИМОЕ: подавление в 10 раз означает, что
    // за 45 мин при 5 В/см форма II 30 кб проходит ~1 мм — это и есть
    // «не входит в гель» (Mickel 1977) на глаз. До перекалибровки V11
    // множитель насаживания был жёстче (0.05), но там и огстоновская
    // ветвь колец была другой; утверждение Mickel качественное.
    assert!(
        rel < 0.15 * lin,
        "форма II 30 кб не встала: {rel:.2e} vs {lin:.3}"
    );
}

#[test]
fn regimes_are_reported_sanely() {
    assert_eq!(MobilityV2::regime(300.0, 1.0, 5.0), MigrationRegime::Ogston);
    assert_eq!(
        MobilityV2::regime(48502.0, 1.0, 5.0),
        MigrationRegime::OrientedReptation
    );
    // Трэппинг после V11 — ОКНО по размеру (Heller H5: 5–20 кб), а не
    // «всё крупное при плотном геле»: 12 кб при 1 % и 0.5 В/см внутри
    // окна, а 48.5 кб при 2 % — уже вне его.
    assert_eq!(
        MobilityV2::regime(12000.0, 1.0, 0.5),
        MigrationRegime::EntropicTrapping
    );
    assert_ne!(
        MobilityV2::regime(48502.0, 2.0, 0.13),
        MigrationRegime::EntropicTrapping
    );
    // Неориентированная рептация: за кроссовером, но до N*.
    // (N* = K_ORIENT/ε^Q_ORIENT — при слабом поле велик.)
    assert!(3000.0 > MobilityV2::crossover_bp(1.0));
    assert!(MobilityV2::blobs(3000.0, 1.0) < MobilityV2::n_star_blobs(0.05, 1.0));
    assert_eq!(
        MobilityV2::regime(3000.0, 1.0, 0.05),
        MigrationRegime::Reptation
    );
}
