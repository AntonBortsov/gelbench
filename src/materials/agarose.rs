//! Радиус инерции червеобразной цепи из корреляции концов подцепи (Doi & Edwards 1986, §2.3)
//! и контурного интеграла для Rg² (Rubinstein & Colby 2003, §2). Для B-ДНК берём 0.34 нм/п.о.
//! (Cantor & Schimmel 1980) и персистентную длину по Hagerman (1988).

use crate::fmath::DetMath;
pub const BP_NM: f64 = 0.34;

/// Персистентная длина дсДНК, нм.
pub const PERSISTENCE_NM: f64 = 50.0;

#[inline]
pub fn contour_length_nm(n_bp: f64) -> f64 {
    BP_NM * n_bp
}

/// Квадрат радиуса инерции, нм². Предел жёсткого стержня равен L²/12, длинного клубка — L·lp/3.
#[inline]
pub fn rg_sq_closed_nm2(n_bp: f64) -> f64 {
    debug_assert!(n_bp > 0.0, "размер популяции должен быть положительным");
    let lambda = contour_length_nm(n_bp) / PERSISTENCE_NM;
    let lp2 = PERSISTENCE_NM * PERSISTENCE_NM;
    lp2 * (lambda / 3.0 - 1.0 + 2.0 / lambda - 2.0 * (1.0 - (-lambda).dexp()) / (lambda * lambda))
}

#[inline]
pub fn rg_nm(n_bp: f64) -> f64 {
    rg_sq_closed_nm2(n_bp).sqrt()
}

/// Независимая квадратура Симпсона для проверки аналитического Rg².
/// Число интервалов растёт с числом персистентных длин в цепи.
pub fn rg_sq_simpson_nm2(n_bp: f64) -> f64 {
    debug_assert!(n_bp > 0.0);
    let l = contour_length_nm(n_bp);
    let lambda = l / PERSISTENCE_NM;
    let want = (8.0 * lambda).ceil() as usize;
    let n_int = if want <= 2048 {
        2048
    } else {
        want.div_ceil(2) * 2
    };
    let r2 = |s: f64| {
        2.0 * PERSISTENCE_NM * s
            - 2.0 * PERSISTENCE_NM * PERSISTENCE_NM * (1.0 - (-s / PERSISTENCE_NM).dexp())
    };
    let h = l / n_int as f64;
    let integrand = |s: f64| (l - s) * r2(s);
    let mut sum = integrand(0.0) + integrand(l);
    for i in 1..n_int {
        let s = h * i as f64;
        let w = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += w * integrand(s);
    }
    sum * h / 3.0 / (l * l)
}

/// Табличные пары (длина в п.о., Rg в нм) для PHYSICS.md §2, округлённые до 0.001 нм.
pub const RG_REFERENCE_TABLE: [(f64, f64); 8] = [
    (100.0, 9.197),
    (300.0, 24.612),
    (500.0, 37.227),
    (1000.0, 61.595),
    (2000.0, 95.781),
    (5000.0, 161.171),
    (10000.0, 232.893),
    (23000.0, 357.582),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wlc_limits_converge() {
        // Стержень: ρ = Rg²/(L_c²/12) → 1 как 1 − λ/5 (+λ²/30).
        // Асимптотический коэффициент 1/5 подтверждён численно (dev/λ → 0.199).
        // Допуск λ/4 покрывает старшие члены ряда; монотонное приближение к 1
        // обязательно при уменьшении N.
        let mut prev_dev = f64::INFINITY;
        for &n in &[200.0, 100.0, 50.0, 25.0, 10.0, 5.0] {
            let l = contour_length_nm(n);
            let lambda = l / PERSISTENCE_NM;
            let dev = (rg_sq_closed_nm2(n) / (l * l / 12.0) - 1.0).abs();
            assert!(
                dev <= lambda / 4.0,
                "стержень N={n}: отклонение {dev:.3e} > допуска {:.3e}",
                lambda / 4.0
            );
            assert!(
                dev <= prev_dev,
                "стержень N={n}: сходимость не монотонна ({dev} против {prev_dev})"
            );
            prev_dev = dev;
        }
        // Клубок: ρ = Rg²/(L_c·lp/3) → 1 как 1 − 3/λ (+O(λ⁻²)).
        // Асимптотика подтверждена численно (dev·λ → 3.000); допуск 3.5/λ.
        // Монотонное приближение при увеличении N.
        let mut prev_dev = f64::INFINITY;
        for &n in &[5000.0, 50000.0, 500000.0, 5000000.0] {
            let l = contour_length_nm(n);
            let lambda = l / PERSISTENCE_NM;
            let dev = (rg_sq_closed_nm2(n) / (l * PERSISTENCE_NM / 3.0) - 1.0).abs();
            assert!(
                dev <= 3.5 / lambda,
                "клубок N={n}: отклонение {dev:.3e} > допуска {:.3e}",
                3.5 / lambda
            );
            assert!(
                dev <= prev_dev,
                "клубок N={n}: сходимость не монотонна ({dev} против {prev_dev})"
            );
            prev_dev = dev;
        }
    }

    #[test]
    fn rg_closed_matches_simpson() {
        // База 2048 интервалов + адаптив 8 узлов на lp держит расхождение
        // замкнутой формы и интеграла ниже 1e-9 во всём диапазоне, включая
        // середину (N≈100–4000), где четвёртая производная максимальна.
        // Проверяем весь диапазон, а не только края.
        for n in [
            10.0, 50.0, 100.0, 250.0, 300.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0,
            6000.0, 8000.0, 15000.0, 23000.0, 40000.0, 60000.0, 100000.0,
        ] {
            let a = rg_sq_closed_nm2(n);
            let b = rg_sq_simpson_nm2(n);
            assert!(
                (a - b).abs() / a < 1e-9,
                "N={n}: замкнутая форма {a:.10} против Симпсона {b:.10}"
            );
        }
    }

    #[test]
    fn rg_reference_table_matches_engine() {
        // Правило документ↔код: опорные значения PHYSICS.md §2 равны движку.
        // Допуск АБСОЛЮТНЫЙ 5e-4 нм — по разрядности печати таблицы
        // (ревью №M2: относительный 0.5% был в 250 раз слабее реального
        // согласия и пропускал заниженные значения таблицы).
        for &(n, rg_doc) in &RG_REFERENCE_TABLE {
            let rg_engine = rg_nm(n);
            assert!(
                (rg_engine - rg_doc).abs() <= 5e-4,
                "Rg({n}): движок {rg_engine:.5}, документ {rg_doc}"
            );
        }
    }
}
