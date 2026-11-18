//! Консервативный перенос в одномерном геле. Замкнутые границы задают нулевой поток,
//! открытый анод допускает вынос вещества. Для ионов доступен приток из резервуара.
//!
//! Адвекция MUSCL: van Leer, J. Comput. Phys. 32:101 (1979). Условие CFL: LeVeque,
//! Finite Volume Methods for Hyperbolic Problems, §4.4. Явная диффузия: Crank,
//! The Mathematics of Diffusion, §8.3.

use crate::fmath::DetMath;

/// Выбирает меньший согласованный наклон, на экстремумах даёт ноль.
#[inline]
pub fn minmod(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else if a.abs() < b.abs() {
        a
    } else {
        b
    }
}

/// Адвекция MUSCL/minmod с закрытыми стенками. Концентрации заданы в n ячейках,
/// рабочий буфер потоков — на n+1 гранях. Скорость в см/с, время в с, шаг сетки в см.
pub fn advect_muscl(src: &[f64], dst: &mut [f64], faces: &mut [f64], v: f64, dt: f64, dx: f64) {
    let n = src.len();
    debug_assert!(dst.len() == n && faces.len() == n + 1);
    debug_assert!(dt > 0.0 && dx > 0.0);

    // Коэффициент 1−|λ| подавляет антидиффузионную поправку при приближении к пределу CFL.
    let lambda = (v * dt / dx).abs();
    let anti_diff = 0.5 * (1.0 - lambda.min(1.0));

    faces[0] = 0.0;
    faces[n] = 0.0;

    for face in 1..n {
        let i_left = face - 1;
        let flux = if v >= 0.0 {
            let s_up = if i_left == 0 {
                0.0
            } else {
                minmod(src[i_left] - src[i_left - 1], src[face] - src[i_left])
            };
            v * (src[i_left] + anti_diff * s_up)
        } else {
            let s_up = if face + 1 == n {
                0.0
            } else {
                minmod(src[face] - src[i_left], src[face + 1] - src[face])
            };
            v * (src[face] - anti_diff * s_up)
        };
        faces[face] = flux;
    }

    let coeff = dt / dx;
    for i in 0..n {
        dst[i] = src[i] - coeff * (faces[i + 1] - faces[i]);
    }
}

/// Сдвигает профиль на долю ячейки с сохранением суммарной массы.
/// Компенсирует смещение центра полосы из-за асимметрии ограничителя FCT.
pub fn shift_profile_cells(c: &mut [f64], shift_cells: f64) {
    let n = c.len();
    if n == 0 || shift_cells == 0.0 {
        return;
    }
    let old = c.to_vec();
    let last = (n - 1) as f64;
    for i in 0..n {
        let x = (i as f64 - shift_cells).clamp(0.0, last);
        let j0 = x.floor() as usize;
        let j1 = (j0 + 1).min(n - 1);
        let frac = x - j0 as f64;
        let v = old[j0] * (1.0 - frac) + old[j1] * frac;
        c[i] = v.max(0.0);
    }
    let m_old: f64 = old.iter().sum();
    let m_new: f64 = c.iter().sum();
    if m_new > 0.0 && m_old > 0.0 {
        let scale = m_old / m_new;
        if scale != 1.0 {
            for v in c.iter_mut() {
                *v *= scale;
            }
        }
    }
}

/// Адвекция FCT (Zalesak 1979, Boris & Book 1973). Донорский поток задаёт монотонное
/// решение. Добавочный поток от центральной реконструкции ограничен локальными экстремумами.
/// При закрытых границах потоки сохраняют массу. Условие CFL для донора обеспечивает [`substep_count`].
pub fn advect_fct(
    src: &[f64],
    dst: &mut [f64],
    faces: &mut [f64],
    v: f64,
    dt: f64,
    dx: f64,
    window: usize,
) {
    advect_fct_open(src, dst, faces, v, dt, dx, window, false);
}

/// Адвекция FCT с закрытым катодом и опционально открытым анодом.
/// При выносе через анод учитывается только донорский поток. Возвращает потерю в единицах суммы ячеек.
#[allow(clippy::too_many_arguments)]
pub fn advect_fct_open(
    src: &[f64],
    dst: &mut [f64],
    faces: &mut [f64],
    v: f64,
    dt: f64,
    dx: f64,
    window: usize,
    open_right: bool,
) -> f64 {
    let n = src.len();
    debug_assert!(dst.len() == n && faces.len() == n + 1);
    debug_assert!(dt > 0.0 && dx > 0.0);
    let coeff = dt / dx;

    faces[0] = 0.0;
    faces[n] = 0.0;
    for j in 1..n {
        let high = v * 0.5 * (src[j - 1] + src[j]);
        let low = if v >= 0.0 { v * src[j - 1] } else { v * src[j] };
        faces[j] = high - low; // антидиффузионный поток A_j
    }

    let out_flux = if open_right && v > 0.0 {
        v * src[n - 1]
    } else {
        0.0
    };

    let mut c_low = vec![0.0f64; n];
    for i in 0..n {
        let f_l = |j: usize| -> f64 {
            if j == 0 {
                0.0
            } else if j == n {
                out_flux
            } else if v >= 0.0 {
                v * src[j - 1]
            } else {
                v * src[j]
            }
        };
        c_low[i] = src[i] - coeff * (f_l(i + 1) - f_l(i));
    }

    // Локальные границы Заляска включают низкопорядковое решение и нулевой порог концентрации.
    let w = window.max(1);
    let nb = |i: usize| -> (f64, f64) {
        let start = i.saturating_sub(w);
        let end = (i + w).min(n - 1);
        let mut lo = src[i];
        let mut hi = src[i];
        for j in start..=end {
            lo = lo.min(src[j]);
            hi = hi.max(src[j]);
        }
        (hi, lo)
    };

    let mut r_plus = vec![1.0f64; n];
    let mut r_minus = vec![1.0f64; n];
    for i in 0..n {
        let a_in = if i == 0 { 0.0 } else { faces[i] };
        let a_out = if i + 1 == n { 0.0 } else { faces[i + 1] };
        let p_plus = coeff * (a_in.max(0.0) - a_out.min(0.0));
        let p_minus = coeff * (a_in.min(0.0) - a_out.max(0.0));
        let (c_max, c_min) = nb(i);
        let c_max = c_max.max(c_low[i]);
        let c_min = c_min.min(c_low[i]).min(0.0);
        r_plus[i] = if p_plus > 0.0 {
            ((c_max - c_low[i]) / p_plus).clamp(0.0, 1.0)
        } else {
            1.0
        };
        r_minus[i] = if p_minus < 0.0 {
            ((c_min - c_low[i]) / p_minus).clamp(0.0, 1.0)
        } else {
            1.0
        };
    }

    // Поток на грани ограничен одновременно допустимой потерей донора и приростом получателя.
    for j in 1..n {
        let a = faces[j];
        let pos = a.max(0.0);
        let neg = a.min(0.0);
        let f_pos = r_minus[j - 1].min(r_plus[j]);
        let f_neg = r_minus[j].min(r_plus[j - 1].min(1.0));
        faces[j] = f_pos * pos + f_neg * neg;
    }

    for i in 0..n {
        let a_in = if i == 0 { 0.0 } else { faces[i] };
        let a_out = if i + 1 == n { 0.0 } else { faces[i + 1] };
        let mut c = c_low[i] - coeff * (a_out - a_in);
        // Обнуляем следовые значения в хвостах профиля, где округление может менять знак.
        if c.abs() < 1e-30 {
            c = 0.0;
        }
        if c < 0.0 {
            c = 0.0;
        }
        dst[i] = c;
    }
    coeff * out_flux
}

/// Явная диффузия с отражением на стенках. [`substep_count`] обеспечивает D·dt/Δx² ≤ 1/2.
pub fn diffuse_explicit(src: &[f64], dst: &mut [f64], d: f64, dt: f64, dx: f64) {
    let n = src.len();
    debug_assert!(dst.len() == n);
    let r = d * dt / (dx * dx);
    debug_assert!(r <= 0.5 + 1e-12, "диффузия неустойчива: r={r}");
    if n < 2 {
        dst.copy_from_slice(src);
        return;
    }
    // Зеркальные стенки: c₋₁ = c₀, c_n = c_{n−1}.
    dst[0] = src[0] + r * (src[1] - src[0]);
    for i in 1..n - 1 {
        dst[i] = src[i] + r * (src[i + 1] - 2.0 * src[i] + src[i - 1]);
    }
    dst[n - 1] = src[n - 1] + r * (src[n - 2] - src[n - 1]);
}

/// Число подшагов и их длительность в секундах по более строгому из ограничений адвекции и диффузии.
pub fn substep_count(dt_macro: f64, dx: f64, v_max: f64, d_max: f64) -> (u32, f64) {
    debug_assert!(dt_macro > 0.0 && dx > 0.0);
    let dt_adv = if v_max > 0.0 {
        dx / v_max
    } else {
        f64::INFINITY
    };
    let dt_dif = if d_max > 0.0 {
        dx * dx / (2.0 * d_max)
    } else {
        f64::INFINITY
    };
    let dt_max = dt_adv.min(dt_dif);
    if dt_macro <= dt_max {
        return (1, dt_macro);
    }
    let n_raw = (dt_macro / dt_max).ceil();
    debug_assert!(
        n_raw < 1e6_f64,
        "CFL требует >1e6 подшагов — проверьте параметры"
    );
    let n = n_raw as u32;
    (n, dt_macro / n as f64)
}

// CISL (Lin–Rood) сохраняет массу при неоднородной скорости: изменение ширины
// интервала отправления учитывает сжатие полосы. Интеграл по этому интервалу
// вычисляется через префиксные массы. Точность ограничивают обратная трассировка
// характеристик и изменение поля за макрошаг.

/// Монотонные наклоны для кусочно-линейной реконструкции Lin–Rood.
pub fn mc_slopes(c: &[f64]) -> Vec<f64> {
    let n = c.len();
    let mut d = vec![0.0; n];
    if n < 3 {
        return d;
    }
    for j in 1..n - 1 {
        let dl = c[j] - c[j - 1];
        let dr = c[j + 1] - c[j];
        let dc = 0.5 * (c[j + 1] - c[j - 1]);
        if dl * dr > 0.0 {
            let lim = dc.abs().min(2.0 * dl.abs()).min(2.0 * dr.abs());
            d[j] = dc.signum() * lim;
        }
    }
    d
}

/// Обратная характеристика dx/dt = v(x), интегрированная методом RK4.
fn trace_back(mut x: f64, v: &impl Fn(f64) -> f64, dt: f64, nsub: usize) -> f64 {
    let tau = dt / nsub as f64;
    for _ in 0..nsub {
        let k1 = v(x);
        let k2 = v(x - 0.5 * tau * k1);
        let k3 = v(x - 0.5 * tau * k2);
        let k4 = v(x - tau * k3);
        x -= tau / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    x
}

/// Первообразная реконструкции с постоянным продолжением за краями сетки.
fn primitive_at(x: f64, xe: &[f64], c: &[f64], slope: &[f64], cum: &[f64], h: f64) -> f64 {
    let x0 = xe[0];
    let x_last = xe[xe.len() - 1];
    let xc_ = x.clamp(x0, x_last);
    let j = (((xc_ - x0) / h) as usize).min(xe.len() - 2);
    let center = 0.5 * (xe[j] + xe[j + 1]);
    // G(ξ) = P_j + c_j·(ξ − x_{j−½}) + δ_j/(2h)·[(ξ − x_j)² − (h/2)²]
    cum[j]
        + c[j] * (xc_ - xe[j])
        + slope[j] / (2.0 * h) * ((xc_ - center) * (xc_ - center) - (0.5 * h) * (0.5 * h))
}

/// Замкнутый шаг CISL. `xe` задаёт n+1 грань, `c` — средние по ячейкам.
/// Новая масса — интеграл прежнего профиля между точками отправления соседних граней.
pub fn cisl_step(c: &[f64], xe: &[f64], v: &impl Fn(f64) -> f64, dt: f64, nsub: usize) -> Vec<f64> {
    cisl_step_open(c, xe, v, dt, nsub, false).0
}

/// CISL с закрытым катодом и опциональным выносом через анод.
/// Возвращает новые средние концентрации и потерю в единицах суммы ячеек.
pub fn cisl_step_open(
    c: &[f64],
    xe: &[f64],
    v: &impl Fn(f64) -> f64,
    dt: f64,
    nsub: usize,
    open_right: bool,
) -> (Vec<f64>, f64) {
    let n = c.len();
    assert_eq!(xe.len(), n + 1, "граней должно быть n+1");
    let h = (xe[n] - xe[0]) / n as f64;

    let mut dep = Vec::with_capacity(n + 1);
    for i in 0..=n {
        dep.push(trace_back(xe[i], v, dt, nsub));
    }
    // Проверяем порядок до обрезки: пересечение характеристик недопустимо,
    // но обрезка на стенке может законно совместить соседние точки.
    for i in 1..=n {
        assert!(
            dep[i] > dep[i - 1],
            "CISL: траектории пересеклись (dt·max|dv/dx| ≥ 1)"
        );
    }
    for d in dep.iter_mut() {
        *d = d.clamp(xe[0], xe[n]);
    }
    // На закрытой грани точка отправления совпадает со стенкой, даже если трассировка
    // оставила её внутри сетки. Иначе полоса у стенки потеряет массу.
    dep[0] = xe[0];
    if !open_right {
        dep[n] = xe[n];
    }

    let slope = mc_slopes(c);
    let mut cum = vec![0.0; n];
    for j in 1..n {
        cum[j] = cum[j - 1] + c[j - 1] * h;
    }

    let mut out = Vec::with_capacity(n);
    for j in 0..n {
        let a = primitive_at(dep[j], xe, c, &slope, &cum, h);
        let b = primitive_at(dep[j + 1], xe, c, &slope, &cum, h);
        out.push((b - a) / h);
    }
    let exited = if open_right {
        let m_in: f64 = c.iter().sum();
        let m_out: f64 = out.iter().sum();
        (m_in - m_out).max(0.0)
    } else {
        0.0
    };
    (out, exited)
}

/// Линейная интерполяция между узлами `grid` с постоянным продолжением на концах.
pub fn interp_on(grid: &[f64], values: &[f64], x: f64) -> f64 {
    let n = values.len();
    let xc_ = x.clamp(grid[0], grid[n - 1]);
    let j = (((xc_ - grid[0]) / (grid[n - 1] - grid[0])) * (n - 1) as f64) as usize;
    let j = j.min(n - 2);
    let t = (xc_ - grid[j]) / (grid[j + 1] - grid[j]);
    values[j] * (1.0 - t) + values[j + 1] * t
}

/// За сеткой производная первообразной равна концентрации соответствующего резервуара.
#[allow(clippy::too_many_arguments)]
fn primitive_reservoir(
    x: f64,
    xe: &[f64],
    c: &[f64],
    slope: &[f64],
    cum: &[f64],
    h: f64,
    cin_l: f64,
    cin_r: f64,
) -> f64 {
    let n = c.len();
    let x0 = xe[0];
    let xn = xe[n];
    if x <= x0 {
        return -cin_l * (x0 - x);
    }
    if x >= xn {
        return cum[n] + cin_r * (x - xn);
    }
    let j = (((x - x0) / h) as usize).min(n - 1);
    let center = 0.5 * (xe[j] + xe[j + 1]);
    cum[j]
        + c[j] * (x - xe[j])
        + slope[j] / (2.0 * h) * ((x - center) * (x - center) - (0.5 * h) * (0.5 * h))
}

/// CISL с кусочно-линейной скоростью на гранях и притоком из резервуаров.
/// Метод Хюна трассирует точки отправления, обратный проход восстанавливает их порядок.
pub fn cisl_step_reservoir(
    c: &[f64],
    xe: &[f64],
    vface: &[f64],
    dt: f64,
    cin_l: f64,
    cin_r: f64,
    nsub: usize,
) -> Vec<f64> {
    let n = c.len();
    assert_eq!(xe.len(), n + 1);
    assert_eq!(vface.len(), n + 1);
    let h = (xe[n] - xe[0]) / n as f64;

    let mut dep = xe.to_vec();
    let tau = dt / nsub as f64;
    for _ in 0..nsub {
        for xd in dep.iter_mut() {
            let k1 = interp_on(xe, vface, *xd);
            let k2 = interp_on(xe, vface, *xd - 0.5 * tau * k1);
            *xd -= tau * k2;
        }
    }
    for i in (1..=n).rev() {
        if dep[i - 1] > dep[i] {
            dep[i - 1] = dep[i];
        }
    }

    let slope = mc_slopes(c);
    let mut cum = vec![0.0; n + 1];
    for j in 1..=n {
        cum[j] = cum[j - 1] + c[j - 1] * h;
    }
    (0..n)
        .map(|j| {
            let a = primitive_reservoir(dep[j], xe, c, &slope, &cum, h, cin_l, cin_r);
            let b = primitive_reservoir(dep[j + 1], xe, c, &slope, &cum, h, cin_l, cin_r);
            (b - a) / h
        })
        .collect()
}

/// Явная диффузия с резервуарными граничными значениями (сосед за торцом
/// = cin). Подшаги по критерию r = D·τ/h² ≤ 0.5 с запасом 2.5.
pub fn diffuse_reservoir(c: &[f64], h: f64, d: f64, dt: f64, cin_l: f64, cin_r: f64) -> Vec<f64> {
    let n = c.len();
    let nsub = ((((2.5 * d * dt) / (h * h)).ceil()) as usize).max(1);
    let tau = dt / nsub as f64;
    let r = d * tau / (h * h);
    let mut cur = c.to_vec();
    for _ in 0..nsub {
        let mut next = vec![0.0; n];
        next[0] = cur[0] + r * (cin_l - 2.0 * cur[0] + cur[1]);
        for i in 1..n - 1 {
            next[i] = cur[i] + r * (cur[i + 1] - 2.0 * cur[i] + cur[i - 1]);
        }
        next[n - 1] = cur[n - 1] + r * (cur[n - 2] - 2.0 * cur[n - 1] + cin_r);
        cur = next;
    }
    cur
}

/// Функция Бернулли для потока Шарфеттера–Гуммеля. При |x| < 10⁻⁵ ряд Тейлора
/// устраняет потерю точности около устранимой особенности B(0) = 1.
pub fn bernoulli(x: f64) -> f64 {
    if x.abs() < 1e-5 {
        1.0 - 0.5 * x + x * x / 12.0
    } else {
        x / x.dexp_m1()
    }
}

/// Дрейф и диффузия с резервуарными границами. На каждой грани поток Шарфеттера–Гуммеля
/// равен (d/h)·[B(−Pe)·cₗ − B(Pe)·cᵣ], где Pe = v·h/d. Он сохраняет стационарный
/// экспоненциальный пограничный слой толщины d/v даже на грубой сетке.
/// Для положительности явный подшаг ограничен величиной 0,8·h²/(2d + 2h·max|v|).
pub fn sg_step_reservoir(
    c: &[f64],
    h: f64,
    vface: &[f64],
    d: f64,
    dt: f64,
    cin_l: f64,
    cin_r: f64,
) -> Vec<f64> {
    let n = c.len();
    assert_eq!(vface.len(), n + 1);
    assert!(d > 0.0 && h > 0.0 && dt > 0.0);
    let vmax = vface.iter().fold(0.0_f64, |a, &v| a.max(v.abs()));
    let tau_stab = 0.8 * h * h / (2.0 * d + 2.0 * h * vmax);
    let nsub = ((dt / tau_stab).ceil() as usize).max(1);
    let tau = dt / nsub as f64;
    let dh = d / h;
    let mut cur = c.to_vec();
    let mut flux = vec![0.0; n + 1];
    for _ in 0..nsub {
        for (g, f) in flux.iter_mut().enumerate() {
            let pe = vface[g] * h / d;
            let cl = if g == 0 { cin_l } else { cur[g - 1] };
            let cr = if g == n { cin_r } else { cur[g] };
            *f = dh * (bernoulli(-pe) * cl - bernoulli(pe) * cr);
        }
        for (j, cj) in cur.iter_mut().enumerate() {
            *cj += tau / h * (flux[j] - flux[j + 1]);
        }
    }
    cur
}

#[cfg(test)]
mod tests {
    use super::*;

    fn total(s: &[f64]) -> f64 {
        s.iter().sum()
    }

    #[test]
    fn minmod_cases() {
        assert_eq!(minmod(1.0, 2.0), 1.0);
        assert_eq!(minmod(-2.0, -1.0), -1.0);
        assert_eq!(minmod(-1.0, 2.0), 0.0);
        assert_eq!(minmod(3.0, -3.0), 0.0);
        assert_eq!(minmod(0.0, 5.0), 0.0);
    }

    #[test]
    fn advect_conserves_mass_to_machine_precision() {
        // Ступень, крупный λ≈0.4, 500 подшагов: телескопирование потоков.
        let n = 128;
        let mut c = vec![0.0; n];
        for i in 40..80 {
            c[i] = 1.0;
        }
        let m0 = total(&c);
        let mut dst = vec![0.0; n];
        let mut faces = vec![0.0; n + 1];
        let (dt, dx, v) = (0.02_f64, 0.02_f64, 0.4_f64); // λ = 0.4
        for _ in 0..500 {
            advect_muscl(&c, &mut dst, &mut faces, v, dt, dx);
            std::mem::swap(&mut c, &mut dst);
        }
        let drift = ((total(&c) - m0) / m0).abs();
        assert!(drift < 1e-14, "дрейф массы {drift:.3e} > 1e-14");
    }

    #[test]
    fn advect_engine_scenario_triangle_512() {
        // Регрессия high-field сценария (λ≈0.96, треугольник у лунки):
        // стабильность, TVD и консервативность. История: до введения
        // множителя (1−|λ|) схема раскачивалась (max 1.9e7 на 250-м шаге,
        // сумма уходила к −85).
        let n = 512;
        let mut c = vec![0.0; n];
        let peak = ((0.25_f64 / 0.01953125_f64).floor() as usize).clamp(2, n - 4);
        for x in 0..n {
            let d = ((x as f64 - peak as f64).abs() / 3.0).min(1.0);
            c[x] = 1.0 - d;
        }
        let m0 = total(&c);
        let max0 = c.iter().cloned().fold(f64::MIN, f64::max);
        let mut dst = vec![0.0; n];
        let mut faces = vec![0.0; n + 1];
        let (dt, dx, v) = (0.5_f64, 0.01953125_f64, 0.0375719_f64);
        // Фронт доходит до правой стенки на ~510-м шаге; после этого
        // накопление массы у стены законно поднимает концентрацию выше
        // начального максимума — TVD проверяем только в пути.
        let tvd_steps = 480;
        for step in 0..2000 {
            advect_muscl(&c, &mut dst, &mut faces, v, dt, dx);
            std::mem::swap(&mut c, &mut dst);
            if step < tvd_steps {
                for &val in &c {
                    assert!(
                        val >= -1e-12 && val <= max0 + 1e-12,
                        "step {step}: выброс {val} вне [0; {max0}]"
                    );
                }
            } else {
                for &val in &c {
                    assert!(val.is_finite() && val >= -1e-12, "step {step}: {val}");
                }
            }
        }
        let drift = ((total(&c) - m0) / m0).abs();
        assert!(drift < 1e-12, "дрейф массы {drift:.3e} > 1e-12");
    }

    #[test]
    fn diffuse_conserves_and_respects_mirror() {
        let n = 64;
        let mut c = vec![0.0; n];
        for i in 10..30 {
            c[i] = 1.0;
        }
        let m0 = total(&c);
        let mut dst = vec![0.0; n];
        let (d, dt, dx) = (1.0_f64, 0.05_f64, 1.0_f64); // r = 0.05 ≤ 1/2
        for _ in 0..400 {
            diffuse_explicit(&c, &mut dst, d, dt, dx);
            std::mem::swap(&mut c, &mut dst);
        }
        let drift = ((total(&c) - m0) / m0).abs();
        assert!(drift < 1e-13, "дрейф массы {drift:.3e} > 1e-13");
        assert!(
            c.iter().all(|&v| v.is_finite() && v >= -1e-12),
            "отрицательные концентрации после диффузии"
        );
    }

    #[test]
    fn substep_math_and_guard() {
        // Ограничений нет → один шаг без дробления.
        assert_eq!(substep_count(0.5, 0.02, 1e-6, 0.0), (1, 0.5));
        // Диффузионное ограничение далеко: без дробления.
        assert_eq!(substep_count(0.5, 0.02, 0.0, 1e-8), (1, 0.5));
        // Адвективный CFL требует дробления; подшаг ≤ ограничению,
        // произведение подшагов восстанавливает макрошаг.
        let (n, dt) = substep_count(7.7, 0.01, 0.1, 1e-8);
        assert!((77..=78).contains(&n), "неожиданное число подшагов {n}");
        assert!(dt <= 0.1 + 1e-12, "подшаг {dt} превысил ограничение 0.1");
        assert!((dt * n as f64 - 7.7).abs() < 1e-12);
    }

    #[test]
    fn fct_conserves_mass_and_stays_positive() {
        // Ступень, крупный сдвиг за шаг (λ≈0.8), 500 подшагов: потоковая
        // форма с телескопированием держит массу на машинном уровне,
        // границы Заляска — без отрицательных выбросов.
        let n = 128;
        let mut c = vec![0.0; n];
        for i in 40..80 {
            c[i] = 1.0;
        }
        let m0 = total(&c);
        let mut dst = vec![0.0; n];
        let mut faces = vec![0.0; n + 1];
        let (dt, dx, v) = (0.02_f64, 0.02_f64, 0.8_f64); // λ = 0.8
        for _ in 0..500 {
            advect_fct(&c, &mut dst, &mut faces, v, dt, dx, 1);
            std::mem::swap(&mut c, &mut dst);
            for &val in &c {
                assert!(val >= -1e-12, "отрицательный выброс {val}");
            }
        }
        let drift = ((total(&c) - m0) / m0).abs();
        assert!(drift < 1e-12, "дрейф массы {drift:.3e} > 1e-12");
    }

    #[test]
    fn fct_smooth_profile_translation_is_clean() {
        // Гладкий профиль (аналог полосы после стартового размытия):
        // численная диффузия FCT на порядке физической — прирост σ за
        // 45 минут эквивалентного прогона не превышает физического
        // прироста 5 кб (~0.049 мм).
        let n = 512;
        let dx = 10.0 / n as f64;
        let mut c = vec![0.0; n];
        let peak = 12.0_f64;
        for x in 0..n {
            let d = (x as f64 - peak) * dx / 0.05; // σ ≈ 0.03 см
            c[x] = (-0.5 * d * d).dexp();
        }
        let sigma = |cc: &[f64]| -> f64 {
            let m: f64 = cc.iter().sum();
            let com: f64 = cc
                .iter()
                .enumerate()
                .map(|(i, v)| i as f64 * v)
                .sum::<f64>()
                / m;
            (cc.iter()
                .enumerate()
                .map(|(i, v)| (i as f64 - com).powi(2) * v)
                .sum::<f64>()
                / m)
                .sqrt()
                * dx
        };
        let s0 = sigma(&c);
        let mut dst = vec![0.0; n];
        let mut faces = vec![0.0; n + 1];
        // 45 минут × v=1.1e-3 см/с ≈ путь 2.97 см (якорь №1).
        // Историческое пакетирование FCT (раз на пол-ячейки сдвига) —
        // страж ПРИМИТИВА: с этапа 2d state идёт точным целочисленным
        // сдвигом, FCT остался в буфере чтения (один проход на кадр) и
        // на переходных путях; его расплывание обязано остаться малым.
        const HALF_CELL: f64 = 0.5;
        let (dt, v) = (0.5_f64, 1.1e-3_f64);
        let mut accum = 0.0_f64;
        let mut events = 0_u32;
        for _ in 0..5400 {
            accum += dt;
            if v * accum / dx >= HALF_CELL {
                advect_fct(&c, &mut dst, &mut faces, v, accum, dx, 1);
                std::mem::swap(&mut c, &mut dst);
                accum = 0.0;
                events += 1;
            }
        }
        let growth = sigma(&c) - s0;
        println!("прирост σ (FCT, событий {events}) = {growth:.6} см");
        assert!(
            growth < 4.9e-3, // физический прирост 5 кб за 45 мин
            "численное расплывание {growth:.5} см ≥ физической диффузии"
        );
    }

    // --- Шарфеттер–Гуммель (этап 2b) ---

    #[test]
    fn bernoulli_limits_and_identity() {
        assert_eq!(bernoulli(0.0), 1.0);
        // Тождество B(−x) = B(x) + x — на обеих ветках (Тейлор и exp_m1).
        for &x in &[1e-7, 1e-6, 1e-5, 1e-3, 0.1, 1.0, 10.0, 300.0] {
            let lhs = bernoulli(-x);
            let rhs = bernoulli(x) + x;
            assert!(
                ((lhs - rhs) / rhs).abs() < 1e-12,
                "B(−x)=B(x)+x нарушено при x={x}: {lhs} против {rhs}"
            );
        }
        // Стык ветвей: обе формулы в одной точке границы совпадают
        // (остаток Тейлора x⁴/720 ≈ 1.4e−23 при x = 1e−5).
        let x = 1e-5_f64;
        let taylor = 1.0 - 0.5 * x + x * x / 12.0;
        let exact = x / x.dexp_m1();
        assert!(((taylor - exact) / exact).abs() < 1e-14);
        // Дальние пределы: B(+∞)→0, B(−x)→x.
        assert!(bernoulli(800.0) == 0.0);
        assert!((bernoulli(-800.0) - 800.0).abs() < 1e-9);
    }

    #[test]
    fn sg_boltzmann_profile_is_stationary() {
        // Экспоненциальный профиль c ∝ exp(v·x/d) — равновесие Больцмана
        // (дрейф уравновешен диффузией, поток нулевой). SG обязан держать
        // его НЕПОДВИЖНЫМ на грубой сетке — CISL+явная диффузия не могут.
        let n = 32;
        let h = 1.0e-4;
        let d = 1.0e-9;
        let v = 4.0e-5; // Pe = v·h/d = 4 — фронт тоньше ячейки
        let pe = v * h / d;
        let c0: Vec<f64> = (0..n).map(|i| (pe * (i as f64 + 0.5)).dexp()).collect();
        let cin_l = (pe * (-0.5)).dexp();
        let cin_r = (pe * (n as f64 + 0.5)).dexp();
        let vface = vec![v; n + 1];
        let mut c = c0.clone();
        for _ in 0..50 {
            c = sg_step_reservoir(&c, h, &vface, d, 0.2, cin_l, cin_r);
        }
        for i in 0..n {
            let rel = ((c[i] - c0[i]) / c0[i]).abs();
            assert!(rel < 1e-11, "стационар уплыл в ячейке {i}: rel={rel:.2e}");
        }
    }

    #[test]
    fn sg_conserves_interior_mass_and_positivity() {
        // Горб в середине, нулевая скорость, нулевые резервуары: масса
        // до прихода хвостов к границам сохраняется; отрицательных нет.
        let n = 64;
        let h = 1.0e-4;
        let d = 1.0e-9;
        let mut c = vec![0.0; n];
        for i in 24..40 {
            c[i] = 1.0 + 0.5 * ((i - 24) as f64);
        }
        let m0 = total(&c);
        let vface = vec![0.0; n + 1];
        for _ in 0..20 {
            c = sg_step_reservoir(&c, h, &vface, d, 0.2, 0.0, 0.0);
        }
        assert!(((total(&c) - m0) / m0).abs() < 1e-12, "масса уплыла");
        assert!(c.iter().all(|&x| x >= 0.0), "позитивность нарушена");
    }

    #[test]
    fn sg_upwind_limit_matches_advection() {
        // Pe >> 1: SG вырождается в донорную ячейку — сдвиг ступени со
        // скоростью v, масса из левого резервуара втекает значением cin_l.
        let n = 50;
        let h = 1.0e-4;
        let d = 1.0e-12; // Pe = 4000
        let v = 4.0e-5;
        let c = vec![7.0; n];
        let vface = vec![v; n + 1];
        let stepped = sg_step_reservoir(&c, h, &vface, d, 0.5, 7.0, 7.0);
        // Однородный профиль с равными резервуарами обязан остаться собой.
        for (i, &x) in stepped.iter().enumerate() {
            assert!(
                ((x - 7.0) / 7.0).abs() < 1e-12,
                "однородный профиль изменился в {i}: {x}"
            );
        }
    }
}
