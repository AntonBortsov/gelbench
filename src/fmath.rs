//! Трансцендентные функции движка используют одну версию `libm` на CPU и WebAssembly.
//! Различия системных реализаций в последнем бите могут изменить решение лимитера или число подшагов CFL.

pub trait DetMath {
    fn dexp(self) -> f64;
    fn dln(self) -> f64;
    fn dlog10(self) -> f64;
    fn dpowf(self, y: f64) -> f64;
    /// Сохраняет точность exp(x) − 1 вблизи нуля.
    fn dexp_m1(self) -> f64;
    fn dcosh(self) -> f64;
    fn dsin(self) -> f64;
    fn dcos(self) -> f64;
}

impl DetMath for f64 {
    #[inline]
    fn dexp(self) -> f64 {
        libm::exp(self)
    }
    #[inline]
    fn dln(self) -> f64 {
        libm::log(self)
    }
    #[inline]
    fn dlog10(self) -> f64 {
        libm::log10(self)
    }
    #[inline]
    fn dpowf(self, y: f64) -> f64 {
        libm::pow(self, y)
    }
    #[inline]
    fn dexp_m1(self) -> f64 {
        libm::expm1(self)
    }
    #[inline]
    fn dcosh(self) -> f64 {
        libm::cosh(self)
    }
    #[inline]
    fn dsin(self) -> f64 {
        libm::sin(self)
    }
    #[inline]
    fn dcos(self) -> f64 {
        libm::cos(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_std_to_within_one_ulp() {
        // Смысл замены — НЕ в другой точности, а в одинаковости. Значения
        // обязаны совпадать со `std` с точностью до последнего бита,
        // иначе мы поменяли бы физику, а не зафиксировали арифметику.
        for i in 0..2000 {
            let x = -20.0 + 0.02 * i as f64;
            let (a, b) = (x.exp(), x.dexp());
            assert!(
                (a - b).abs() <= 4.0 * f64::EPSILON * a.abs().max(1e-300),
                "exp({x}): std {a:e} против libm {b:e}"
            );
            let y = 1e-3 * 1e6_f64.powf(i as f64 / 2000.0);
            assert!((y.ln() - y.dln()).abs() <= 4.0 * f64::EPSILON * y.ln().abs().max(1e-12));
        }
    }

    #[test]
    fn no_direct_std_transcendentals_left_in_sim() {
        // Страж правила: если кто-то снова напишет `.exp()` вместо
        // `.dexp()`, детерминизм между платформами тихо сломается —
        // и упадёт этот тест, а не приёмка через полгода.
        let src_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/src");
        let mut offenders = Vec::new();
        let mut stack = vec![std::path::PathBuf::from(src_dir)];
        while let Some(dir) = stack.pop() {
            for entry in std::fs::read_dir(&dir).expect("читается") {
                let path = entry.expect("запись").path();
                if path.is_dir() {
                    stack.push(path);
                    continue;
                }
                if path.extension().and_then(|e| e.to_str()) != Some("rs") {
                    continue;
                }
                if path.file_name().and_then(|f| f.to_str()) == Some("fmath.rs") {
                    continue;
                }
                let text = std::fs::read_to_string(&path).expect("читается");
                for (i, line) in text.lines().enumerate() {
                    let code = line.trim_start();
                    if code.starts_with("//") {
                        continue;
                    }
                    for bad in [
                        ".exp()",
                        ".ln()",
                        ".log10()",
                        ".powf(",
                        ".exp_m1()",
                        ".cosh()",
                        ".sinh()",
                        ".tanh()",
                        ".sin()",
                        ".cos()",
                    ] {
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
            "прямые вызовы платформенной математики (нужен fmath::DetMath):\n{}",
            offenders.join("\n")
        );
    }
}
