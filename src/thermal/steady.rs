//! Steady transverse temperature rise under face cooling and lateral edge losses.

use super::{EdgeKind, SteadyProfileError, ThermalConfig, K_GEL};

impl ThermalConfig {
    /// Решает k·ΔT″ − s·ΔT + q_v = 0 прогонкой Томаса. Возвращает профиль перегрева и самосогласованный h.
    /// На боковых краях задаётся температура среды или условие Робина первого порядка.
    pub fn steady_profile(&self) -> Result<(Vec<f64>, f64), SteadyProfileError> {
        let n = self.n_points.max(3);
        let dy = self.width_m / (n - 1) as f64;
        let h = self.solve_h_selfconsistent();
        let s = self.sink_coeff(h);
        let q_v0 = self.q_v();
        let k_dy2 = K_GEL / (dy * dy);
        let k_dy = K_GEL / dy;

        // При обратной связи пересчитываем проводимость по среднему перегреву профиля.
        let mut kappa_scale = 1.0f64;
        let mut iters = 0usize;
        loop {
            let q_v = q_v0 * kappa_scale;
            let mut a = vec![0.0; n];
            let mut b = vec![0.0; n];
            let mut c = vec![0.0; n];
            let mut d = vec![0.0; n];
            for i in 1..n - 1 {
                a[i] = -k_dy2;
                b[i] = 2.0 * k_dy2 + s;
                c[i] = -k_dy2;
                d[i] = q_v; // Ta нормирована к нулю; диагональ приведена положительной
            }
            match self.edge {
                EdgeKind::Pinned => {
                    b[0] = 1.0;
                    c[0] = 0.0;
                    d[0] = 0.0;
                    a[n - 1] = 0.0;
                    b[n - 1] = 1.0;
                    d[n - 1] = 0.0;
                }
                EdgeKind::Convective => {
                    b[0] = k_dy + h;
                    c[0] = -k_dy;
                    d[0] = 0.0;
                    a[n - 1] = -k_dy;
                    b[n - 1] = k_dy + h;
                    d[n - 1] = 0.0;
                }
            }

            let mut cp = vec![0.0; n];
            let mut dp = vec![0.0; n];
            cp[0] = c[0] / b[0];
            dp[0] = d[0] / b[0];
            for i in 1..n {
                let m = b[i] - a[i] * cp[i - 1];
                cp[i] = c[i] / m;
                dp[i] = (d[i] - a[i] * dp[i - 1]) / m;
            }
            let mut t = vec![0.0; n];
            t[n - 1] = dp[n - 1];
            for i in (0..n - 1).rev() {
                t[i] = dp[i] - cp[i] * t[i + 1];
            }
            if !self.kappa_feedback {
                break Ok((t, h));
            }
            let t_mean = t.iter().sum::<f64>() / n as f64;
            let new_scale = 1.0 + crate::circuit::K_KAPPA_PER_C * t_mean;
            if (new_scale - kappa_scale).abs() / new_scale < 1e-12 {
                break Ok((t, h));
            }
            let loop_gain = crate::circuit::K_KAPPA_PER_C * t_mean / kappa_scale;
            if !loop_gain.is_finite() || !new_scale.is_finite() {
                break Err(SteadyProfileError::NonFiniteKappaFeedback { iteration: iters });
            }
            if loop_gain >= 1.0 {
                break Err(SteadyProfileError::KappaFeedbackRunaway { loop_gain });
            }
            iters += 1;
            if iters > 10_000 {
                break Err(SteadyProfileError::IterationLimit { iterations: iters });
            }
            kappa_scale = new_scale;
        }
    }

    /// Фиксированный сток с торцами при температуре среды. При нулевом стоке пик равен q_v·W²/(8k).
    pub fn steady_profile_explicit_sink_pinned(&self, s_sink: f64) -> Vec<f64> {
        let n = self.n_points.max(3);
        let dy = self.width_m / (n - 1) as f64;
        let q_v = self.q_v();
        let k_dy2 = K_GEL / (dy * dy);
        let mut b = vec![0.0; n];
        let mut c = vec![0.0; n];
        let mut d = vec![0.0; n];
        b[0] = 1.0;
        c[0] = 0.0;
        d[0] = 0.0;
        for i in 1..n - 1 {
            c[i] = -k_dy2;
            b[i] = 2.0 * k_dy2 + s_sink;
            d[i] = q_v;
        }
        c[n - 1] = 0.0;
        b[n - 1] = 1.0;
        d[n - 1] = 0.0;
        let mut cp = vec![0.0; n];
        let mut dp = vec![0.0; n];
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];
        for i in 1..n {
            let aa = if i == n - 1 { 0.0 } else { -k_dy2 };
            let m = b[i] - aa * cp[i - 1];
            cp[i] = c[i] / m;
            dp[i] = (d[i] - aa * dp[i - 1]) / m;
        }
        let mut t = vec![0.0; n];
        t[n - 1] = dp[n - 1];
        for i in (0..n - 1).rev() {
            t[i] = dp[i] - cp[i] * t[i + 1];
        }
        t
    }
}
