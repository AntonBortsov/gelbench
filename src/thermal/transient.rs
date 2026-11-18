//! Implicit transverse heat diffusion and the time constant for approach to equilibrium.
//! Temperatures are in K, time in s and volumetric power density in W/m³.

use super::{EdgeKind, ThermalConfig, K_GEL, RHO_CP};

impl ThermalConfig {
    /// Шаг Эйлера назад для перегрева относительно температуры среды.
    ///
    /// ```text
    /// ρ·c_p·(Tⁿ⁺¹−Tⁿ)/dt = k_g·T''ⁿ⁺¹ − s·(Tⁿ⁺¹−T_amb) + q_v,  s = 2h/t_g
    /// ```
    ///
    /// Теплоотдачу и мощность задаёт вызывающая сторона на каждом макрошаге.
    pub fn transient_step(
        &self,
        t_prev: &[f64],
        dt_s: f64,
        h_w_m2k: f64,
        q_v_w_m3: f64,
    ) -> Vec<f64> {
        let n = t_prev.len().max(3);
        let dy = self.width_m / (n - 1) as f64;
        let s = 2.0 * h_w_m2k / self.t_gel_m;
        let k_dy2 = K_GEL / (dy * dy);
        let k_dy = K_GEL / dy;
        let rc_dt = RHO_CP / dt_s;

        let mut a = vec![0.0; n];
        let mut b = vec![0.0; n];
        let mut c = vec![0.0; n];
        let mut d = vec![0.0; n];
        for i in 1..n - 1 {
            a[i] = -k_dy2;
            b[i] = rc_dt + 2.0 * k_dy2 + s;
            c[i] = -k_dy2;
            d[i] = rc_dt * t_prev[i] + q_v_w_m3;
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
                b[0] = k_dy + h_w_m2k;
                c[0] = -k_dy;
                d[0] = 0.0;
                a[n - 1] = -k_dy;
                b[n - 1] = k_dy + h_w_m2k;
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
        t
    }

    /// Возвращает конечный и усреднённый по времени профили холодного старта и h.
    /// Теплоотдача заморожена на самосогласованном стационарном значении.
    pub fn transient_run(&self, t_final_s: f64, dt_s: f64) -> (Vec<f64>, Vec<f64>, f64) {
        let h = self.solve_h_selfconsistent();
        let q_v0 = self.q_v();
        let n = self.n_points.max(3);
        let n_steps = ((t_final_s / dt_s).ceil() as usize).max(1);
        let dt = t_final_s / n_steps as f64;
        let mut t = vec![0.0; n];
        let mut mean = vec![0.0; n];
        for _ in 0..n_steps {
            let q_v = if self.kappa_feedback {
                let t_mean = t.iter().sum::<f64>() / n as f64;
                q_v0 * (1.0 + crate::circuit::K_KAPPA_PER_C * t_mean)
            } else {
                q_v0
            };
            t = self.transient_step(&t, dt, h, q_v);
            for (m, v) in mean.iter_mut().zip(t.iter()) {
                *m += v;
            }
        }
        for m in mean.iter_mut() {
            *m /= n_steps as f64;
        }
        (t, mean, h)
    }

    /// Постоянная времени прогрева τ = ρ·c_p·t_g/(2h), с.
    pub fn heating_time_constant_s(&self) -> f64 {
        let h = self.solve_h_selfconsistent();
        RHO_CP * self.t_gel_m / (2.0 * h)
    }
}
