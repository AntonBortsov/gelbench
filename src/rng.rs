//! Один поток xoshiro256** на сценарий с посевом SplitMix64 (Blackman & Vigna 2019).
//! Сериализация состояния сохраняет последовательность при возобновлении прогона.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut sm = SplitMix64 { state: seed };
        let mut s = [0u64; 4];
        for slot in s.iter_mut() {
            *slot = sm.next_u64();
        }
        // Нулевой вектор — поглощающее состояние xoshiro.
        if s == [0, 0, 0, 0] {
            s[0] = 0x9E3779B97F4A7C15;
        }
        Rng { s }
    }

    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// Равномерное [0; 1) с 53 битами мантиссы.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

/// Инициализация четырёх слов xoshiro из одного сида (SplitMix64, Vigna 2015).
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_same_sequence() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..1000 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn different_seeds_diverge_immediately() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(43);
        assert_ne!(a.next_u64(), b.next_u64());
    }

    #[test]
    fn serde_roundtrip_continues_identically() {
        // Пауза/восстановление §3.6: сериализованное состояние обязано
        // продолжать ту же последовательность бит-в-бит.
        let mut live = Rng::new(7);
        live.next_u64();
        live.next_f64();
        let saved = serde_json::to_string(&live).expect("сериализация состояния ГСЧ");
        let mut restored: Rng = serde_json::from_str(&saved).unwrap();
        for _ in 0..100 {
            assert_eq!(live.next_u64(), restored.next_u64());
        }
    }
}
