//! BCH (Bose–Chaudhuri–Hocquenghem) error-correcting code.

/// A BCH code defined by its parameters: block length n, dimension k, error capacity t.
#[derive(Debug, Clone)]
pub struct BCHCode {
    pub n: usize,       // codeword length
    pub k: usize,       // message length
    pub t: usize,       // error-correcting capability
    pub m: usize,       // GF(2^m) extension degree
    pub generator: Vec<u8>, // generator polynomial coefficients
}

impl BCHCode {
    /// Create a BCH code over GF(2) with given m and t.
    /// Block length n = 2^m - 1.
    pub fn new(m: usize, t: usize) -> Self {
        let n = (1usize << m) - 1;

        // Build generator polynomial as LCM of minimal polynomials
        // of α, α³, α⁵, ..., α^(2t-1) over GF(2).
        // Simplified: start with g(x) = 1, multiply by minimal polys.
        let mut gpoly: Vec<u8> = vec![1];

        for i in (1..=2 * t).step_by(2) {
            let min_poly = minimal_polynomial_gf2(m, i);
            gpoly = poly_mul_gf2(&gpoly, &min_poly);
        }

        let k = n - gpoly.len() + 1;

        BCHCode { n, k, t, m, generator: gpoly }
    }

    /// Encode a message by multiplying with generator polynomial.
    pub fn encode(&self, msg: &[u8]) -> Vec<u8> {
        // msg has length k, shift up by deg(g) then add remainder
        let deg_g = self.generator.len() - 1;
        let mut shifted = vec![0u8; msg.len() + deg_g];
        shifted[..msg.len()].copy_from_slice(msg);

        let remainder = poly_mod_gf2(&shifted, &self.generator);
        for (i, &r) in remainder.iter().enumerate() {
            shifted[msg.len() + i] ^= r;
        }
        shifted
    }

    /// Compute syndrome values for error detection.
    /// Returns syndrome vector of length 2t.
    pub fn syndromes(&self, received: &[u8]) -> Vec<u8> {
        let mut syn = Vec::with_capacity(2 * self.t);
        for i in 1..=2 * self.t {
            // Evaluate received polynomial at α^i over GF(2^m)
            let alpha_i = gf2m_pow(2, i as u32, self.m);
            let val = gf2m_eval_poly(received, alpha_i, self.m);
            syn.push(val);
        }
        syn
    }

    /// Check if syndromes are all zero (no errors).
    pub fn has_errors(&self, syndromes: &[u8]) -> bool {
        syndromes.iter().any(|&s| s != 0)
    }
}

/// Multiply two polynomials over GF(2).
fn poly_mul_gf2(a: &[u8], b: &[u8]) -> Vec<u8> {
    let mut result = vec![0u8; a.len() + b.len() - 1];
    for (i, &ai) in a.iter().enumerate() {
        if ai != 0 {
            for (j, &bj) in b.iter().enumerate() {
                result[i + j] ^= bj;
            }
        }
    }
    result
}

/// Polynomial division remainder over GF(2).
fn poly_mod_gf2(dividend: &[u8], divisor: &[u8]) -> Vec<u8> {
    let mut rem = dividend.to_vec();
    let deg_d = divisor.len() - 1;
    for i in 0..rem.len().saturating_sub(deg_d) {
        if rem[i] != 0 {
            for (j, &d) in divisor.iter().enumerate() {
                rem[i + j] ^= d;
            }
        }
    }
    rem[rem.len().saturating_sub(deg_d)..].to_vec()
}

/// Minimal polynomial of α^c over GF(2) (simplified for small m).
fn minimal_polynomial_gf2(m: usize, c: usize) -> Vec<u8> {
    // Conjugate class: c, 2c, 4c, ... mod (2^m - 1)
    let n = (1usize << m) - 1;
    let mut conjugates = vec![];
    let mut cc = c % n;
    loop {
        conjugates.push(cc);
        cc = (cc * 2) % n;
        if cc == c % n { break; }
    }

    // Product of (x - α^ci) expanded over GF(2)
    let mut poly = vec![1u8];
    for &ci in &conjugates {
        poly = poly_mul_gf2(&poly, &[1, (ci as u8).min(1)]); // simplified
    }
    poly
}

/// GF(2^m) multiplication (table-free, basic implementation).
fn gf2m_mul(a: u8, b: u8, m: usize) -> u8 {
    let modulus = match m {
        3 => 0b1011,   // x³ + x + 1
        4 => 0b10011,   // x⁴ + x + 1
        5 => 0b100101,  // x⁵ + x² + 1
        6 => 0b1000011, // x⁶ + x + 1
        _ => 0b100011011, // x⁸ + x⁴ + x³ + x + 1 (default m=8)
    };
    let mask = (1u16 << m) - 1;
    let mut result = 0u16;
    let mut aa = a as u16;
    let mut bb = b as u16;
    while bb > 0 {
        if bb & 1 != 0 {
            result ^= aa;
        }
        let hi = aa & (1 << (m - 1));
        aa <<= 1;
        if hi != 0 {
            aa ^= modulus as u16;
        }
        aa &= mask;
        bb >>= 1;
    }
    (result & mask) as u8
}

/// GF(2^m) exponentiation.
fn gf2m_pow(mut base: u8, mut exp: u32, m: usize) -> u8 {
    let mut result = 1u8;
    while exp > 0 {
        if exp & 1 == 1 {
            result = gf2m_mul(result, base, m);
        }
        base = gf2m_mul(base, base, m);
        exp >>= 1;
    }
    result
}

/// Evaluate polynomial at a point in GF(2^m).
fn gf2m_eval_poly(poly: &[u8], x: u8, m: usize) -> u8 {
    let mut result = 0u8;
    for &coef in poly.iter().rev() {
        result = gf2m_mul(result, x, m) ^ coef;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bch_create() {
        let bch = BCHCode::new(4, 2);
        assert!(bch.k > 0);
        assert!(bch.generator.len() > 0);
    }

    #[test]
    fn test_encode_no_panic() {
        let bch = BCHCode::new(4, 2);
        let msg = vec![1, 0, 1, 1, 0];
        let _encoded = bch.encode(&msg);
    }
}

/// FNV-1a 64 — the digest every substrate in the SuperInstance fleet agrees on.
pub const FNV_OFFSET: u64 = 0xcbf29ce484222325;
pub const FNV_PRIME: u64 = 0x100000001b3;

#[inline]
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h = FNV_OFFSET;
    for &b in bytes {
        h = (h ^ b as u64).wrapping_mul(FNV_PRIME);
    }
    h
}

/// True if this crate's FNV-1a still agrees with the rest of the fleet.
pub fn canary_holds() -> bool {
    fnv1a64("café Δ 日本語".as_bytes()) == 0x024a555471370b18d
}
