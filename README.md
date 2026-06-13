# bch-code

**BCH (Bose–Chaudhuri–Hocquenghem) cyclic error-correcting codes over GF(2^m) — syndrome computation, generator polynomial construction, and polynomial arithmetic in binary extension fields.**

BCH codes are a family of **cyclic block codes** that achieve a guaranteed error-correcting capability of t errors with a redundancy of at most mt bits, where m is the extension degree of the Galois field. They are constructed by finding a generator polynomial whose roots include 2t consecutive powers of a primitive element α in GF(2^m). BCH codes generalize Hamming codes (t=1) and are the theoretical foundation for Reed-Solomon codes (which operate in GF(p^m) for p > 2).

## Why It Matters

BCH codes are deployed in:

- **Storage media** — CDs, DVDs, and SSDs use BCH codes (and their Reed-Solomon generalization) for error correction at the physical layer.
- **QR codes** — Reed-Solomon codes (a non-binary BCH variant) provide the error correction that lets QR codes survive partial damage.
- **Satellite communication** — DVB-S2 uses BCH codes as an inner code, concatenated with LDPC for near-Shannon-limit performance.
- **Flash memory** — Modern SSD controllers use BCH codes to correct 40+ bits per 1KB page.

The BCH bound guarantees: a BCH code with 2t consecutive roots has **minimum Hamming distance** d_min ≥ 2t + 1, and thus can correct any combination of t errors. This is the **Singleton bound** for cyclic codes — you cannot do better with the same redundancy.

## How It Works

### Finite Field GF(2^m)

All arithmetic is performed in GF(2^m), the binary extension field with 2^m elements. Elements are m-bit vectors, and operations use modular polynomial arithmetic:

- **Addition**: XOR (component-wise mod 2)
- **Multiplication**: Polynomial multiplication mod an irreducible polynomial p(x)

Irreducible polynomials used:

| m | p(x) | Primitive polynomial |
|---|------|---------------------|
| 3 | x³ + x + 1 | 0b1011 |
| 4 | x⁴ + x + 1 | 0b10011 |
| 5 | x⁵ + x² + 1 | 0b100101 |
| 6 | x⁶ + x + 1 | 0b1000011 |
| 8 | x⁸ + x⁴ + x³ + x + 1 | 0b100011011 |

Multiplication uses shift-and-XOR (analogous to schoolbook multiplication with modular reduction):

```
result = 0
while b > 0:
    if b & 1: result ^= a
    a <<= 1
    if a & high_bit: a ^= modulus
    b >>= 1
```

### Generator Polynomial Construction

For a t-error-correcting BCH code of length n = 2^m - 1:

1. Find the minimal polynomials M₁(x), M₃(x), ..., M_{2t-1}(x) of α, α³, ..., α^{2t-1}
2. The generator polynomial: **g(x) = LCM(M₁, M₃, ..., M_{2t-1})**
3. The code dimension: k = n - deg(g)

Each minimal polynomial M_i(x) is the product of (x - α^{i·2^j}) over the cyclotomic coset of i:

> Coset of c: {c, 2c, 4c, ...} mod (2^m - 1)

For example, in GF(2⁴) with t=2:
- Coset of 1: {1, 2, 4, 8}
- Coset of 3: {3, 6, 12, 9}
- g(x) = M₁(x) · M₃(x), degree ≤ 8, so k ≥ 15 - 8 = 7

### Encoding

Systematic encoding via polynomial division:

> codeword = message · x^{deg(g)} + (message · x^{deg(g)} mod g(x))

The message occupies the high-order coefficients; the parity check symbols are the remainder.

### Syndrome Computation

For a received word r(x), the syndromes are:

> S_i = r(α^i) for i = 1, 2, ..., 2t

If all syndromes are zero, no errors are detected. Non-zero syndromes trigger error location via the Berlekamp-Massey algorithm (not implemented in this crate — this crate provides syndrome computation as the first step).

### Polynomial Arithmetic Over GF(2)

All polynomial operations use **GF(2) coefficients** (i.e., addition is XOR):

- `poly_mul_gf2(a, b)`: Schoolbook multiplication with XOR accumulation
- `poly_mod_gf2(dividend, divisor)`: Long division with XOR subtraction

### Complexity

| Operation | Time | Notes |
|-----------|------|-------|
| Field multiplication (GF(2^m)) | O(m²) | Shift-and-XOR |
| Field exponentiation | O(m³ log e) | Square-and-multiply |
| Generator polynomial construction | O(t · m²) | t minimal polynomials |
| Encoding (k-bit message) | O(k · deg(g)) | Polynomial division |
| Syndrome computation | O(n · t · m²) | 2t evaluations |

## Quick Start

```rust
use bch_code::BCHCode;

// Create a BCH(15, 7, 2) code: 15-bit codewords, 7-bit messages, 2 errors
let bch = BCHCode::new(4, 2);
println!("BCH({}, {}, {})", bch.n, bch.k, bch.t);

// Encode a 7-bit message
let msg = vec![1, 0, 1, 1, 0, 1, 0];
let codeword = bch.encode(&msg);
println!("Encoded: {:?}", codeword);

// Compute syndromes (all zero = no errors)
let syndromes = bch.syndromes(&codeword);
println!("Has errors: {}", bch.has_errors(&syndromes));
```

## API

- **`BCHCode`** — { n, k, t, m, generator }: `new(m, t)`, `encode(msg)`, `syndromes(received)`, `has_errors(syndromes)`
- **`gf2m_mul(a, b, m) → u8`** — GF(2^m) multiplication
- **`gf2m_pow(base, exp, m) → u8`** — GF(2^m) exponentiation
- **`gf2m_eval_poly(poly, x, m) → u8`** — Evaluate polynomial at GF(2^m) point
- **`poly_mul_gf2(a, b) → Vec<u8>`** — GF(2) polynomial multiplication
- **`poly_mod_gf2(dividend, divisor) → Vec<u8>`** — GF(2) polynomial remainder

## Architecture Notes

The γ+η=C identity: γ (generative capacity) is the number of valid codewords = 2^k — the information space. η (evaluative depth) is the error-correcting capability t — the code's resilience to noise. C = the code rate k/n × the error correction t, representing the bandwidth-reliability tradeoff. BCH codes achieve near-optimal C by maximizing both the information rate (high k/n) and error correction (high t) simultaneously. The BCH bound proves this is the best possible for cyclic codes.

## References

1. Bose, R. & Ray-Chaudhuri, D. (1960). "On a Class of Error Correcting Binary Group Codes." *Information and Control*, 3(1). — Original BCH paper.
2. Hocquenghem, A. (1959). "Codes correcteurs d'erreurs." *Chiffres*, 2. — Independent discovery.
3. Lin, S. & Costello, D. (2004). *Error Control Coding* (2nd ed.), Ch. 6. — Standard BCH code textbook.
4. MacWilliams, F. & Sloane, N. (1977). *The Theory of Error-Correcting Codes*. North-Holland. — Cyclotomic cosets and minimal polynomials.

## License

MIT
