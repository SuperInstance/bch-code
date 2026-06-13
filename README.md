# BCH Code

**A Rust library for BCH (Bose–Chaudhuri–Hocquenghem) error-correcting codes** — cyclic block codes over GF(2) that can detect and correct multiple bit errors per codeword.

## Why It Matters

BCH codes are the workhorses of digital communication and storage. They're used in:

- **NAND flash memory** — BCH is the standard ECC for SSD controllers
- **QR codes** — BCH codes encode format information
- **Satellite communication** — forward error correction in noisy channels
- **Optical discs** — DVD and Blu-ray use BCH-based codes

The key advantage of BCH codes: you can **design them for any error-correcting capability *t***. For a code over GF(2^m) with block length n = 2^m − 1, choosing parameter *t* gives you a code that corrects up to *t* bit errors per codeword. The Reed-Solomon code is a special (non-binary) case of BCH.

The encoding overhead is *n − k* parity bits (where *k* is the message length), and both encoding and decoding are polynomial-time: encoding is O(n·k), syndrome computation is O(n·t), and error location via the Berlekamp-Massey algorithm is O(t²).

## How It Works

**Code construction**: Given extension degree *m* and error capacity *t*, the code length is n = 2^m − 1. The generator polynomial g(x) is the LCM of the minimal polynomials of α¹, α³, α⁵, ..., α^(2t−1) over GF(2), where α is a primitive element of GF(2^m). The degree of g(x) determines the number of parity bits (n − k).

**GF(2^m) arithmetic**: The library implements multiplication in the Galois field using the standard shift-and-XOR algorithm with irreducible polynomials for m = 3, 4, 5, 6, 8. For example, GF(2⁴) uses the modulus x⁴ + x + 1 = 0b10011. Exponentiation uses repeated squaring.

**Encoding**: Systematic encoding via polynomial division — shift the message up by deg(g), compute the remainder mod g(x), and append it as parity bits.

**Syndrome computation**: For error detection, evaluate the received polynomial at α^i for i = 1..2t. All-zero syndromes mean no errors.

## Quick Start

```rust
use bch_code::BCHCode;

// Create a BCH(15, 7, 2) code: 15-bit codewords, 7-bit messages, 2-error-correcting
let bch = BCHCode::new(4, 2);
println!("Code: ({}, {}, t={})", bch.n, bch.k, bch.t);

// Encode a message
let msg = vec![1, 0, 1, 1, 0, 0, 1];
let codeword = bch.encode(&msg);

// Check for errors
let syndromes = bch.syndromes(&codeword);
if !bch.has_errors(&syndromes) {
    println!("No errors detected");
}
```

## API

- **`BCHCode`** — BCH code with parameters n, k, t, m, and generator polynomial
  - `new(m, t)` — Construct code over GF(2^m) correcting t errors
  - `encode(msg)` — Systematic encoding via polynomial division
  - `syndromes(received)` — Compute 2t syndromes for error detection
  - `has_errors(syndromes)` — Check if any syndrome is non-zero

## Architecture Notes

Provides the forward error correction layer for SuperInstance data reliability pipelines. The GF(2^m) arithmetic primitives are reusable for Reed-Solomon and other algebraic codes. See the [architecture overview](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## License

MIT
