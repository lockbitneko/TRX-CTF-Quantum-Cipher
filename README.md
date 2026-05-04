# TRX CTF 2025 | Quantum Cipher Solver

> Coordinate descent key recovery against a "quantum" cipher that's entirely classical.

## Challenge

- **Name:** Quantum Cipher
- **Category:** Crypto
- **Author:** theromanxpl0it
- **Server:** `nc quantumcipher.ctf.theromanxpl0.it 9099`
- **Flag:** `REDACTED`

## How It Works

The cipher encrypts each flag byte independently (ECB mode) using 8-qubit quantum circuits through 12 rounds. Each round applies:

1. **CRY layer** - key-dependent RY rotations on selected qubits
2. **ISWAP permutation** - fixed qubit swap pattern
3. **RZX interaction** - fixed two-qubit rotation gates

The server returns the **full 256-element statevector** per byte, not a quantum measurement, but the complete complex amplitude vector. This leaks enough information to classically simulate and invert the entire encryption.

## Attack

Since ISWAP + RZX layers are key-independent, they're precomputed into a single 256×256 unitary matrix **V**. Each round becomes `state = V · CRY(state, key_byte)`.

With 5 known plaintexts from the flag format `TRX{...}`, we recover the 12-byte key using **coordinate descent**:

- For each round, try all 256 key byte values
- Pick the one minimizing `Σ ||enc(b) - target||²` over known pairs
- Repeat until convergence (~2 iterations)

Result: 256¹² ≈ 10²⁹ brute-force → ~6,000 evaluations. Solved in ~8 seconds.

## Files

| File | Description |
|------|-------------|
| `solve.py` | Solver script - builds V, recovers key via coordinate descent, decrypts flag |
| `encrypted_flag.npy` | Encrypted flag captured from the remote server (9728 complex numbers = 38 flag bytes × 256 amplitudes) |

## Usage

```bash
pip3 install numpy qiskit pycryptodome --break-system-packages
python3 solve.py
```

The solver also supports connecting directly to the remote server if `encrypted_flag.npy` is not present.

## Output

```bash
[1] Building fixed unitary V (ISWAP + RZX) ...
    shape = (256, 256)   |VV-I| = 4.44e-16

[2] Getting encrypted flag ...
[*] Loaded encrypted_flag.npy: 9728 complex numbers  ->  38 flag bytes

[3] Recovering encryption key ...
[*] Known plaintexts: [('T', 0), ('R', 1), ('X', 2), ('{', 3), ('}', 37)]
[*] Starting coordinate descent key recovery ...
    -- total loss = 3.4882e-30  (7.9s)
    -- CONVERGED
[*] Recovered key: [230, 136, 229, 30, 16, 211, 2, 190, 158, 211, 151, 214]

[4] Verifying key against known plaintexts ...
    OK  'T' (pos  0)  err = 1.26e-16
    OK  'R' (pos  1)  err = 1.19e-16
    OK  'X' (pos  2)  err = 1.16e-16
    OK  '{' (pos  3)  err = 1.52e-16
    OK  '}' (pos 37)  err = 1.21e-16

[5] Decrypting flag ...

============================================================
  FLAG: TRX{***************************}
============================================================
```

## Dependencies

- Python 3.8+
- numpy
- qiskit
- pycryptodome
