#!/usr/bin/env python3

import numpy as np
import time, sys, os, re, warnings

warnings.filterwarnings("ignore")

NQUBITS = 8
ROUNDS = 12
KEY1_INDEX = [3, 1, 6, 0, 5, 4, 7, 2]
PERMUTATION = [5, 2, 7, 0, 3, 1, 4, 6]
RXZ_INTERACTION = [(0, 3), (4, 5), (7, 2), (1, 6)]
THETA_CRY = np.pi / 15
THETA_RZX = np.pi / np.e

def byte_to_state(bv):
    s = np.zeros(256, dtype=np.complex128)
    s[bv] = 1.0
    return s

def apply_ry(state, qubit, theta):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    gate = np.array([[c, -s], [s, c]], dtype=np.complex128)
    shape = tuple([2] * NQUBITS)
    r = state.reshape(shape)
    r = np.moveaxis(r, qubit, -1).reshape(128, 2)
    r = r @ gate.T
    r = r.reshape(shape)
    r = np.moveaxis(r, -1, qubit)
    return r.reshape(-1)

def apply_cry(state, key_byte):
    for i in range(NQUBITS):
        if (key_byte >> KEY1_INDEX[i]) & 1:
            state = apply_ry(state, i, THETA_CRY)
    return state

def one_round(state, kb, V):
    return V @ apply_cry(state, kb)

def encrypt_full(byte_val, key_bytes, V):
    state = byte_to_state(byte_val)
    for r in range(ROUNDS):
        state = one_round(state, key_bytes[r], V)
    return state

def build_V():
    """Build the fixed unitary V = ISWAP_perm . RZX_interaction"""
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Operator
    qc = QuantumCircuit(NQUBITS)
    for i in range(NQUBITS):
        qc.iswap(PERMUTATION[i], PERMUTATION[(i + 1) % NQUBITS])
    for q1, q2 in RXZ_INTERACTION:
        qc.rzx(THETA_RZX, q1, q2)
    return np.array(Operator(qc).data, dtype=np.complex128)

def fetch_from_remote(host, port):
    """Connect to the challenge server and grab the encrypted flag."""
    import socket
    print(f"[*] Connecting to {host}:{port} ...")
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.settimeout(600)
    s.connect((host, port))
    data = b''
    while b'> ' not in data:
        chunk = s.recv(1048576)
        if not chunk:
            break
        data += chunk

    menu = data.decode('utf-8', errors='replace')
    print("[*] Connected. Server says:")
    for line in menu.strip().split('\n')[:6]:
        print(f"    {line}")

    print("[*] Sending option 2 (get encrypted flag)...")
    s.send(b'2\n')

    print("[*] Waiting for encrypted flag (this takes a while)...")
    data2 = b''
    while b'> ' not in data2:
        chunk = s.recv(1048576)
        if not chunk:
            break
        data2 += chunk
        dots = '.' * (min(len(data2) // 50000, 50))
        sys.stdout.write(f"\r    Received {len(data2):,} bytes {dots}")
        sys.stdout.flush()
    print()

    s.close()

    text = data2.decode('utf-8', errors='replace')
    numbers = re.findall(r'np\.complex128\(([^)]+)\)', text)
    ct = np.array([complex(n.replace(' ', '')) for n in numbers], dtype=np.complex128)
    print(f"[*] Parsed {len(ct)} complex numbers  ->  {len(ct)//256} flag bytes\n")
    return ct

def fetch_from_npy(path):
    """Load pre-saved encrypted flag from .npy file."""
    ct = np.load(path).astype(np.complex128)
    print(f"[*] Loaded {path}: {len(ct)} complex numbers  ->  {len(ct)//256} flag bytes\n")
    return ct

def fetch_from_chall(filepath):
    """Import QCipher from chall.py and encrypt flag in-memory (offline only)."""
    import importlib.util
    print(f"[*] Importing QCipher from {filepath} ...")
    spec = importlib.util.spec_from_file_location("chall", filepath)
    chall = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(chall)
    qcipher = chall.QCipher(chall.NQUBITS, chall.ROUNDS)
    flag = chall.FLAG
    print(f"[*] Flag in chall.py: {flag}")
    print(f"[*] NQUBITS={chall.NQUBITS}, ROUNDS={chall.ROUNDS}")
    print("[*] Encrypting flag ...")
    ct_list = qcipher.encrypt(flag)
    ct = np.array(ct_list, dtype=np.complex128)
    print(f"[*] Got {len(ct)} complex numbers  ->  {len(ct)//256} flag bytes\n")
    return ct, flag.decode() if isinstance(flag, bytes) else flag

def solve(flag_ct, V):
    """Recover the 12-byte key using coordinate descent."""
    flag_len = len(flag_ct) // 256

    known = [(0, ord('T')), (1, ord('R')), (2, ord('X')),
             (3, ord('{')), (flag_len - 1, ord('}'))]
    targets = {bv: flag_ct[p * 256:(p + 1) * 256] for p, bv in known}
    print(f"[*] Known plaintexts: {[(chr(b), p) for p, b in known]}")

    print("[*] Starting coordinate descent key recovery ...")
    key = [0] * ROUNDS
    t0 = time.time()

    for iteration in range(10):
        improved = False
        for r in range(ROUNDS):
            pre = {}
            for bv in targets:
                st = byte_to_state(bv)
                for rr in range(r):
                    st = one_round(st, key[rr], V)
                pre[bv] = st

            best_kr, best_loss = key[r], float('inf')
            for kr in range(256):
                loss = 0.0
                for bv, target in targets.items():
                    st = one_round(pre[bv].copy(), kr, V)
                    for rr in range(r + 1, ROUNDS):
                        st = one_round(st, key[rr], V)
                    loss += np.sum(np.abs(st - target) ** 2)
                if loss < best_loss:
                    best_loss, best_kr = loss, kr

            if best_kr != key[r]:
                improved = True
            key[r] = best_kr

            bar = '#' * int(min((best_loss / 2.0) * 30, 30))
            sys.stdout.write(
                f"\r    iter {iteration}  round {r:2d}/11  "
                f"key[{r:2d}]={best_kr:3d}  loss={best_loss:.4e}  {bar}"
            )
            sys.stdout.flush()

        total = sum(np.sum(np.abs(encrypt_full(bv, key, V) - t) ** 2)
                    for bv, t in targets.items())
        print(f"\n    -- total loss = {total:.4e}  ({time.time() - t0:.1f}s)")
        if not improved:
            print("    -- CONVERGED\n")
            break

    print(f"[*] Recovered key: {key}")
    return key

def decrypt_flag(flag_ct, key, V):
    """Decrypt the flag by trying all 256 possible bytes for each position."""
    flag_len = len(flag_ct) // 256
    flag = ''
    for pos in range(flag_len):
        block = flag_ct[pos * 256:(pos + 1) * 256]
        found = False
        for bv in range(256):
            if np.allclose(encrypt_full(bv, key, V), block, atol=1e-8):
                flag += chr(bv)
                found = True
                break
        if not found:
            best_bv, best_dist = 0, float('inf')
            for bv in range(32, 127):
                dist = np.max(np.abs(encrypt_full(bv, key, V) - block))
                if dist < best_dist:
                    best_dist, best_bv = dist, bv
            flag += chr(best_bv)
            print(f"    WARNING pos {pos}: closest '{chr(best_bv)}' dist={best_dist:.2e}")
    return flag

if __name__ == '__main__':
    t_start = time.time()
    print("=" * 60)
    print("  Quantum Cipher CTF Solver")
    print("=" * 60 + "\n")

    print("[1] Building fixed unitary V (ISWAP + RZX) ...")
    V = build_V()
    err = np.max(np.abs(V @ V.conj().T - np.eye(256)))
    print(f"    shape = {V.shape}   |VV-I| = {err:.2e}\n")

    print("[2] Getting encrypted flag ...")
    flag_ct = None
    true_flag = None

    if os.path.exists('encrypted_flag.npy'):
        try:
            flag_ct = fetch_from_npy('encrypted_flag.npy')
        except Exception as e:
            print(f"    Failed to load .npy: {e}\n")

    if flag_ct is None:
        try:
            flag_ct = fetch_from_remote('quantumcipher.ctf.theromanxpl0.it', 9099)
        except Exception as e:
            print(f"    Remote connection failed: {e}\n")

    if flag_ct is None and os.path.exists('chall.py'):
        try:
            flag_ct, true_flag = fetch_from_chall('chall.py')
        except Exception as e:
            print(f"    Local import failed: {e}\n")

    if flag_ct is None:
        print("ERROR: Could not get encrypted flag.")
        print("  - Put encrypted_flag.npy in the same directory, OR")
        print("  - Make sure the remote server is up, OR")
        print("  - Put chall.py in the same directory (for offline testing)")
        sys.exit(1)

    print("[3] Recovering encryption key ...")
    key = solve(flag_ct, V)

    print("[4] Verifying key against known plaintexts ...")
    flag_len = len(flag_ct) // 256
    known = [(0, ord('T')), (1, ord('R')), (2, ord('X')),
             (3, ord('{')), (flag_len - 1, ord('}'))]
    targets = {bv: flag_ct[p * 256:(p + 1) * 256] for p, bv in known}
    for pos, bv in known:
        comp = encrypt_full(bv, key, V)
        e = np.max(np.abs(comp - targets[bv]))
        ok = 'OK' if e < 1e-6 else 'FAIL'
        print(f"    {ok}  '{chr(bv)}' (pos {pos:2d})  err = {e:.2e}")
    print()

    print("[5] Decrypting flag ...")
    flag = decrypt_flag(flag_ct, key, V)

    print(f"\n{'=' * 60}")
    print(f"  FLAG: {flag}")
    print(f"{'=' * 60}")
    if true_flag:
        match = "MATCH" if flag == true_flag else "MISMATCH"
        print(f"  Verification: {match}")
    print(f"  Time: {time.time() - t_start:.1f}s")
    print(f"{'=' * 60}")
