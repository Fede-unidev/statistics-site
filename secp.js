// Toy ECIES-style scheme on the REAL secp256k1 curve (BigInt). The private key k is restricted to a small range.
// Needs a global keccak256 (js-sha3).
const SECP = (() => {
  const P = 2n ** 256n - 2n ** 32n - 977n;
  const N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141n;
  const G = [0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798n,
             0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8n];
  const mod = a => { const r = a % P; return r >= 0n ? r : r + P; };
  const pw = (b, e) => { let r = 1n; b = mod(b); while (e > 0n) { if (e & 1n) r = r * b % P; b = b * b % P; e >>= 1n; } return r; };
  function inv(a) { let r0 = P, r1 = mod(a), s0 = 0n, s1 = 1n; while (r1) { const q = r0 / r1; [r0, r1] = [r1, r0 - q * r1]; [s0, s1] = [s1, s0 - q * s1]; } return mod(s0); }
  // affine arithmetic (null = point at infinity)
  function add(A, B) {
    if (!A) return B; if (!B) return A;
    let l;
    if (A[0] === B[0]) { if (mod(A[1] + B[1]) === 0n) return null; l = mod(3n * A[0] * A[0] * inv(2n * A[1])); }
    else l = mod((B[1] - A[1]) * inv(B[0] - A[0]));
    const x = mod(l * l - A[0] - B[0]);
    return [x, mod(l * (A[0] - x) - A[1])];
  }
  const neg = A => A && [A[0], mod(-A[1])];
  // Jacobian arithmetic (no inversions): fast scalar multiplication and brute-force walking
  const JINF = [0n, 1n, 0n];
  function dbl(A) {
    const [X, Y, Z] = A; if (Z === 0n || Y === 0n) return JINF;
    const a = X * X % P, b = Y * Y % P, c = b * b % P, d = mod(2n * ((X + b) * (X + b) - a - c)), e = 3n * a, f = e * e % P, x3 = mod(f - 2n * d);
    return [x3, mod(e * (d - x3) - 8n * c), 2n * Y * Z % P];
  }
  function addJ(A, B) {
    if (A[2] === 0n) return B; if (B[2] === 0n) return A;
    const [X1, Y1, Z1] = A, [X2, Y2, Z2] = B;
    const z1 = Z1 * Z1 % P, z2 = Z2 * Z2 % P, u1 = X1 * z2 % P, u2 = X2 * z1 % P, s1 = Y1 * Z2 % P * z2 % P, s2 = Y2 * Z1 % P * z1 % P;
    const h = mod(u2 - u1), r0 = mod(s2 - s1);
    if (h === 0n) return r0 === 0n ? dbl(A) : JINF;
    const t = 2n * h % P, i = t * t % P, j = h * i % P, r = 2n * r0 % P, v = u1 * i % P, x3 = mod(r * r - j - 2n * v);
    return [x3, mod(r * (v - x3) - 2n * s1 * j), mod(((Z1 + Z2) * (Z1 + Z2) - z1 - z2) * h)];
  }
  const toJ = A => [A[0], A[1], 1n];
  function toA(A) { if (A[2] === 0n) return null; const zi = inv(A[2]), z2 = zi * zi % P; return [A[0] * z2 % P, A[1] * z2 % P * zi % P]; }
  const eqJA = (J, A) => { if (J[2] === 0n) return false; const z2 = J[2] * J[2] % P; return J[0] === A[0] * z2 % P && J[1] === A[1] * z2 % P * J[2] % P; };
  function mul(k, A) { let R = JINF, Q = toJ(A); while (k > 0n) { if (k & 1n) R = addJ(R, Q); Q = dbl(Q); k >>= 1n; } return toA(R); }
  // encodings
  const hex = (n, w = 64) => n.toString(16).padStart(w, '0');
  const b2h = b => [...b].map(v => v.toString(16).padStart(2, '0')).join('');
  const h2b = h => (h.match(/../g) || []).map(x => parseInt(x, 16));
  const comp = A => (A[1] & 1n ? '03' : '02') + hex(A[0]);
  function decomp(h) {
    const x = BigInt('0x' + h.slice(2)); let y = pw(mod(x * x * x + 7n), (P + 1n) / 4n);
    if ((y & 1n) !== BigInt(parseInt(h.slice(0, 2), 16) - 2)) y = P - y; return [x, y];
  }
  const on = A => mod(A[1] * A[1]) === mod(A[0] * A[0] * A[0] + 7n);
  // keystream = Keccak-256(S.x || counter), block by block
  function stream(S, len) {
    const out = [], sx = h2b(hex(S[0])); let c = 0;
    while (out.length < len) { out.push(...keccak256.array(sx.concat([0, 0, 0, c++]))); }
    return out.slice(0, len);
  }
  const xor = (a, b) => a.map((v, i) => v ^ b[i]);
  function encrypt(bytes, Q, r) { const R = mul(r, G), S = mul(r, Q); return { R, c: xor(bytes, stream(S, bytes.length)) }; }
  function decrypt(R, c, k) { const S = mul(k, R); return xor(c, stream(S, c.length)); }
  const rand = bits => { const a = crypto.getRandomValues(new Uint8Array(Math.ceil(bits / 8))); let v = 0n; a.forEach(x => v = (v << 8n) | BigInt(x)); return (v & ((1n << BigInt(bits)) - 1n)) || 1n; };
  const tick = () => new Promise(r => setTimeout(r));
  // Attack 1: exhaustive search over k in [1, 2^bits): walk kG and compare with Q
  async function brute(Q, bits, prog) {
    const lim = 1n << BigInt(bits), Gj = toJ(G); let cur = Gj, k = 1n;
    while (k < lim) {
      for (let i = 0; i < 1500 && k < lim; i++, k++) { if (eqJA(cur, Q)) return k; cur = addJ(cur, Gj); }
      if (prog) prog(k); await tick();
    }
    return null;
  }
  // Attack 2: baby-step giant-step, O(sqrt(2^bits)) time and memory
  async function bsgs(Q, bits, prog) {
    const m = 1n << BigInt((bits + 1) >> 1), mm = Number(m), tab = new Map(); let B = null;
    for (let j = 1; j <= mm; j++) { B = add(B, G); tab.set(B[0], [j, B[1]]); if (j % 2000 === 0) { if (prog) prog('baby ' + j + '/' + mm); await tick(); } }
    const step = neg(mul(m, G)); let cur = Q;
    for (let i = 0; i <= mm; i++) {
      if (cur === null) return BigInt(i) * m;
      const hit = tab.get(cur[0]);
      if (hit) { const k = hit[1] === cur[1] ? BigInt(i) * m + BigInt(hit[0]) : BigInt(i) * m - BigInt(hit[0]); if (k > 0n && eqA(mul(k, G), Q)) return k; }
      cur = add(cur, step); if (i % 2000 === 0) { if (prog) prog('giant ' + i + '/' + mm); await tick(); }
    }
    return null;
  }
  const eqA = (A, B) => A && B && A[0] === B[0] && A[1] === B[1];
  return { P, N, G, add, mul, comp, decomp, on, hex, b2h, h2b, stream, encrypt, decrypt, rand, brute, bsgs, eqA };
})();
if (typeof module !== 'undefined') module.exports = SECP;
