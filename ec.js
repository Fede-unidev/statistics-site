// Curve y² = x³ + b (mod p), small odd prime p (p < ~1e7 to stay within Number precision)
function makeCurve(p, b) {
  const mod = n => ((n % p) + p) % p;
  const pw = (x, e) => { let r = 1; x = mod(x); while (e > 0) { if (e & 1) r = r * x % p; x = x * x % p; e = Math.floor(e / 2); } return r; };
  const inv = n => pw(n, p - 2);                       // Fermat: n^(p-2) = n^-1
  const rhs = x => mod(mod(x * x) * x + b);
  const on = P => mod(P[1] * P[1]) === rhs(P[0]);
  const has = x => { const r = rhs(x); return r === 0 ? 1 : (pw(r, (p - 1) / 2) === 1 ? 2 : 0); }; // number of points with abscissa x
  function add(P, Q) {                                  // null = point at infinity O
    if (!P) return Q; if (!Q) return P;
    if (P[0] === Q[0] && mod(P[1] + Q[1]) === 0) return null;
    const l = P[0] === Q[0] ? mod(3 * P[0] * P[0]) * inv(2 * P[1]) % p
                            : mod(Q[1] - P[1]) * inv(Q[0] - P[0]) % p;
    const x = mod(l * l - P[0] - Q[0]);
    return [x, mod(l * (P[0] - x) - P[1])];
  }
  function mul(k, G) { let R = null, A = G; while (k > 0) { if (k & 1) R = add(R, A); A = add(A, A); k = Math.floor(k / 2); } return R; }
  function order(G) { let R = G, n = 1; while (R) { R = add(R, G); n++; } return n; }
  function find() { for (let x = 1; x < p; x++) if (has(x)) for (let y = 1; y < p; y++) if (mod(y * y) === rhs(x)) return [x, y]; }
  function dlog(Q, G) { let R = null, j = 0; do { R = add(R, G); j++; } while (!(R && Q && R[0] === Q[0] && R[1] === Q[1]) && R); return j; } // brute force
  return { p, b, add, mul, order, find, on, has, dlog };
}
if (typeof module !== 'undefined') module.exports = { makeCurve };
