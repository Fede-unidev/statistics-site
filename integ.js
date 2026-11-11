// Riemann- and Lebesgue-style numerical integration (HW6, Ex. 1)
const INT = (() => {
  const F = {
    sin:     { name: 'sin(x), [0, π]', f: Math.sin, a: 0, b: Math.PI, exact: 2, p: 1 },
    sq:      { name: 'x², [0, 1]', f: x => x * x, a: 0, b: 1, exact: 1 / 3, p: 1 },
    step3:   { name: '1 if x < 0.5, 3 if x ≥ 0.5, [0, 1]', f: x => x < 0.5 ? 1 : 3, a: 0, b: 1, exact: 2, p: 1 },
    abs:     { name: '|x − 0.5|, [0, 1]', f: x => Math.abs(x - 0.5), a: 0, b: 1, exact: 0.25, p: 1 },
    invsqrt: { name: '1/√x, (0, 1]  (unbounded)', f: x => 1 / Math.sqrt(x), a: 0, b: 1, exact: 2, p: 3 },
    step1:   { name: '1 if x < 0.5, 0 if x ≥ 0.5, [0, 1]', f: x => x < 0.5 ? 1 : 0, a: 0, b: 1, exact: 0.5, p: 1 },
    sinc:    { name: 'sin(x)/x, [0, π]', f: x => x === 0 ? 1 : Math.sin(x) / x, a: 0, b: Math.PI, exact: 1.8519370519824662, p: 1 },
  };
  // exact Lebesgue measure of the sublevel set G(t) = mu{x in [a,b] : f(x) < t}, and the range [ymin, ymax] of f
  const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
  F.sin.G = t => 2 * Math.asin(clamp(t, 0, 1)); F.sin.ymin = 0; F.sin.ymax = 1;
  F.sq.G = t => Math.sqrt(clamp(t, 0, 1)); F.sq.ymin = 0; F.sq.ymax = 1;
  F.step3.G = t => t <= 1 ? 0 : t <= 3 ? 0.5 : 1; F.step3.ymin = 1; F.step3.ymax = 3;
  F.abs.G = t => 2 * clamp(t, 0, 0.5); F.abs.ymin = 0; F.abs.ymax = 0.5;
  F.step1.G = t => t <= 0 ? 0 : t <= 1 ? 0.5 : 1; F.step1.ymin = 0; F.step1.ymax = 1;
  F.invsqrt.G = t => t <= 1 ? 0 : 1 - 1 / (t * t); F.invsqrt.ymin = 1; F.invsqrt.ymaxOf = K => K * K; // truncated at K^2
  F.sinc.G = t => { if (t <= 0) return 0; if (t >= 1) return Math.PI; let lo = 0, hi = Math.PI; for (let i = 0; i < 60; i++) { const m = (lo + hi) / 2; (Math.sin(m) / m > t) ? lo = m : hi = m; } return Math.PI - lo; }; F.sinc.ymin = 0; F.sinc.ymax = 1;
  // Lebesgue sum with EXACT preimage measures: slab j = {y_j <= f < y_{j+1}} has measure G(y_{j+1}) - G(y_j); cost O(K), no function samples
  function lebesgueExact(fn, K) {
    const { a, b, p, G, ymin } = fn, ymax = fn.ymaxOf ? fn.ymaxOf(K) : fn.ymax; let lw = 0, mi = 0, up = 0, prev = G(ymin);
    for (let j = 0; j < K; j++) {
      const y0 = ymin + (ymax - ymin) * Math.pow(j / K, p), y1 = ymin + (ymax - ymin) * Math.pow((j + 1) / K, p), g1 = j === K - 1 ? (b - a) : G(y1), m = g1 - prev; prev = g1;
      lw += y0 * m; up += y1 * m; mi += 0.5 * (y0 + y1) * m;
    }
    return { lower: lw, mid: mi, upper: up };
  }
  // Riemann sums on a uniform partition of [a,b] into N subintervals
  function riemann(fn, N) {
    const { f, a, b } = fn, h = (b - a) / N; let L = 0, M = 0, R = f(b) ;
    for (let i = 0; i < N; i++) { L += f(a + i * h); M += f(a + (i + 0.5) * h); if (i > 0) R += f(a + i * h); }
    return { left: L * h, right: R * h, mid: M * h, trap: (L + R) * h / 2 };
  }
  // Lebesgue-style sum: partition the RANGE into K levels y_0 < ... < y_K (spacing power p: y_j = ymin + (ymax-ymin)(j/K)^p),
  // estimate mu(f^{-1}[y_j, y_{j+1})) by the fraction of the N midpoint samples whose value falls in that slab.
  function lebesgue(fn, N, K) {
    const { f, a, b, p } = fn, h = (b - a) / N, v = new Float64Array(N); let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < N; i++) { const y = f(a + (i + 0.5) * h); v[i] = y; if (y < lo) lo = y; if (y > hi) hi = y; }
    if (hi === lo) return { lower: lo * (b - a), mid: lo * (b - a), upper: lo * (b - a) };
    const cnt = new Float64Array(K), w = hi - lo, ip = 1 / p;
    for (let i = 0; i < N; i++) { let j = Math.floor(K * Math.pow((v[i] - lo) / w, ip)); if (j >= K) j = K - 1; cnt[j]++; }
    let lw = 0, mi = 0, up = 0;
    for (let j = 0; j < K; j++) { const y0 = lo + w * Math.pow(j / K, p), y1 = lo + w * Math.pow((j + 1) / K, p), m = cnt[j] / N * (b - a); lw += y0 * m; up += y1 * m; mi += 0.5 * (y0 + y1) * m; }
    return { lower: lw, mid: mi, upper: up };
  }
  return { F, riemann, lebesgue, lebesgueExact };
})();
if (typeof module !== 'undefined') module.exports = INT;
