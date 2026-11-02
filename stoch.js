// Random-walk simulation helpers for HW4
const STO = (() => {
  const params = (model, p) => model === 'rad' ? { mu: 0, s2: 1, step: 2 } : { mu: p, s2: p * (1 - p), step: 1 };
  // Phi(z): standard normal CDF (Abramowitz-Stegun 7.1.26 erf, |error| < 1.5e-7)
  function Phi(z) {
    const x = Math.abs(z) / Math.SQRT2, t = 1 / (1 + 0.3275911 * x);
    const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
    return z >= 0 ? 0.5 * (1 + y) : 0.5 * (1 - y);
  }
  // m independent attackers, each probes N = max(ns) systems; record S_n at every n in ns (sorted, ascending)
  function simulate(model, p, m, ns) {
    const N = ns[ns.length - 1], out = ns.map(() => new Float64Array(m)), rad = model === 'rad';
    for (let j = 0; j < m; j++) {
      let S = 0, idx = 0, next = ns[0];
      for (let i = 1; i <= N; i++) {
        S += rad ? (Math.random() < 0.5 ? 1 : -1) : (Math.random() < p ? 1 : 0);
        if (i === next) { out[idx][j] = S; idx++; next = ns[idx]; }
      }
    }
    return out;
  }
  // one trajectory, S_k recorded at checkpoints pts (sorted)
  function path(model, p, pts) {
    const r = [], rad = model === 'rad'; let S = 0, idx = 0;
    for (let i = 1; i <= pts[pts.length - 1]; i++) {
      S += rad ? (Math.random() < 0.5 ? 1 : -1) : (Math.random() < p ? 1 : 0);
      if (i === pts[idx]) { r.push(S); idx++; }
    }
    return r;
  }
  const checkpoints = N => { const s = new Set(); for (let j = 0; j < 300; j++) s.add(Math.max(1, Math.round(Math.pow(N, j / 299)))); return [...s].sort((a, b) => a - b); };
  const mean = a => a.reduce((s, v) => s + v, 0) / a.length;
  const variance = a => { const m = mean(a); return a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1); };
  function ks(z) { const s = Float64Array.from(z).sort(), n = s.length; let d = 0; for (let i = 0; i < n; i++) { const f = Phi(s[i]); d = Math.max(d, (i + 1) / n - f, f - i / n); } return d; }
  return { params, Phi, simulate, path, checkpoints, mean, variance, ks };
})();
if (typeof module !== 'undefined') module.exports = STO;
