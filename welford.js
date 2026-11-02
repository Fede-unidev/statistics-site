// Welford's online algorithm + helpers for HW5
const WEL = (() => {
  // Ex. 1: the whole algorithm. Core = the 3 statements in the loop.
  function welford(a) {
    let M = 0, S = 0, k = 0;
    for (const x of a) { k++; const d = x - M; M += d / k; S += d * (x - M); }
    return { k, mean: M, ss: S, s2: k > 1 ? S / (k - 1) : NaN };
  }
  // Ex. 2: score distributions on {0,...,10}
  const binom = (n, k) => { let r = 1; for (let i = 1; i <= k; i++) r = r * (n - k + i) / i; return r; };
  function pmfOf(dist, p) {
    const w = Array.from({ length: 11 }, (_, k) => dist === 'uniform' ? 1 : dist === 'binom' ? binom(10, k) * p ** k * (1 - p) ** (10 - k) : 0.7 ** k);
    const t = w.reduce((s, v) => s + v, 0); return w.map(v => v / t);
  }
  function moments(pmf) {
    const mu = pmf.reduce((s, q, k) => s + q * k, 0), c = e => pmf.reduce((s, q, k) => s + q * (k - mu) ** e, 0);
    return { mu, s2: c(2), mu4: c(4) };
  }
  const sampler = pmf => { const cdf = []; pmf.reduce((s, q, k) => (cdf[k] = s + q), 0); return () => { const u = Math.random(); for (let k = 0; k < 10; k++) if (u < cdf[k]) return k; return 10; }; };
  // m hackers x n systems: Welford per hacker -> sample mean and sample variance of the n scores
  function experiment(pmf, m, n) {
    const draw = sampler(pmf), means = new Float64Array(m), vars = new Float64Array(m);
    for (let j = 0; j < m; j++) {
      let M = 0, S = 0;
      for (let k = 1; k <= n; k++) { const x = draw(); const d = x - M; M += d / k; S += d * (x - M); }
      means[j] = M; vars[j] = S / (n - 1);
    }
    return { means, vars };
  }
  // exact theory: Var(mean) = s2/n ; Var(S^2) = mu4/n - s2^2 (n-3)/(n(n-1))
  const theory = (mom, n) => ({ sdMean: Math.sqrt(mom.s2 / n), sdVar: Math.sqrt(mom.mu4 / n - mom.s2 ** 2 * (n - 3) / (n * (n - 1))) });
  const mean = a => { let s = 0; for (const v of a) s += v; return s / a.length; };
  const sd = a => { const m = mean(a); let s = 0; for (const v of a) s += (v - m) ** 2; return Math.sqrt(s / (a.length - 1)); };
  return { welford, pmfOf, moments, sampler, experiment, theory, mean, sd };
})();
if (typeof module !== 'undefined') module.exports = WEL;
