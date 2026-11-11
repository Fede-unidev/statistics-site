// Standard normal CDF: several numerical methods (HW6, Ex. 2)
const NCDF = (() => {
  const SQ2PI = Math.sqrt(2 * Math.PI), SQPI = Math.sqrt(Math.PI);
  const pdf = t => Math.exp(-t * t / 2) / SQ2PI;
  // 1. numerical integration:  Phi(x) = 1/2 +- int_0^|x| pdf(t) dt
  function trapezoid(x, N) {
    const a = Math.abs(x), h = a / N; let s = 0.5 * (pdf(0) + pdf(a));
    for (let i = 1; i < N; i++) s += pdf(i * h);
    const I = s * h; return x >= 0 ? 0.5 + I : 0.5 - I;
  }
  function simpson(x, N) {
    N += N & 1; const a = Math.abs(x), h = a / N; let s = pdf(0) + pdf(a);
    for (let i = 1; i < N; i++) s += pdf(i * h) * (i & 1 ? 4 : 2);
    const I = s * h / 3; return x >= 0 ? 0.5 + I : 0.5 - I;
  }
  // 2a. Maclaurin series of erf (alternating terms: cancellation for large |x|)
  function erfSeries(z) {
    let term = z, sum = z;
    for (let n = 1; n < 1000; n++) { term *= -z * z / n; const add = term / (2 * n + 1); sum += add; if (Math.abs(add) < 1e-17 * Math.abs(sum)) break; }
    return 2 / SQPI * sum;
  }
  const viaErfSeries = x => 0.5 * (1 + erfSeries(x / Math.SQRT2));
  // 2b. series with positive terms: Phi(x) = 1/2 + pdf(x) (x + x^3/3 + x^5/(3.5) + ...)
  function positiveSeries(x) {
    if (x > 40) return 1; if (x < -40) return 0;
    let t = x, s = x; for (let i = 3; i < 2000; i += 2) { t *= x * x / i; s += t; if (Math.abs(t) < 1e-17 * Math.abs(s)) break; }
    return 0.5 + s * pdf(x);
  }
  // 3. continued fraction for the tail:  Q(x) = pdf(x) / (x + 1/(x + 2/(x + 3/(x + ...)))), good for x >= 3
  function tailCF(x, K) { let f = x; for (let k = K; k >= 1; k--) f = x + k / f; return pdf(x) / f; }
  function continuedFraction(x, K = 60) { const q = tailCF(Math.abs(x), K); return x >= 0 ? 1 - q : q; }
  // 4. Abramowitz & Stegun 7.1.26: erf(z) ~ 1 - (a1 t + ... + a5 t^5) e^{-z^2}, |error| < 1.5e-7
  function abramowitzStegun(x) {
    const z = Math.abs(x) / Math.SQRT2, t = 1 / (1 + 0.3275911 * z),
      e = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-z * z);
    return 0.5 * (1 + (x >= 0 ? e : -e));
  }
  // 5. rational approximation in double precision (Hart 1968, as implemented by West 2005): small tail computed directly
  function rational(x) {
    const ax = Math.abs(x); let c;
    if (ax > 37) c = 0; else {
      const ex = Math.exp(-ax * ax / 2);
      if (ax < 7.07106781186547) {
        let b = 3.52624965998911e-2 * ax + 0.700383064443688; b = b * ax + 6.37396220353165; b = b * ax + 33.912866078383;
        b = b * ax + 112.079291497871; b = b * ax + 221.213596169931; b = b * ax + 220.206867912376; c = ex * b;
        b = 8.83883476483184e-2 * ax + 1.75566716318264; b = b * ax + 16.064177579207; b = b * ax + 86.7807322029461;
        b = b * ax + 296.564248779674; b = b * ax + 637.333633378831; b = b * ax + 793.826512519948; b = b * ax + 440.413735824752; c /= b;
      } else {
        let b = ax + 0.65; b = ax + 4 / b; b = ax + 3 / b; b = ax + 2 / b; b = ax + 1 / b; c = ex / b / 2.506628274631;
      }
    }
    return x > 0 ? 1 - c : c;
  }
  return { pdf, trapezoid, simpson, viaErfSeries, positiveSeries, continuedFraction, abramowitzStegun, rational };
})();
if (typeof module !== 'undefined') module.exports = NCDF;
