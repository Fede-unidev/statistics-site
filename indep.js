// Contingency tables and the chi-square test of independence (HW7)
const IND = (() => {
  // ---- parsing and table construction ----
  function parse(text, header) {
    let lines = text.split(/\r?\n/).map(s => s.trim()).filter(Boolean); if (header) lines.shift();
    return lines.map(l => l.split(/[,;\t]/).map(s => s.trim())).filter(a => a.length >= 2 && a[0] !== '' && a[1] !== '').map(a => [a[0], a[1]]);
  }
  function build(pairs) {
    const xs = [...new Set(pairs.map(p => p[0]))], fs = [...new Set(pairs.map(p => p[1]))];
    const O = xs.map(() => fs.map(() => 0)); const xi = new Map(xs.map((v, i) => [v, i])), fi = new Map(fs.map((v, i) => [v, i]));
    for (const [x, f] of pairs) O[xi.get(x)][fi.get(f)]++;
    return { xs, fs, O };
  }
  // ---- statistics ----
  function analyze(t) {
    const { xs, fs, O } = t, r = xs.length, c = fs.length, rows = O.map(row => row.reduce((s, v) => s + v, 0)), cols = fs.map((_, j) => O.reduce((s, row) => s + row[j], 0)), n = rows.reduce((s, v) => s + v, 0);
    const E = O.map((row, i) => row.map((_, j) => rows[i] * cols[j] / n));
    const contrib = O.map((row, i) => row.map((o, j) => (o - E[i][j]) ** 2 / E[i][j]));
    const resid = O.map((row, i) => row.map((o, j) => (o - E[i][j]) / Math.sqrt(E[i][j])));
    const chi2 = contrib.reduce((s, row) => s + row.reduce((a, v) => a + v, 0), 0), df = (r - 1) * (c - 1);
    const lowE = E.flat().filter(v => v < 5).length / (r * c);
    return { r, c, n, rows, cols, E, contrib, resid, chi2, df, p: df > 0 ? chi2p(chi2, df) : NaN, V: Math.sqrt(chi2 / (n * Math.max(1, Math.min(r - 1, c - 1)))), lowE, minE: Math.min(...E.flat()) };
  }
  // ---- chi-square distribution: survival function via the regularized incomplete gamma ----
  function lgamma(z) {
    const g = 7, C = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
    z -= 1; let x = C[0]; for (let i = 1; i < g + 2; i++) x += C[i] / (z + i); const t = z + g + 0.5;
    return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
  }
  function gammaQ(a, x) { // upper regularized gamma Q(a,x)
    if (x <= 0) return 1;
    if (x < a + 1) { let ap = a, sum = 1 / a, del = sum; for (let n = 0; n < 1000; n++) { ap++; del *= x / ap; sum += del; if (Math.abs(del) < Math.abs(sum) * 1e-16) break; } return Math.max(0, 1 - sum * Math.exp(-x + a * Math.log(x) - lgamma(a))); }
    let b = x + 1 - a, cc = 1 / 1e-300, d = 1 / b, h = d;
    for (let i = 1; i < 1000; i++) { const an = -i * (i - a); b += 2; d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300; cc = b + an / cc; if (Math.abs(cc) < 1e-300) cc = 1e-300; d = 1 / d; const del = d * cc; h *= del; if (Math.abs(del - 1) < 1e-16) break; }
    return Math.exp(-x + a * Math.log(x) - lgamma(a)) * h;
  }
  const chi2p = (x, df) => gammaQ(df / 2, x / 2);
  function chi2crit(df, alpha) { let lo = 0, hi = 1; while (chi2p(hi, df) > alpha) hi *= 2; for (let i = 0; i < 100; i++) { const m = (lo + hi) / 2; chi2p(m, df) > alpha ? lo = m : hi = m; } return (lo + hi) / 2; }
  // ---- permutation test: shuffle the F labels, recompute chi-square ----
  function permP(pairs, chi2obs, B) {
    const xs = [...new Set(pairs.map(p => p[0]))], fs = [...new Set(pairs.map(p => p[1]))], n = pairs.length;
    const X = Int32Array.from(pairs, p => xs.indexOf(p[0])), F = Int32Array.from(pairs, p => fs.indexOf(p[1])), r = xs.length, c = fs.length;
    const rows = new Float64Array(r), cols = new Float64Array(c); for (let k = 0; k < n; k++) { rows[X[k]]++; cols[F[k]]++; }
    let ge = 0; const cnt = new Float64Array(r * c);
    for (let b = 0; b < B; b++) {
      for (let i = n - 1; i > 0; i--) { const j = Math.floor(Math.random() * (i + 1)); const t = F[i]; F[i] = F[j]; F[j] = t; }
      cnt.fill(0); for (let k = 0; k < n; k++) cnt[X[k] * c + F[k]]++;
      let s = 0; for (let i = 0; i < r; i++) for (let j = 0; j < c; j++) { const e = rows[i] * cols[j] / n; s += (cnt[i * c + j] - e) ** 2 / e; }
      if (s >= chi2obs - 1e-9) ge++;
    }
    return (ge + 1) / (B + 1);
  }
  // ---- data ----
  const mulberry = a => () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
  const toText = pairs => 'X,F\n' + pairs.map(p => p.join(',')).join('\n');
  // X uniform on k categories; P(F = i | X = i) gets extra weight s (s = 0: independent)
  function synthetic(n, s, seed, k = 3) {
    const rnd = mulberry(seed), xl = ['low', 'mid', 'high', 'top', 'max'].slice(0, k), fl = ['A', 'B', 'C', 'D', 'E'].slice(0, k), out = [];
    for (let m = 0; m < n; m++) { const i = Math.floor(rnd() * k); let j; if (rnd() < s) j = i; else j = Math.floor(rnd() * k); out.push([xl[i], fl[j]]); }
    return out;
  }
  function fromCounts(rowNames, colNames, M) { const out = []; M.forEach((row, i) => row.forEach((v, j) => { for (let k = 0; k < v; k++) out.push([rowNames[i], colNames[j]]); })); return out; }
  const hairEye = () => fromCounts(['Black', 'Brown', 'Red', 'Blond'], ['Brown', 'Blue', 'Hazel', 'Green'], [[68, 20, 15, 5], [119, 84, 54, 29], [26, 17, 14, 14], [7, 94, 10, 16]]);
  return { parse, build, analyze, gammaQ, chi2p, chi2crit, permP, synthetic, hairEye, toText, mulberry };
})();
if (typeof module !== 'undefined') module.exports = IND;
