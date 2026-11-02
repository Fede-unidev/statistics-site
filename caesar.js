// Caesar cipher + statistical cryptanalysis (letters A-Z only)
const CAE = (() => {
  const mod26 = x => ((x % 26) + 26) % 26;
  function shift(text, k) { return text.replace(/[a-z]/gi, c => { const b = c <= 'Z' ? 65 : 97; return String.fromCharCode(mod26(c.charCodeAt(0) - b + k) + b); }); }
  const enc = (t, k) => shift(t, k), dec = (t, k) => shift(t, -k);
  const letters = t => [...t.toUpperCase().replace(/[^A-Z]/g, '')].map(c => c.charCodeAt(0) - 65);
  // reference model built from a training corpus
  function model(text) {
    const L = letters(text), u = Array(26).fill(0), b = Array(676).fill(0), t = Array(17576).fill(0);
    L.forEach((x, i) => { u[x]++; if (i > 0) b[L[i - 1] * 26 + x]++; if (i > 1) t[(L[i - 2] * 26 + L[i - 1]) * 26 + x]++; });
    const n = L.length, lg = (a, s) => { const tot = a.reduce((p, v) => p + v, 0); return a.map(v => Math.log((v + s) / (tot + s * a.length))); };
    return { n, p: u.map(v => v / n), lu: lg(u, 0.5), lb: lg(b, 0.5), lt: lg(t, 0.1) };
  }
  const WORDS = ['THE', 'AND', 'THAT', 'HAVE', 'WITH', 'YOU', 'FOR', 'NOT', 'THIS', 'BUT', 'HIS', 'ARE'];
  // each method scores a candidate plaintext (array of 0..25); higher = better
  const scorers = {
    chi2: (P, M) => { const o = Array(26).fill(0); P.forEach(x => o[x]++); let s = 0; for (let i = 0; i < 26; i++) { const e = P.length * M.p[i]; s += (o[i] - e) ** 2 / e; } return -s; },
    unigram: (P, M) => P.reduce((s, x) => s + M.lu[x], 0),
    bigram: (P, M) => { let s = 0; for (let i = 1; i < P.length; i++) s += M.lb[P[i - 1] * 26 + P[i]]; return s; },
    trigram: (P, M) => { let s = 0; for (let i = 2; i < P.length; i++) s += M.lt[(P[i - 2] * 26 + P[i - 1]) * 26 + P[i]]; return s; },
    words: P => { const s = P.map(x => String.fromCharCode(65 + x)).join(''); return WORDS.reduce((a, w) => a + (s.split(w).length - 1), 0); },
  };
  // returns the recovered key (shift used to ENCRYPT) for the chosen method
  function crack(C, M, method) {
    if (method === 'peak') { const o = Array(26).fill(0); C.forEach(x => o[x]++); return mod26(o.indexOf(Math.max(...o)) - 4); } // most frequent letter assumed to be E
    let best = -Infinity, bk = 0;
    for (let k = 0; k < 26; k++) { const s = scorers[method](C.map(x => mod26(x - k)), M); if (s > best) { best = s; bk = k; } }
    return bk;
  }
  const methods = ['peak', 'chi2', 'unigram', 'bigram', 'trigram', 'words'];
  return { enc, dec, letters, model, crack, methods, scorers, mod26 };
})();
if (typeof module !== 'undefined') module.exports = CAE;
