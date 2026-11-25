# giac latex printer (gen2markup) -- Python prototype notes (1.6/T4.c1, CAS-615)
Source markup.cc (C0009) (giac 2.0.0.21); port markup.py (latex arms), parse.py, check.py: 301/301
rows = giac latex() via libgiac (probe.cas -> probe.out -> table.tsv). markup.py cites markup.cc per function.
## Entry and block
- export_latex :3629: flags TOPLEVEL|LATEX|ELEMPOW; gen2markup; prepend_minus. TOPLEVEL/FACTOR stripped (:1243).
- MarkupBlock :179: priority (0 = atom), type bits :83-101, neg (pending minus, not yet in text), latex,
  split_pos/split_pos_tex (elementary-function split point), appl (implicit application `\sin x`).
- Priorities :65-80. Wrap rule: child wrapped (`\left(..\right)`, priority -> 0) when child.priority >=
  parent's, at each call site (comparisons differ per site: pow uses `>=EXP || left.appl`; limit/sum use `>`).
- prepend_minus :972: only if neg; parenthesize first when priority > MUL; "-" prefix; priority=MUL,
  type=NEG (other type bits lost), neg=false. Parents call it on children, except prod/inv/division/
  neg/abs/cplx-im which carry `neg` up instead.
## Atoms
- Int :1323: |v| text, neg = v<0, type NUMERIC_EXACT|LEADING_DIGIT. Bool :1308: idnt2markup("true")
  -> `\mathrm{true}`, type SUBTYPE_IDNT|LEADING_DIGIT.
- Double :1362: text = giac print of |v| (12 sig. digits, `1e-05`). If text has `e` with nonzero
  mantissa and exponent: strip exponent leading zeros, add `.0` to an integer mantissa,
  `mant\times10^{-ex}`, MUL: `1e-10`->`1.0\times10^{-10}`, `1.23e+20`->`1.23\times10^{20}`; `12345678901.0`->`12345678901`.
- Cplx :1408: im 0 -> re. re 0: |im|=1 -> `\mathrm{i}` (neg if -1); else im block (FACTOR) wrapped when
  priority>=MUL, `im \mathrm{i}`, priority MUL. Else priority ADD: re (minus prepended) + `+`/`-` +
  (|im|=1 ? "" : |im| block + " ") + `\mathrm{i}`; `1+2 \mathrm{i}`, `2-3 \mathrm{i}`.
- Frac :1463: `\frac{|num|}{|den|}`, FRACTION, MUL, neg by sign. Gaussian num -> division(num,den).
- String :408/:333: trim spaces/tabs/newlines, escape `\ % _ $ { } ^` with `\`, `\text{``s''}`;
  several lines -> `\begin{array}{l}\text{..}\\\text{..}\end{array}`, blank inner lines dropped.
- Idnt :1486: infinity `\infty ` (trailing space), pi `\pi `, euler_gamma `\gamma `; else idnt2markup;
  HAS_SUBSCRIPT if `_{` in text; IDNT_NAME if name length > 1. undef -> `\mathrm{undef}`.
## idnt2markup :621-934 (tex, non-unit)
1. `greek_` -> `\greek\_`. 2. trailing digits split as index when stem is `log`, one char, greek or
double letter (`x12` -> `x_{12}`, `alpha1` -> `\alpha_{1}`). 3. else first `_` splits when stem is
one char/greek/double or index all digits (`abc_1` -> `\mathrm{abc}_{1}`, `a_bc` -> `a_{\mathrm{bc}}`).
4. greek stem/index -> `\name` (phi, epsilon -> `\varphi`, `\varepsilon`); double letter `xx` ->
`\mathbf{x}`; stem length > 1 -> `\mathrm{..}`. 5. index split on `_`; all tokens numeric -> joined `,`
(`x_1_2` -> `x_{1,2}`); index kept raw if one char/greek or atof != 0, else `\mathrm{..}`. 6. no index:
`_` -> `\_`. atof is C atof (prefix; accepts inf/nan). Greek list :141 (no `iota`, `upsilon`, `omicron`).
## Vectors :1658-1773
- Matrix (ckmatrix: rows same length > 0, no nested vect): `\left[\begin{array}{cc}a&b\\c&d\end{array}\right]`.
- Else items minus-prepended, joined `,`: list/poly1 `\left[..\right]`, set `\left\{..\right\}`,
  empty set `\emptyset `, seq `\left(..\right)` (bare at TOPLEVEL). Empty list `\left[\right]`.
## Operators (_SYMB, dispatch order :1777-3611)
- neg :1777: child block with neg flipped (no text). Distributes visually: `-x-1`, `-x y`.
- plus :1793: n-ary (flattened), priority ADD; `+` only before non-neg terms; each term prepend_minus
  then wrapped if >= ADD (`-\left(x+y\right)` when a term is neg of a sum). Unary plus (non-vect arg)
  :1827: `+` prefix if atom, priority ADD (`x\to +\infty `).
- pointminus :1848: `a\ominus b` (right.neg branch is dead: get_leftright already prepended it).
- prod :1869-2032 (the hard one): flatten; scanning from the END, ints and fracs -> ncnst, pi/euler_gamma
  -> cnsts (both therefore reversed), pure-imaginary c (im != 1): im -> cnsts (im=-1 flips ml.neg) and c
  replaced by i in place. args = ncnst+cnsts+rest (`sqrt(3)*pi` -> `\pi  \sqrt{3}`, two spaces).
  Leading fraction :1897-1937: if the product has exactly one inv(...) factor, its arg an int,
  and num_count > 1 (`x*inv(2)*inv(3)` -> `\frac{x}{2\cdot 3}`): remove it and the first int factor (if any) and put fraction(n,d) (NOT reduced:
  `2*inv(4)*x` -> `\frac{2}{4} x`) in front. Each factor: inv -> denominator; child neg counted (odd count
  -> ml.neg), never printed in place; wrapped when its side has > 1 factor and (priority > MUL or ==MUL
  and not cplx) and not FRACTION and (not SUMPROD or not last on its side). Separator `\cdot ` when the
  factor is LEADING_DIGIT/SUBTYPE_IDNT/FACTORIAL/IDNT_NAME, or the previous one on that side is
  IDNT_NAME/SUBTYPE_IDNT, or both FRACTION, or previous is ELEMAPP with appl; else " ".
  (`x\cdot 2^{y}`, `\sin x\cdot x`, `x \sin x`, `k\mathclose{!}\cdot \left(n-k\right)\mathclose{!}`).
  Empty numerator `1`; with a denominator `\frac{num}{den}` (FRACTION). Same-side order kept.
- inv :2082: `\frac1{X}` for ANY X (`\frac1{x+1}`, `\frac1{2 x}`), neg carried, FRACTION, MUL.
- division :2097: `\frac{a}{b}`, neg = a.neg xor b.neg, FRACTION, MUL.
- exp :2148: exp(1) -> `\mathrm{e}` (priority 0); else `\mathrm{e}^{X}` (minus prepended), EXP, POWER.
- sqrt :2191: `\sqrt{X}`, UNARY, ROOT. pow :2224: left minus-prepended; exponent with 2*e == 1 (Frac
  1/2 or double .5, NOT a symbolic 1/2) -> `\sqrt{}`; ELEMPOW :2248: pow, right priority 0, |e| != 1,
  left ELEMAPP -> insert `^{e}` at split_pos_tex: `\cos ^{2}x`, `\sin ^{2}\left(x+1\right)`; result ELEMAPP
  with left's appl. Else `left^{right}` (right never wrapped: `x^{y^{z}}`, `\sqrt{x}^{-1}`), left wrapped
  if appl or priority >= EXP (`\left(\sin x\right)^{n+1}`, `\left(-x\right)^{2}`); POWER|LEADING_DIGIT of left.
- abs/floor/ceil :2281: `\left|X\right|`, `\left\lfloor X\right\rfloor `, `\left\lceil X\right\rceil `
  (trailing spaces). floor/ceil prepend the child's minus; abs DROPS it: quote(abs(-x)) -> `\left|x\right|`.
- integrate/sum/product :2375: priority MUL, SUMPROD. var from args[1] (or `k=a..b` form); integral
  var -> `\,\mathrm{d}x` appended after the body; bounds `_{lb}^{ub}` with `\int `/`\sum `/`\prod `
  (space before `_`): `\int _{0}^{1}f\left(x\right)\,\mathrm{d}x`, `\sum _{k=1}^{n}..`, `\sum _{k}..`.
  sum/int carry body neg up, wrap body if NEG type or priority > MUL; product prepends, wraps if
  not FRACTION and >= MUL.
- limit :2488: `\lim_{x\to a}` + `^{+}`/`^{-}` for dir +-1 (4th arg, or 3rd after `x=a`); body wrapped
  if priority > MUL (`\lim_{x\to 0}\left(f\left(x\right)+1\right)`).
- derive :2560: (a) partial `derive(f,k)` chains of an idnt -> `f^{\left(k1+1,..\right)}` (seq vect),
  type DERIVATIVE so `of` does not wrap it: `f^{\left(1,2\right)}\left(x,y\right)`. (b) unary derive
  nesting n: `f'`, `f''`, `f'''`, `f^{(4)}`. (c) `derive(e,vars..)`: get_derive_vars :1069 merges
  repeats; one var once: `\frac{\mathrm{d}}{\mathrm{d}x}e`, or `\frac{\mathrm{d}e}{\mathrm{d}x}` when e
  is an atom (priority 0); else `\frac{D^{N}}{D x^{k} D y}e` with D=`\mathrm{d}` (1 var) or `\partial `
  (>1 vars), factors joined " ", atom e goes in the numerator.
- equal `=` / same `==` / different `\neq ` / sto `:=` (args reversed: sto(value,var)) :2728; priorities
  EQ/COMP/ASSGN; sides wrapped when >=. Inequalities :2770 `<`, `\leq `, `>`, `\geq ` (INEQ).
- interval :2942 `a\ldots b` (SET). factorial :3317 `X\mathclose{!}`, wrap if >= UNARY. at :3434 `x_{1}`.
- General :3532: priority APPLY. `of(f,args)` (unknown f) prints f's block; named functions use
  func2markup :1128 (`\sin `.. with trailing space, `\operatorname{arccot}`, `\operatorname{arsinh}`,
  Dirac `\operatorname{\delta}`, Heaviside `\operatorname{\theta}`, `\Gamma`, `\Psi`; others via
  idnt2markup: `\mathrm{erf}`, `\mathrm{rootof}`; one-char result stays bare). Several args (seq):
  `\left(a,b\right)` and the args' minus is NOT prepended (giac bug: f(-x,y) -> `f\left(x,y\right)`).
  One arg: minus prepended; elementary functions (:1106 list) set ELEMAPP, split_pos_tex = name length,
  appl = arg is priority-0, not a number, not a boolean, not a multi-letter name -> no parens: `\sin x`,
  `\ln x`; else `\sin \left(x^{2}\right)`, `\sin \left(2\right)`.
## Parser shapes (parse.py) -- what the Rust builder must produce for printing to match
- Evaluated forms: `n/m` -> Frac; `1/b` -> inv(b); `a/b` -> prod(a,inv(b)) (flattened by prod);
  `-x` -> neg(x), numeric minus folded into the number; `a-b` -> plus(a,neg(b)); `2*i` -> Cplx(0,2),
  `1+2*i` -> Cplx(1,2); f(x) unknown -> of(f,x); `x[1]` -> at(x,1); `set[..]`, `poly1[..]` subtypes.
- Precedence in giac's parser: `/` binds tighter than `*` and unary minus tighter than `*`
  (quote(-x*y) prints `(-x)*y`; quote(2*x/3) -> `2 \frac{x}{3}`), `^` tighter than unary minus.
- quote(...) forms: `/` is division(a,b) (no Frac, no inv): quote(x^(1/2)) -> `x^{\frac{1}{2}}` (not sqrt),
  quote(1/x) -> `\frac{1}{x}`.
## Shape overrides (printed form is not the shape latex() sees; check.py SHAPE)
- caseval prints a top-level `x=1`, `a=b+1`, `[x,y]=[1,2]` as the rhs; latex() sees equal(..)
  (`[x,y]=[1,2]` -> list of equations `\left[x=1,y=2\right]`). caseval joins string newlines with ` | `.
- quote() prints prod(a,inv(b)) and division(a,b) both as `a/b`, but they print differently:
  quote(2*x*inv(3)) -> `\frac{2}{3} x` vs quote(2*x/3) -> `2 \frac{x}{3}` (C0037 trap); the typed shape
  is used for the 6 quote(..inv..) rows.
- `x*i/2` (printed `x*i/2`) gives `\frac{1}{2} x \mathrm{i}`: reproduced by prod(x,i,inv(2)) (leading
  fraction) or prod(x,Cplx(0,1/2)), NOT by prod(x,Frac(i,2)) (that gives `\frac{\mathrm{i}}{2} x`).
  So parse.py keeps Cplx/int as prod(Cplx,inv(int)); `i/2` -> `\frac{\mathrm{i}}{2}`, `(3+2*i)/6` ->
  `\frac{3+2 \mathrm{i}}{6}` either way. The Frac-with-gaussian-num arm (:1466) is ported, unexercised.
- giac eval rewrites (not the printer): x^(2/3)->(x^(1/3))^2, x^(-1/2)->sqrt(x)^-1, ceil->-floor(-x),
  asinh..->ln forms, binomial->factorials, x<y->y>x, epsilon->1e-12, x!=y->true, 1e10->integer.
## Quirks that matter for Rust
- Trailing spaces are output: `\pi `, `\infty `, `\gamma `, `\sin `, `\rfloor `, `\leq `, `\geq `, `\neq `, `\cdot `,
  `\ldots `, `\emptyset `, `\int `; `\pi  x` has two spaces.
- `\frac1{..}` (no braces on 1) only from inv; Frac/division/prod use `\frac{1}{..}`.
- Booleans count as SUBTYPE_IDNT|LEADING_DIGIT; numbers break implicit application (`\sin \left(2\right)`, from source, not probed).
- Leading-fraction gathering builds an unreduced fraction and only uses the FIRST int numerator.
## rootof / algebraic extensions (out of scope; answers only)
- latex(rootof([1,0],[1,0,1,1])) = `\mathrm{rootof}\left(\left[1,0\right],\left[1,0,1,1\right]\right)`, type
  expression (symbolic: generic-function arm, matched). factor(x^3+x+1,rootof(..)), normal(rootof(..)^3)
  also give symbolic `\mathrm{rootof}\left(\left[\left[-1,-1\right],..\right]\right)`; no _EXT reached (:3783 -> legacy tex.cc).
