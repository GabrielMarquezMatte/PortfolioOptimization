# PortOpt — Plano 2: objetivos de QP (UTILITY, IMPORTANCE, KELLY, TRACKING_ERROR, MAX_RETURN)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Os cinco objetivos que se resolvem com o QP do plano 1 — `UTILITY`, `IMPORTANCE`, `KELLY`, `TRACKING_ERROR` e `MAX_RETURN` (vol-alvo) — de ponta a ponta (núcleo → C → Python), com a matriz objetivo × setter da spec §3.3 validada no solve e as métricas `utility` e `tracking_error`.

**Architecture:** Um construtor de média-variância (`G = aΣ` com fator √a·L reaproveitado de `Data`, `g = −b(μ − rf)`) substitui o de `MIN_VARIANCE` e serve a MIN_VARIANCE (a=2, b=0), UTILITY (a=δ, b=1), IMPORTANCE (a=2(1−λ), b=λ) e KELLY (a=1/(1+rf), b=1). `TRACKING_ERROR` tem construtor próprio sobre os retornos em excesso. `MAX_RETURN` percorre o caminho paramétrico do QP `½w'Σw − t(μ−rf)'w`, que é afim por partes em t, e acha a vol-alvo exatamente no segmento certo. A validação de requisitos e da matriz §3.3 sai de `solve.cpp` para `problems/requirements.*`.

**Tech Stack:** o do plano 1 (C++23, CMake + presets, GoogleTest, Python + numpy + pytest + mypy, uv).

**Spec:** `docs/superpowers/specs/2026-11-15-portopt-design.md` (§3.2 setters, §3.3 matriz e notas por objetivo, §4.2 métodos, §5.2–5.3 contrato e métricas).

## Roteiro dos planos

1. Fundação + `MIN_VARIANCE` (concluído).
2. **Objetivos de QP** (este).
3. Homogeneização e leverage: `MAX_SHARPE`, `MAX_DIVERSIFICATION`, `UNBOUNDED`, família implícita do leverage; métricas `sharpe_ratio` e `diversification_ratio`.
4. `RISK_PARITY`: Spinu, Lagrangiano aumentado, `target_volatility`, `risk_budget`; métrica `risk_contributions`.
5. `cpp/portopt.hpp`, benchmarks, wheels/PyPI, release.

## Global Constraints

- Todas as do plano 1 continuam valendo (`docs/superpowers/plans/2026-11-15-portopt-01-fundacao.md`, seção Global Constraints), em especial: interface C só cresce por acréscimo (função nova, valor de enum novo, nunca parâmetro novo); `po_` / `PO_`; nenhuma exceção cruza a fronteira; sem gambiarra numérica; estilo do `STYLEGUIDE.md`; commits **sem** trailer `Co-Authored-By`.
- Positividade definida decidida uma vez: toda Hessiana da forma cΣ recebe `hessian_factor = √c·L`, com L = `Data::covariance_factor()`. Nenhum construtor chama `linalg::cholesky` sobre cΣ.
- Multiplicadores: KKT do problema **exatamente como escrito** abaixo, em w, convenção do header (∇f = λ_budget·1 + λ_ret·∇retorno + λ_vol·∇vol + Σβᵢeᵢ − β₀·1; cada multiplicador = ∂f*/∂(lado direito)). O objetivo QP de cada construtor é f menos uma constante, exceto no `MAX_RETURN`, que reescala (Task 3).
- Formas de minimização f (vão para o header e para `math.md`), com μ̃ = μ − rf·1 (rf = 0 sem caixa):
  - `UTILITY`: f = (δ/2) w'Σw − μ̃'w; δ padrão 5.
  - `IMPORTANCE`: f = (1−λ) w'Σw − λ μ̃'w; λ ∈ [0, 1), padrão 0.5.
  - `KELLY`: f = (1/(2(1+rf))) w'Σw − μ̃'w (a UTILITY com δ = 1/(1+rf); exige 1 + rf > 0).
  - `TRACKING_ERROR`: f = ‖R̃w + rf·1 − b‖², R̃ = R − rf·1·1' (sem caixa, R̃ = R e o termo rf some).
  - `MAX_RETURN`: f = −μ̃'w, sujeito a √(w'Σw) ≤ s.
- Matriz objetivo × setter (spec §3.3, sem as colunas de leverage e de risk parity, que são dos planos 3 e 4). Qualquer outra combinação é `Failure{InvalidArgument}`:

  | Objetivo | exige | aceita além de bounds |
  |---|---|---|
  | MIN_VARIANCE | covariance (+ mu se target_return) | target_return |
  | MAX_RETURN | mu, covariance, max_volatility | — |
  | UTILITY | mu, covariance | risk_aversion |
  | IMPORTANCE | mu, covariance | importance |
  | KELLY | mu, covariance; 1 + rf > 0 | — |
  | TRACKING_ERROR | returns (+ mu se target_return) | target_return |

  Mensagens: `"<setter> does not apply to <OBJECTIVE>"`, `"<OBJECTIVE> requires <item>"` (setters pelos nomes da API C sem prefixo: `target_return`, `risk_aversion`, `importance`, `max_volatility`; itens `mu`, `covariance`, `returns`, `max_volatility`; objetivos pelo sufixo do enum, ex. `UTILITY`); `"KELLY requires risk_free > -1"`.
- Ordem de validação em `solve`: objetivo implementado (os dos planos 3–4 seguem `Unsupported`) → setters fora da matriz → dados exigidos → tamanhos de bounds.

## Review Focus

1. **Escala diária no MAX_RETURN** (Σ ~ 1e-4, μ ~ 5e-4, s ~ 0.01): o caminho em t precisa achar o segmento e a raiz sem depender da escala. Teste na Task 3 (`MaxReturn.DailyScale`).
2. **Todos os μ iguais, sem caixa, no MAX_RETURN**: o retorno é constante no plano do orçamento; com s ≥ σ_min a resposta é a mínima variância, `OPTIMAL`, λ_vol = 0. Teste na Task 3 (`MaxReturn.EqualMeansIsMinVariance`).
3. **Retornos sem posto coluna completo no TRACKING_ERROR** (T < N, coluna duplicada, coluna constante igual a rf com caixa): `INVALID_ARGUMENT`, nunca `INTERNAL_ERROR`. Testes na Task 2.
4. **IMPORTANCE com λ perto de 1** (0.999999): Hessiana 2e-6·Σ, ainda `OPTIMAL` com KKT; λ = 1 é rejeitado no setter. Teste na Task 1 (`Importance.NearOne`).
5. **Python com `returns` sem `benchmark` (ou formatos errados)**: `ValueError` antes da chamada C. Teste na Task 5.

---

### Task 1: Construtor de média-variância, UTILITY, IMPORTANCE, KELLY e a matriz de requisitos

**Files:**
- Create: `src/problems/mean_variance.hpp`, `src/problems/mean_variance.cpp` (substituem `min_variance.*`, que são removidos)
- Create: `src/problems/requirements.hpp`, `src/problems/requirements.cpp` (os checks hoje em `solve.cpp` mais a matriz)
- Modify: `src/core/spec.hpp`, `src/problems/solve.cpp`, `CMakeLists.txt`
- Create: `tests/unit/user_kkt.hpp` (extraído de `min_variance_test.cpp`, generalizado)
- Modify: `tests/unit/problems/min_variance_test.cpp` (usa o helper e o construtor novo; nenhum teste muda de valor)
- Test: `tests/unit/problems/utility_test.cpp`, `tests/unit/problems/requirements_test.cpp`

**Interfaces:**
- Consumes: `po::Data` (`covariance_factor()`, `risk_free()`), `po::problems::add_constraint_rows`, `to_result`, `Translation` (plano 1).
- Produces:
  ```cpp
  // core/spec.hpp — campos novos
  std::optional<double> risk_aversion;   // UTILITY
  std::optional<double> importance;      // IMPORTANCE
  std::optional<double> max_volatility;  // MAX_RETURN (usado na Task 3)

  namespace po::problems {
  // G = a Sigma (factor sqrt(a) L), g = -b (mu - rf 1); constraint rows from add_constraint_rows.
  struct MeanVariance { double quadratic; double linear; };
  MeanVariance mean_variance_of(const Spec& spec, const Data& data);   // a, b of MIN_VARIANCE/UTILITY/IMPORTANCE/KELLY
  Translation mean_variance_problem(const Spec& spec, const Data& data, MeanVariance scales);
  Result solve_mean_variance(const Spec& spec, const Data& data);
  // requirements.hpp
  std::optional<Failure> check_request(const Spec& spec, const Data& data);  // implemented, matrix, data, sizes
  }
  namespace po::testing {   // tests/unit/user_kkt.hpp
  // Checks, independently of the QP: equalities, bounds (per-row relative complementarity, as in plan 1), and
  // grad f = budget*1 + target_return*(mu - rf) + volatility*Sigma w / vol + beta - beta0*1, where each term is
  // skipped when its multiplier is NaN. objective_gradient is grad f(w) over the risky weights, computed by the test.
  void expect_user_kkt(const Spec& spec, const Data& data, const Result& result,
                       std::span<const double> objective_gradient, double tolerance = 1e-10);
  }
  ```
  `mean_variance_of`: MIN_VARIANCE (2, 0); UTILITY (δ, 1); IMPORTANCE (2(1−λ), λ); KELLY (1/(1+rf), 1). Defaults δ = 5, λ = 0.5, rf = 0 sem caixa. `MIN_VARIANCE` passa a ser resolvido por `solve_mean_variance` (mesmo QP do plano 1, bit a bit).

- [ ] **Step 1: Testes** (Σ, μ de `random_covariance`/`random_vector`; todo caso `OPTIMAL` chama `expect_user_kkt` com ∇f calculado no teste)
  - `Utility.UnconstrainedMatchesClosedForm`: n = 6, δ = 3, sem caixa → w = (1/δ)Σ⁻¹μ + ((1 − 1'Σ⁻¹μ/δ)/1'Σ⁻¹1)Σ⁻¹1 em 1e-12 relativo.
  - `Utility.CashMatchesClosedForm`: rf = 0.02, δ = 4 → pesos de risco Σ⁻¹(μ − rf)/δ em 1e-12; caixa = 1 − soma; `budget_multiplier` NaN.
  - `Utility.DefaultRiskAversionIsFive`: sem `risk_aversion` igual a `risk_aversion = 5` bit a bit.
  - `Utility.LongOnlyKkt`: `lower = 0`, `upper = 1`, n = 20, sementes 1..20.
  - `Utility.LargeRiskAversionApproachesMinVariance`: δ = 1e8 → pesos da MIN_VARIANCE em 1e-6.
  - `Importance.ZeroIsMinVariance`: λ = 0 → mesmo resultado da MIN_VARIANCE bit a bit (mesmo QP).
  - `Importance.MatchesUtility`: λ = 0.3 → pesos da UTILITY com δ = 2(1−λ)/λ em 1e-12; multiplicadores = λ × os da UTILITY em 1e-10.
  - `Importance.NearOne`: λ = 0.999999, long-only → `OPTIMAL`, KKT.
  - `Kelly.WithoutCashIsUtilityOne` e `Kelly.WithCashIsUtilityOfRate`: pesos iguais aos da UTILITY com δ = 1 e δ = 1/(1+rf) em 1e-12; com caixa e sem bounds, pesos de risco = (1+rf)Σ⁻¹(μ − rf) em 1e-12.
  - `Requirements.*`: um teste por linha da matriz e por mensagem — `MinVarianceRejectsRiskAversion`, `UtilityRejectsTargetReturn`, `UtilityRequiresMu`, `ImportanceRejectsMaxVolatility`, `KellyRejectsImportance`, `KellyRequiresRateAboveMinusOne` (rf = −1), `UnimplementedObjectivesStayUnsupported` (MAX_SHARPE, MAX_DIVERSIFICATION, RISK_PARITY), e mensagens exatas da Global Constraints.
  - Os testes `MinVariance.*` existentes passam sem mudar valores.
- [ ] **Step 2:** Rodar, ver falhar.
- [ ] **Step 3:** Implementar; `solve.cpp` fica só com o despacho por objetivo.
- [ ] **Step 4:** `ctest --preset gcc-release` e `gcc-debug` PASS.
- [ ] **Step 5: Commit** — `feat: mean-variance builder with UTILITY, IMPORTANCE and KELLY; request checks`

---

### Task 2: TRACKING_ERROR

**Files:**
- Modify: `src/core/data.hpp`, `src/core/data.cpp`
- Create: `src/problems/tracking_error.hpp`, `src/problems/tracking_error.cpp`
- Modify: `src/problems/solve.cpp`, `src/problems/requirements.cpp`, `CMakeLists.txt`
- Test: `tests/unit/core/data_test.cpp`, `tests/unit/problems/tracking_error_test.cpp`

**Interfaces:**
- Consumes: Task 1 (`check_request`, `expect_user_kkt`, `Translation`, `add_constraint_rows`, `to_result`); `linalg::cholesky` (pivô relativo do plano 1).
- Produces:
  ```cpp
  // core/data.hpp
  // Copies returns (t x n, row-major) and benchmark (t). Failure when t = 0, sizes mismatch or a value is not
  // finite. Replaces earlier returns.
  std::expected<void, Failure> set_returns(std::span<const double> returns, std::span<const double> benchmark,
                                           std::size_t t);
  bool has_returns() const;
  std::size_t observation_count() const;          // T
  std::span<const double> returns() const;        // T x N
  std::span<const double> benchmark() const;      // T
  namespace po::problems {
  // G = 2 R~'R~, g = -2 R~'(b - rf 1), R~ = R - rf 1 1'. Failure{InvalidArgument,
  // "returns minus the risk-free rate do not have full column rank"} when cholesky(R~'R~) fails.
  std::expected<Translation, Failure> tracking_error_problem(const Spec& spec, const Data& data);
  std::expected<Result, Failure> solve_tracking_error(const Spec& spec, const Data& data);
  }
  ```
  O fator passado ao QP é √2 vezes o de R̃'R̃ (uma única decisão de posto, por `linalg::cholesky`). A covariância não é usada; `Data` pode ter `cov` ausente. `solve` passa a devolver o `Failure` do construtor.

- [ ] **Step 1: Testes** (R de `random_vector` em T×N, sementes fixas)
  - `Data.SetReturnsCopiesAndValidates`: copia; t = 0, tamanhos errados e NaN → `Failure`.
  - `TrackingError.RecoversExactReplication`: T = 60, N = 5, w* com 1'w* = 1, b = Rw* → pesos = w* em 1e-10; f* ≈ 0.
  - `TrackingError.MatchesConstrainedLeastSquares`: b com ruído → w = M⁻¹R'b + ((1 − 1'M⁻¹R'b)/1'M⁻¹1)M⁻¹1, M = R'R, em 1e-10.
  - `TrackingError.CashReplication`: rf = 0.001, b = R̃w* + rf·1 → pesos de risco = w* em 1e-10, caixa = 1 − 1'w*.
  - `TrackingError.CashWithZeroRateIsNotSingular`: rf = 0 com caixa → `OPTIMAL` (o bug da coluna constante não volta).
  - `TrackingError.TargetReturnAndLongOnlyKkt`: μ, `target_return`, `lower = 0` → KKT com ∇f = 2R̃'(R̃w + rf·1 − b).
  - `TrackingError.RankDeficientIsInvalid`: T < N; uma coluna duplicada; com caixa, uma coluna constante igual a rf → `Failure{InvalidArgument}` com a mensagem acima.
  - `TrackingError.WithoutCovariance`: `Data::create(n, {}, {})` + `set_returns` → `OPTIMAL`.
  - `Requirements.TrackingErrorRequiresReturns`, `Requirements.TrackingErrorTargetReturnRequiresMu`.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** PASS em ambos os presets.
- [ ] **Step 5: Commit** — `feat: TRACKING_ERROR over excess returns`

---

### Task 3: MAX_RETURN com vol-alvo

**Files:**
- Modify: `src/qp/problem.hpp`, `src/qp/goldfarb_idnani.cpp` (o conjunto ativo na solução)
- Modify: `src/problems/constraints.hpp`, `src/problems/constraints.cpp` (escala dos multiplicadores)
- Create: `src/problems/max_return.hpp`, `src/problems/max_return.cpp`
- Modify: `src/problems/solve.cpp`, `src/problems/requirements.cpp`, `CMakeLists.txt`
- Test: `tests/unit/qp/goldfarb_idnani_test.cpp`, `tests/unit/problems/max_return_test.cpp`

**Interfaces:**
- Consumes: Tasks 1–2; `Spec::max_volatility`.
- Produces:
  ```cpp
  // qp::Solution — campo novo: inequality rows active at the solution, ascending; empty unless Optimal.
  std::vector<std::size_t> active_inequalities;
  // problems::to_result — parâmetro novo: the QP multipliers are divided by multiplier_scale (the QP objective
  // is multiplier_scale * f plus a constant).
  Result to_result(const Data& data, const ConstraintRows& rows, const qp::Solution& solution,
                   double multiplier_scale = 1.0);
  namespace po::problems { Result solve_max_return(const Spec& spec, const Data& data); }
  ```

  Algoritmo (o QP paramétrico é Q(t): min ½w'Σw − t μ̃'w com as linhas de `add_constraint_rows`, G = Σ, fator L; w(t) é afim por partes e √(w(t)'Σw(t)) é não decrescente em t):
  1. Q(0) é a mínima variância. `Infeasible` → status `INFEASIBLE`. σ_min = vol(w(0)). Se s < σ_min·(1 − 1e-12) → `INFEASIBLE`. Se |s − σ_min| ≤ 1e-12·s → resposta w(0), λ_vol = 0 (vale também com todos os μ̃ iguais no plano do orçamento).
  2. Dobrar t a partir de t₁ = σ_min² / max(‖μ̃‖∞, ε) (t tem unidade de variância sobre retorno) até vol(w(t)) ≥ s. Se dois t consecutivos têm o mesmo conjunto ativo e o mesmo w (relativo 1e-14), a restrição de vol é inativa: a resposta é esse w, com λ_vol = 0 e os demais multiplicadores = limite de u(t)/t, obtido de duas soluções no mesmo segmento: (u(t₂) − u(t₁))/(t₂ − t₁). Limite de 200 dobras → `MAX_ITER`.
  3. Com o colchete [t_a, t_b], bissectar em t até que t_a e t_b tenham o mesmo conjunto ativo (limite 200 passos → `MAX_ITER`).
  4. No segmento, w(t) = w_a + (t − t_a)(w_b − w_a)/(t_b − t_a); resolver vol(w(t))² = s² (quadrática em t, raiz em [t_a, t_b]) e resolver Q(t*) uma última vez para o ponto e os multiplicadores.
  5. Multiplicadores em f = −μ̃'w: os do QP divididos por t* (`multiplier_scale = t*`) e λ_vol = −s/t* (≤ 0, = ∂f*/∂s).
  6. `iterations` = soma das iterações de todos os QPs resolvidos.

- [ ] **Step 1: Testes** (todo `OPTIMAL` chama `expect_user_kkt` com ∇f = −μ̃)
  - `GoldfarbIdnani.ReportsActiveInequalities`: problema com 2 bounds ativos → exatamente esses índices; vazio quando não `Optimal`.
  - `MaxReturn.UnconstrainedIsUpperFrontier`: sem caixa nem bounds, s acima de σ_min → retorno r* com σ²(r*) = s² e r* ≥ B/A (parábola do plano 1) em 1e-10; pesos iguais aos da MIN_VARIANCE com `target_return = r*` em 1e-9.
  - `MaxReturn.CashIsCapitalMarketLine`: rf = 0.02 → pesos de risco (s/√(μ̃'Σ⁻¹μ̃))Σ⁻¹μ̃ em 1e-10.
  - `MaxReturn.BelowMinimumIsInfeasible`, `MaxReturn.AtMinimumIsMinVariance`.
  - `MaxReturn.InactiveVolatilityLongOnly`: `lower = 0`, s = 10·max σᵢ → todo o peso no ativo de maior μ, `OPTIMAL`, λ_vol = 0, KKT.
  - `MaxReturn.MultiplierIsDerivativeOfOptimum`: long-only com vol ativa; (f*(s+h) − f*(s−h))/(2h), h = 1e-6·s, igual a λ_vol em 1e-5 relativo.
  - `MaxReturn.DailyScale`: Σ·1e-4, μ·1e-2, s = 1.2·σ_min → vol(w) = s em 1e-12 relativo, KKT.
  - `MaxReturn.EqualMeansIsMinVariance`: μ = 0.1·1, sem caixa, s = 2σ_min → pesos da MIN_VARIANCE, λ_vol = 0.
  - `MaxReturn.RandomLongOnly`: sementes 1..30, n em 3..25, s entre σ_min e max σᵢ → vol = s (ou λ_vol = 0 com vol < s), KKT.
  - `Requirements.MaxReturnRequiresMaxVolatility`, `Requirements.MaxVolatilityOnlyForMaxReturn`.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** PASS em ambos os presets.
- [ ] **Step 5: Commit** — `feat: MAX_RETURN under a volatility cap along the parametric QP path`

---

### Task 4: Interface C e métricas

**Files:**
- Modify: `include/portopt/portopt.h`, `src/capi/spec_api.cpp`, `src/capi/data_api.cpp`, `src/capi/metrics_api.cpp`, `src/metrics/portfolio.hpp`, `src/metrics/portfolio.cpp`
- Test: `tests/unit/capi/spec_test.cpp`, `tests/unit/capi/data_test.cpp`, `tests/unit/capi/solve_test.cpp`, `tests/unit/capi/metrics_test.cpp`, `tests/unit/metrics/portfolio_test.cpp`

**Interfaces:**
- Consumes: Tasks 1–3.
- Produces (header; nada existente muda de assinatura):
  ```c
  PO_API po_status po_spec_set_risk_aversion(po_spec* s, double delta);   /* delta > 0, finite */
  PO_API po_status po_spec_set_importance(po_spec* s, double lambda);     /* 0 <= lambda < 1 */
  PO_API po_status po_spec_set_max_volatility(po_spec* s, double sigma);  /* sigma > 0, finite */
  /* Copies returns (t x n, row-major, n = asset count) and benchmark (t). Replaces earlier returns. */
  PO_API po_status po_data_set_returns(po_data* d, const double* returns, const double* benchmark, size_t t);
  PO_API po_status po_utility(const po_data* d, const double* w, size_t n, double delta, double* out);
  /* sqrt(sum_t (r_t'w + w_0 rf - b_t)^2 / T), the root mean square the TRACKING_ERROR objective minimizes. */
  PO_API po_status po_tracking_error(const po_data* d, const double* w, size_t n, double* out);
  ```
  ```cpp
  namespace po::metrics {
  double utility(const Data& data, std::span<const double> w, double delta);   // return - delta/2 w'Sigma w
  double tracking_error(const Data& data, std::span<const double> w);
  }
  ```
  O header passa a listar os objetivos resolvidos (MIN_VARIANCE, MAX_RETURN, UTILITY, IMPORTANCE, KELLY, TRACKING_ERROR; `PO_UNSUPPORTED` para os demais), a forma f de cada um (Global Constraints) para a convenção de multiplicadores, o termo λ_vol·Σw/vol com `PO_CONSTRAINT_VOLATILITY`, e a matriz de setters por objetivo.

- [ ] **Step 1: Testes**
  - `CapiSpec.RiskAversionDomain`, `CapiSpec.ImportanceDomain` (0 aceito; 1, −0.1, NaN rejeitados), `CapiSpec.MaxVolatilityDomain` — valor inválido deixa a spec inalterada.
  - `CapiData.SetReturnsCopies` (alterar o array depois não muda o solve), `CapiData.SetReturnsRejectsNaNAndZeroT`, `CapiData.NullReturnsIsInvalid`.
  - `CapiSolve.EachQpObjectiveEndToEnd`: os cinco objetivos novos dão o mesmo resultado que `po::solve` bit a bit; `PO_CONSTRAINT_VOLATILITY` não NaN só no MAX_RETURN com vol ativa.
  - `CapiSolve.SetterOutsideMatrixIsInvalid`: UTILITY + `target_return` → `PO_INVALID_ARGUMENT`, mensagem da Global Constraints.
  - `Metrics.UtilityAndTrackingError` (valores à mão; caixa inclui w₀·rf), `CapiMetrics.TrackingErrorWithoutReturnsIsInvalid`, `Metrics.TrackingErrorAtOptimumIsRootOfObjective` (= √(f*/T)).
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** `ctest` PASS (inclui `c_header_check` e `style`).
- [ ] **Step 5: Commit** — `feat(capi): QP objectives, returns data and utility/tracking-error metrics`

---

### Task 5: Python e documentos

**Files:**
- Modify: `python/src/portopt/_native.py`, `_spec.py`, `_data.py`, `_result.py` (multiplicador de vol)
- Test: `python/tests/test_qp_objectives.py`
- Modify: `docs/guide/math.md`, `docs/guide/validation.md`, `ARCHITECTURE.md`, `README.md`, `CHANGELOG.md`

**Interfaces:**
- Consumes: Task 4.
- Produces:
  ```python
  @dataclass(frozen=True)
  class Spec:
      objective: Objective
      _: KW_ONLY
      bounds: tuple[float | ArrayLike, float | ArrayLike] | None = None
      target_return: float | None = None
      risk_aversion: float | None = None
      importance: float | None = None
      max_volatility: float | None = None
      # __post_init__ builds the po_spec handle through the setters (an invalid value raises ValueError at
      # construction) and keeps it in a private field(init=False, repr=False, compare=False) set with
      # object.__setattr__; dataclasses.replace(spec, target_return=r) builds a new, validated spec.
  class Data:
      def __init__(self, mu, cov, *, risk_free: float | None = None,
                   returns: ArrayLike | None = None, benchmark: ArrayLike | None = None) -> None
      def utility(self, weights: ArrayLike, risk_aversion: float) -> float
      def tracking_error(self, weights: ArrayLike) -> float
  @dataclass(frozen=True)
  class Result: ...; volatility_multiplier: float
  ```
  `Spec` vira dataclass congelada (decisão do autor, 2026-11-15); o exemplo de fronteira do README e de `test_min_variance.py` passa a usar `dataclasses.replace`. `Data` continua classe (carrega arrays grandes). `returns` 2-D (T × N, N = asset count) e `benchmark` 1-D de tamanho T, os dois ou nenhum; senão `ValueError` antes da chamada C. Com `mu` e `cov` ausentes, N vem de `returns.shape[1]`.

- [ ] **Step 1: Testes** — `test_spec_is_frozen_dataclass` (atribuição levanta `FrozenInstanceError`; duas specs com os mesmos campos são iguais; `repr` mostra os campos e não o handle), `test_spec_replace_builds_validated_spec` (`replace(spec, target_return=r)` resolve; `replace(spec, importance=1.0)` levanta `ValueError` na construção), `test_utility_matches_closed_form` (numpy), `test_importance_zero_is_min_variance`, `test_kelly_is_utility_of_rate`, `test_tracking_error_replication` (sem `cov`), `test_max_return_on_capital_market_line`, `test_max_return_below_minimum_is_infeasible`, `test_setter_outside_matrix_raises_value_error`, `test_returns_without_benchmark_raises`, `test_returns_wrong_width_raises`, `test_metrics_utility_and_tracking_error`.
- [ ] **Step 2:** Rodar, ver falhar: `PORTOPT_LIB=$PWD/build/gcc-release/libportopt.dll uv run --project python pytest python/tests -q`.
- [ ] **Step 3:** Implementar o pacote.
- [ ] **Step 4:** `math.md`: uma seção por objetivo com f, o QP construído (G, g, fator), a matriz de setters e, no MAX_RETURN, o caminho paramétrico e λ_vol = −s/t*. `validation.md`: os testes de consistência novos (IMPORTANCE↔UTILITY, KELLY↔UTILITY, MAX_RETURN↔fronteira e CML, derivada de f* por diferenças finitas). `ARCHITECTURE.md`: `mean_variance`, `tracking_error`, `max_return`, `requirements` na tabela e na anatomia do solve. `README.md`: um exemplo de MAX_RETURN. `CHANGELOG.md` em `[Unreleased]`.
- [ ] **Step 5:** Verificar: pytest PASS, `uv run --project python mypy --config-file python/pyproject.toml python/src/portopt` limpo, `python tests/style/check_style.py` limpo, `ctest --preset gcc-release` PASS.
- [ ] **Step 6: Commit** — `feat(python): QP objectives; docs for plan 2`
