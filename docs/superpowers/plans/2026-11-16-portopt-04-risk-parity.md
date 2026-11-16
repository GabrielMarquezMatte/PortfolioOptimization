# PortOpt — Plano 4: RISK_PARITY, target_volatility e risk_budget

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Entregar `RISK_PARITY` de ponta a ponta (núcleo → interface C → Python): Newton de Spinu para o caso sem restrições extras, Lagrangiano aumentado com bounds por L-BFGS projetado para `target_return`, `target_volatility` e bounds ativos, os setters `risk_budget` e `target_volatility`, o resíduo e a métrica `risk_contributions`, validados por KKT, formas fechadas e referências independentes por varredura.

**Architecture:**
- **Pasta nova `src/riskparity/`** (`spinu`, `objective`, `box_lbfgs`, `augmented_lagrangian`, `risk_parity`). `src/problems/` já tem 8 módulos (limite do STYLEGUIDE), então `solve.*` e `requirements.*` saem dela para `src/solve/`. Ordem das dependências: `core` ← `qp` ← `problems` ← `riskparity` ← `solve` ← `capi`.
- **Caminho rápido:** sem `target_return`, sem `target_volatility` e com a solução de Spinu dentro dos bounds, o resultado é a de Spinu, `OPTIMAL`, todos os multiplicadores 0 e resíduo no piso do arredondamento.
- **Caminho geral:** viabilidade e σ_min pela mínima variância (`problems::solve_mean_variance`), depois Lagrangiano aumentado partindo da solução de Spinu. As igualdades são orçamento, `target_return` e vol(w) = s. Os bounds ficam no subproblema, resolvido por L-BFGS projetado. A viabilidade das igualdades é restaurada por Gauss-Newton e os multiplicadores saem por mínimos quadrados (QR) no ponto final.
- **Local, não global:** o problema é não convexo. `OPTIMAL` significa ponto de KKT atingido (spec §3.3), partindo de Spinu; não se promete o mínimo global.

**Tech Stack:** o dos planos 1–3. Nenhuma dependência nova.

**Spec:** `docs/superpowers/specs/2026-11-15-portopt-design.md`. As seções que importam são §3.2 (`risk_budget`, `target_volatility`), §3.3 (linha e nota do RISK_PARITY), §4.2 (Spinu e Lagrangiano aumentado), §5.2 (contrato de status, `po_result_residual`, multiplicadores) e §5.3 (`po_risk_contributions`).

## Roteiro dos planos

1. Fundação + `MIN_VARIANCE` (concluído).
2. Objetivos de QP (concluído).
3. Correção do MAX_RETURN, leverage, homogeneização (concluído).
4. **`RISK_PARITY`** (este).
5. `cpp/portopt.hpp`, benchmarks, wheels/PyPI, release, preparação da 1.0.

## Global Constraints

- Todas as dos planos 1 a 3 continuam valendo (seções Global Constraints de `…-01-fundacao.md`, `…-02-objetivos-qp.md` e `…-03-homogeneizacao-leverage.md`). Em especial:
  - a interface C só cresce por acréscimo, e `PO_ABI_VERSION` fica em 0;
  - uma convenção de multiplicadores: cada multiplicador é ∂f*/∂(lado direito), e uma restrição "≤" tem multiplicador ≤ 0;
  - nenhuma exceção atravessa a interface C;
  - o estilo do `STYLEGUIDE.md` (nunca `else`, aninhamento ≤ 3, função ≤ 60 linhas, ≤ 5 parâmetros, ≤ 8 módulos por pasta);
  - warnings são erros;
  - commits **sem** trailer `Co-Authored-By`;
  - código e documentos do repositório em inglês, documentos de `docs/superpowers/` em português.
- **Tolerância é medida, nunca chutada.**
  - Toda tolerância de decisão nova no código fica relativa ao ruído do problema e justificada em `docs/guide/validation.md` com a medição que a sustenta.
  - Os fatores iniciais deste plano (Newton, KKT, tetos de iteração) são pontos de partida: a Task 5 os substitui pelo valor medido.
  - Um teste que não atinge a sua faixa é reportado com os resíduos. Nunca se afrouxa uma faixa.
- **Faixas de teste por condicionamento** (KKT e concordância com a referência, relativas), as mesmas do plano 3:

  | cond(Σ) | tolerância |
  |---|---|
  | ≤ 1e4 | 1e-10 |
  | 1e6 | 1e-9 |
  | 1e8 | 1e-7 |
- **Rodar testes só com o preset `gcc-release`.** Não rodar `gcc-debug` (preferência registrada do autor; o plano 3 pedia os dois, este não).
- **Formas.**
  - c(w) = w ∘ (Σw) / (w'Σw) e f(w) = Σᵢ (cᵢ − bᵢ)², com Σbᵢ = 1.
  - ∇f = 2 J'(c − b), com J = (diag(Σw) + diag(w)Σ − 2 c (Σw)') / (w'Σw).
  - c é homogênea de grau 0, então w'∇f(w) = 0. Com as restrições, isso dá a identidade λ_orçamento = −λ_retorno·r − λ_vol·s − Σᵢ βᵢwᵢ.
  - Spinu: F(y) = ½ y'Σy − b' log y, y > 0, ∇F = Σy − b/y e ∇²F = Σ + diag(b/y²), estritamente convexa e coerciva. O mínimo satisfaz yᵢ(Σy)ᵢ = bᵢ, logo c = b, e w = y/1'y. A busca em linha é Armijo com y > 0 mantido por fração até a fronteira (a barreira com bᵢ < 1 não é autoconcordante padrão, então não se usa o passo amortecido de Newton).
  - Igualdades do Lagrangiano aumentado, nesta ordem e só as presentes: h₁ = 1'w − 1, h₂ = μ'w − r, h₃ = √(w'Σw) − s. Estacionaridade: ∇f(w) = Σₖ λₖ ∇hₖ + Σᵢ βᵢeᵢ, que é a convenção do header, com ∇h₃ = Σw/vol(w). O L-BFGS projetado minimiza L = f − λ'h + (ρ/2)‖h‖² no box, e o multiplicador é atualizado por λ ← λ − ρh.
- **Semântica do RISK_PARITY** (decisões deste plano, para conferir; as da spec estão em §3.3):
  - Sem caixa: `risk_free` nos dados é `INVALID_ARGUMENT`, `"risk_free does not apply to RISK_PARITY"`.
  - `leverage` é `INVALID_ARGUMENT`, `"leverage does not apply to RISK_PARITY"`.
  - `target_return` exige μ.
  - `risk_budget`: n entradas (conferido no solve, `"risk_budget has k entries; expected n"`), cada bᵢ > 0 finito, |Σbᵢ − 1| ≤ n·ε. A cota é a do erro de arredondamento da soma de n termos, não um chute. Padrão: 1/n.
  - `target_volatility`: s > 0 finito, só no RISK_PARITY.
  - `INFEASIBLE` só quando certificado: restrições lineares incompatíveis (a QP de mínima variância), `target_return` fora da faixa, ou s < σ_min.
    - s é comparado a σ_min com a constante 1e-12 relativa do MAX_RETURN, promovida a `problems::kMinimumVolatilityTolerance` (hoje `kMinimumTolerance` em `max_return.cpp`).
    - Com |s − σ_min| dentro dela, o único ponto viável é a carteira de mínima variância: `OPTIMAL`, multiplicadores de mínima norma.
    - Um s acima da volatilidade máxima de um conjunto viável limitado não é certificável de forma barata: o Lagrangiano aumentado estaciona com ‖h‖ > 0 e o status é `MAX_ITER`.
  - **`MAX_ITER` sem ponto, para todo objetivo:** pesos, multiplicadores e resíduo são NaN. Isso fecha a pendência do header ("documentado quando o RISK_PARITY for acrescentado") e uniformiza a regra.
  - `iterations` = iterações de Newton do Spinu + soma das internas do L-BFGS.
  - Multiplicadores: o de orçamento é finito (sem caixa); `target_return` e volatilidade são NaN quando o setter não foi chamado; leverage é NaN sempre.
  - `po_risk_contributions`: `out` tem n entradas (n = o número de pesos) e a do caixa é 0. Soma 1. É NaN quando a volatilidade é 0.
- **Matriz objetivo × setter**, acréscimo:

  | Objetivo | Exige | Aceita além de bounds |
  |---|---|---|
  | RISK_PARITY | covariance (mu com target_return) | target_return, target_volatility, risk_budget |

## Review Focus

1. **Conjunto viável de um ponto:** n = 2 com `target_return` e μ distintos fixa w; s = σ_min fixa a carteira de mínima variância. Os dois são `OPTIMAL`, sem laço infinito no Lagrangiano aumentado. Testes: `RiskParity.TwoAssetsWithTargetReturnIsThePoint` e `RiskParity.VolatilityEqualMinimumIsMinimumVariance` (Task 4).
2. **Médias iguais e `target_return` redundante ou inconsistente.** O redundante é `OPTIMAL` e igual à solução de Spinu. O inconsistente é `INFEASIBLE`. Testes: `RiskParity.EqualMeansRedundantTarget` e `RiskParity.EqualMeansInconsistentTarget` (Task 4).
3. **Ortante.** Com `lower < 0` e a solução de Spinu dentro do box, o resultado é a de Spinu (pesos positivos), não outra solução por ortante. Teste: `RiskParity.NegativeLowerBoundKeepsPositiveOrthant` (Task 2).
4. **Orçamentos muito desiguais (bᵢ de 1e-6 a 1) com cond(Σ) até 1e8.** Testes: `Spinu.UnevenBudgetsConditioningSweep` (Task 2) e a bateria `RiskParityFuzz` (Task 5).
5. **Não convergência nunca vira `OPTIMAL` fora de KKT.** s inalcançável num conjunto limitado dá `MAX_ITER` com NaN, e a bateria checa KKT e viabilidade em todo `OPTIMAL`. Testes: `RiskParity.UnreachableVolatilityIsMaxIter` (Task 4) e `RiskParityFuzz.NeverSilentlyWrong` (Task 5).

---

### Task 0: Referências independentes e KKT do usuário para RISK_PARITY

**Files:**
- Create: `tests/unit/fuzz/risk_parity_reference.hpp`, `tests/unit/fuzz/risk_parity_reference_test.cpp`
- Modify: `src/core/spec.hpp` (os dois campos novos, sem uso no solver ainda), `tests/unit/user_kkt.hpp`, `CMakeLists.txt` (o teste novo)

**Interfaces:**
- Consumes: `random_covariance`, `random_vector`, `spread_mu` (planos 1–3), `kkt.hpp`.
- Produces:
  ```cpp
  // core/spec.hpp
  std::vector<double> risk_budget;               // empty means 1/n
  std::optional<double> target_volatility;

  ```
  ```cpp
  namespace po::testing {
  // f(w) = sum (c_i - b_i)^2, plain loops in long double: independent of po::metrics.
  double risk_parity_value(const po::Data& data, std::span<const double> budget, std::span<const double> w);
  // grad f by the closed form of the plan header, coded here (checked against central differences in the self-test).
  std::vector<double> risk_parity_gradient(const po::Data& data, std::span<const double> budget, std::span<const double> w);
  // Positive-orthant solution of two assets (the root in (0,1) of c1 b2 = c2 b1 with w1 + w2 = 1).
  std::vector<double> two_asset_risk_parity(double sigma1, double sigma2, double rho, double b1);
  struct ScanResult { double value; std::vector<double> weights; };   // value NaN when no scanned point is feasible
  // Global minimum by scan, n = 3 only. Feasible set a curve: budget + target_return (a line, parametrized by one
  // weight) or budget + target_volatility (an ellipse around the minimum-variance point, by angle). 20001 points,
  // then golden-section on the best bracket; points outside the box are discarded.
  ScanResult scan_reference(const po::Spec& spec, const po::Data& data, std::span<const double> budget);
  struct RiskParityShape
  {
      std::size_t n;
      double condition;
      double scale;                     // 1 (annual) or 1e-4 (daily) on Sigma; sqrt on mu
      enum class Budget { Equal, Uneven } budget;       // Uneven: log-uniform b_i in [1e-6, 1], normalized
      enum class Box { None, LongOnly, Capped } box;    // Capped: lower 0, upper 1.2 / n
      bool target_return;               // r = a convex combination of two means
      bool target_volatility;           // s between sigma_min and 1.5 sigma_min of the linear set
      std::uint64_t seed;
  };
  struct GeneratedRiskParity { po::Spec spec; po::Data data; };
  GeneratedRiskParity generate(const RiskParityShape& shape);
  std::vector<RiskParityShape> risk_parity_shapes(std::size_t count, std::uint64_t seed);   // every combination
  // KKT in the user's variables: equalities, bounds, sign and complementarity of beta, stationarity with
  // risk_parity_gradient and the volatility term lambda_vol Sigma w / vol, NaN pattern of the multipliers, and
  // the identity w'grad f = 0.
  void expect_risk_parity_kkt(const po::Spec&, const po::Data&, const po::Result&, double tolerance);
  }
  ```
  `generate` preenche `risk_budget` e `target_volatility` no `Spec`. O teste de autoverificação do `scan_reference` com `target_volatility` usa um ponto de partida conhecido: Σ diagonal, s = vol da solução de Spinu e o valor da varredura é 0 em 1e-9.

- [ ] **Step 1: Testes** (de si mesmos, sem solver)
  - `RiskParityReference.GradientMatchesCentralDifferences`: n = 5 e 12, sementes 1..20, passo 1e-6, relativo 1e-6.
  - `RiskParityReference.DiagonalClosedForm`: Σ diagonal, wᵢ ∝ √bᵢ/σᵢ, valor 0 em 1e-14.
  - `RiskParityReference.TwoAssetsClosedFormHasZeroValue`: ρ ∈ {−0.5, 0, 0.9}, b₁ ∈ {0.5, 0.2, 0.01}.
  - `RiskParityReference.ScanFindsKnownMinimum`: n = 3, `target_return` igual ao retorno da solução diagonal de Spinu, o valor da varredura é 0 em 1e-9.
  - `RiskParityReference.ContributionsSumToOne`.
- [ ] **Step 2:** Rodar, ver falhar (cabeçalho inexistente). **Step 3:** Implementar os helpers. **Step 4:** PASS em `gcc-release`.
- [ ] **Step 5: Commit** — `test: independent risk parity references and user-variable KKT`

---

### Task 1: Métricas `risk_contributions` e resíduo (núcleo)

**Files:**
- Modify: `src/metrics/portfolio.hpp`, `src/metrics/portfolio.cpp`
- Test: `tests/unit/metrics/portfolio_test.cpp`

**Interfaces:**
- Produces:
  ```cpp
  namespace po::metrics {
  // c_i = w_i (Sigma w)_i / w'Sigma w over the risky weights; the cash entry (last, when present) is 0.
  // Sums to 1; all NaN when the volatility is 0. Needs the covariance.
  std::vector<double> risk_contributions(const Data& data, std::span<const double> weights);
  // sum_i (c_i - b_i)^2 over the risky weights; NaN when the volatility is 0. budget has asset_count entries.
  double risk_budget_residual(const Data& data, std::span<const double> weights, std::span<const double> budget);
  }
  ```

- [ ] **Step 1: Testes**
  - `Metrics.RiskContributionsByHand`: n = 2, Σ = [[0.04, 0.006], [0.006, 0.09]], w = (0.6, 0.4), c calculado à mão em 1e-15.
  - `Metrics.RiskContributionsSumToOneAndScaleInvariant`: c(tw) = c(w).
  - `Metrics.RiskContributionsCashEntryIsZero`: com `risk_free`, n + 1 entradas, a última é 0.
  - `Metrics.RiskContributionsZeroVolatilityIsNaN`: todo o peso no caixa.
  - `Metrics.RiskBudgetResidualZeroAtEqualContribution`: o ponto de `two_asset_risk_parity` dá resíduo < 1e-28.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar em `portfolio.cpp`. **Step 4:** PASS.
- [ ] **Step 5: Commit** — `feat(metrics): risk contributions and risk budget residual`

---

### Task 2: Spec, requisitos, pasta `solve`, Newton de Spinu e o caminho rápido

**Files:**
- Move: `src/problems/solve.*` e `src/problems/requirements.*` para `src/solve/`; `tests/unit/problems/{solve,requirements}_test.cpp` para `tests/unit/solve/` (`git mv`; incluir os caminhos novos em `CMakeLists.txt`, `ARCHITECTURE.md` e `tests/style`, se este enumerar pastas)
- Modify: `src/solve/requirements.cpp`, `src/solve/solve.cpp`, `CMakeLists.txt`
- Create: `src/riskparity/spinu.hpp`, `src/riskparity/spinu.cpp`, `src/riskparity/risk_parity.hpp`, `src/riskparity/risk_parity.cpp`
- Test: `tests/unit/riskparity/spinu_test.cpp`, `tests/unit/riskparity/risk_parity_test.cpp`, `tests/unit/solve/requirements_test.cpp`

**Interfaces:**
- Consumes: `metrics::risk_contributions` e `risk_budget_residual` (Task 1), `linalg::cholesky`, `linalg::solve_lower`.
- Produces:
  ```cpp
  namespace po::riskparity {
  struct SpinuResult { std::vector<double> weights; int iterations; bool converged; };
  // Positive-orthant solution w = y / 1'y of min 1/2 y'Sigma y - b' log y, by damped Newton with Armijo.
  SpinuResult spinu(const Data& data, std::span<const double> budget);
  std::vector<double> risk_budget_of(const Spec& spec, const Data& data);   // spec.risk_budget or 1/n
  std::expected<Result, Failure> solve_risk_parity(const Spec& spec, const Data& data);
  }
  ```
  `solve_risk_parity` nesta task cobre só o caminho rápido. Para `target_return`, `target_volatility` ou solução de Spinu fora do box ela devolve `Failure{FailureKind::Unsupported, "RISK_PARITY with active constraints not implemented yet"}`, e a Task 4 remove isso. `solve_request` despacha `Objective::RiskParity` para ela.

- [ ] **Step 1: Testes**
  - `Requirements.RiskParityMatrix`:
    - aceita `target_return` (exige mu), `target_volatility` e `risk_budget`;
    - `leverage` dá `"leverage does not apply to RISK_PARITY"`;
    - `max_volatility`, `risk_aversion` e `importance` dão `"<setter> does not apply to RISK_PARITY"`;
    - `risk_free` dá `"risk_free does not apply to RISK_PARITY"`;
    - `risk_budget` de tamanho errado dá `"risk_budget has 2 entries; expected 3"`.
  - `Requirements.TargetVolatilityOnlyForRiskParity`: nos outros objetivos, `"target_volatility does not apply to <OBJECTIVE>"`.
  - `Requirements.UnimplementedObjectivesStayUnsupported` sai, porque nenhum objetivo fica sem solver, e entra `Requirements.EveryObjectiveIsSolved`.
  - `Spinu.DiagonalClosedForm`: wᵢ ∝ √bᵢ/σᵢ em 1e-12, orçamento igual e desigual.
  - `Spinu.TwoAssetsClosedForm` contra `two_asset_risk_parity` em 1e-12.
  - `Spinu.ContributionsEqualBudget`: `risk_budget_residual` ≤ piso do arredondamento, para n ∈ {1, 2, 10, 100}, sementes 1..20, cond ∈ {1e2, 1e4, 1e6}.
  - `Spinu.UnevenBudgetsConditioningSweep`: bᵢ de 1e-6 a 1, cond ∈ {1e2, 1e4, 1e6, 1e8}, escala diária e anual; converge, w > 0 e as contribuições batem com b na faixa do cond.
  - `RiskParity.BudgetOnlyIsSpinu`: `OPTIMAL`, multiplicador de orçamento 0, `bound_multipliers` 0, `target_return_multiplier` e `volatility_multiplier` NaN, `leverage_multiplier` NaN, `expect_risk_parity_kkt` (Task 0).
  - `RiskParity.NegativeLowerBoundKeepsPositiveOrthant`: `bounds(-0.5, 1)` devolve o mesmo vetor de `bounds` aberto.
  - `RiskParity.SingleAssetIsTheWholeBudget`: n = 1, w = (1).
  - `RiskParity.ActiveConstraintsUnsupportedForNow`: `target_return` dá `Unsupported`.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar:
  - O move.
  - A linha de `requirements.cpp` (campos novos em `Requirements` e `Setter`, a checagem de caixa e a de tamanho do `risk_budget`).
  - `spinu`. Critério de parada: decremento de Newton λ² ≤ 16·n·ε·(1 + |F|), ou F deixa de decrescer; o fator é medido no Step 4. A Hessiana é fatorada por Cholesky e o passo vem de dois solves triangulares, sem inversa.
  - `solve_risk_parity` e o despacho.
- [ ] **Step 4:** PASS em `gcc-release`. Registrar em `validation.md` (seção nova "Numerical decisions of RISK_PARITY") o número de iterações de Newton por n e a medição que fixa o fator do critério.
- [ ] **Step 5: Commit** — `feat: RISK_PARITY by Spinu Newton; solve folder; risk_budget and target_volatility in the spec`

---

### Task 3: Objetivo, gradiente, igualdades e L-BFGS com bounds

**Files:**
- Create: `src/riskparity/objective.hpp`, `src/riskparity/objective.cpp`, `src/riskparity/box_lbfgs.hpp`, `src/riskparity/box_lbfgs.cpp`
- Test: `tests/unit/riskparity/objective_test.cpp`, `tests/unit/riskparity/box_lbfgs_test.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: `Data`, `Spec`, `risk_budget_of` (Task 2).
- Produces:
  ```cpp
  namespace po::riskparity {
  // f(w) and, when gradient is non-empty, grad f(w) (formula of the plan header). w'Sigma w must be positive.
  double risk_parity_objective(const Data& data, std::span<const double> budget, std::span<const double> w,
                               std::span<double> gradient);

  // The equalities of the spec, in the order budget, target_return, volatility (those present).
  class Equalities
  {
  public:
      static Equalities of(const Spec& spec, const Data& data);
      std::size_t count() const;
      // h (count) and jacobian (count x n, row-major) at w.
      void evaluate(std::span<const double> w, std::span<double> h, std::span<double> jacobian) const;
      // True for the volatility row; used to map multipliers to the named constraints.
      enum class Kind { Budget, TargetReturn, Volatility };
      std::span<const Kind> kinds() const;
  };

  using SmoothFunction = std::function<double(std::span<const double> x, std::span<double> gradient)>;
  struct BoxOptions { std::size_t memory; int max_iterations; double tolerance; };
  struct BoxResult { std::vector<double> x; std::vector<double> gradient; double value; int iterations; bool converged; };
  // L-BFGS on the free variables with projection onto [lower, upper] and a projected Armijo line search.
  // converged: ||x - P(x - gradient)||_inf <= tolerance.
  BoxResult minimize_box(const SmoothFunction& function, std::span<const double> lower, std::span<const double> upper,
                         std::vector<double> start, const BoxOptions& options);
  }
  ```

- [ ] **Step 1: Testes**
  - `RiskParityObjective.GradientMatchesIndependentReference`: contra `risk_parity_gradient` (Task 0) em 1e-12 relativo, n ∈ {2, 7, 30}.
  - `RiskParityObjective.ValueMatchesMetric`: igual a `risk_budget_residual` (Task 1).
  - `Equalities.JacobianByFiniteDifferences`: as três linhas, passo 1e-6, relativo 1e-6.
  - `BoxLbfgs.AgainstQpSolver`: QPs convexos com bounds (G = Σ aleatória, g aleatório, `lower`/`upper` apertados, n ∈ {5, 30}, sementes 1..30), `converged`, e x concorda com `qp::solve` em 1e-8 relativo.
  - `BoxLbfgs.ReturnsExactBoundValues`: nos pesos ativos, x é exatamente o bound (a projeção, não um arredondamento).
  - `BoxLbfgs.StartOutsideBoxIsProjected`.
  - `BoxLbfgs.IterationCapIsReported`: com `max_iterations = 1`, `converged = false`.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** PASS. Registrar em `validation.md` a concordância com o QP.
- [ ] **Step 5: Commit** — `feat: risk parity objective, equalities and projected L-BFGS`

---

### Task 4: Lagrangiano aumentado, viabilidade, multiplicadores e status

**Files:**
- Create: `src/riskparity/augmented_lagrangian.hpp`, `src/riskparity/augmented_lagrangian.cpp`
- Modify: `src/riskparity/risk_parity.cpp`, `src/problems/max_return.cpp` e `src/problems/constraints.hpp` (a constante compartilhada), `CMakeLists.txt`
- Test: `tests/unit/riskparity/augmented_lagrangian_test.cpp`, `tests/unit/riskparity/risk_parity_test.cpp`

**Interfaces:**
- Consumes: `spinu`, `risk_parity_objective`, `Equalities`, `minimize_box` (Tasks 2–3), `problems::solve_mean_variance`, `linalg::qr`.
- Produces:
  ```cpp
  namespace problems { inline constexpr double kMinimumVolatilityTolerance = 1e-12; }   // was kMinimumTolerance in max_return.cpp

  namespace po::riskparity {
  struct AugmentedResult
  {
      std::vector<double> weights;
      std::vector<double> multipliers;          // one per Equalities row, convention of the header
      int outer_iterations;
      int inner_iterations;
      bool converged;                           // KKT residual and ||h|| under tolerance after restoration
  };
  AugmentedResult augmented_lagrangian(const Spec& spec, const Data& data, std::span<const double> budget,
                                       std::span<const double> start);
  }
  ```
  Política em `solve_risk_parity` (a ordem é o contrato):
  1. Spinu; sem convergência → `MAX_ITER`.
  2. Caminho rápido (Task 2).
  3. A QP de mínima variância com os bounds e o `target_return` do spec (`problems::solve_mean_variance` sobre uma cópia com `objective = MinVariance`): `Infeasible` → `INFEASIBLE`. Com `target_volatility`: s < σ_min·(1 − kMinimumVolatilityTolerance) → `INFEASIBLE`, e |s − σ_min| ≤ kMinimumVolatilityTolerance·s → esse ponto.
  4. `augmented_lagrangian` partindo da solução de Spinu. As iterações externas aumentam ρ ×10 quando ‖h‖ não cai por um fator 4 (Nocedal-Wright §17.4). O subproblema é `minimize_box`. No fim, restauração por Gauss-Newton em h sobre as variáveis livres, e os multiplicadores (igualdades e β dos bounds ativos, por QR) vêm do gradiente.
  5. Verificação final de KKT e viabilidade, relativa ao ruído do problema. Se falhar, o status é `MAX_ITER`.
  6. Um β de sinal errado tira o bound do conjunto ativo e continua o Lagrangiano aumentado.
  7. `MAX_ITER` devolve pesos, multiplicadores e resíduo NaN.

- [ ] **Step 1: Testes** (todo `OPTIMAL` chama `expect_risk_parity_kkt`)
  - `RiskParity.TwoAssetsWithTargetReturnIsThePoint`: n = 2, μ distintos, w vem do sistema 2×2 (orçamento e retorno), resíduo igual a `risk_parity_value`.
  - `RiskParity.TargetReturnLineAgainstScan`: n = 3, `target_return` entre os μ, resíduo ≤ varredura·(1 + 1e-9) e pesos a 1e-6 quando o mínimo é único.
  - `RiskParity.TargetVolatilityEllipseAgainstScan`: n = 3, s entre σ_min e 1.5σ_min.
  - `RiskParity.VolatilityEqualMinimumIsMinimumVariance`: s = σ_min, w igual ao de `MIN_VARIANCE`, e o `expect_risk_parity_kkt` aceita multiplicadores de mínima norma.
  - `RiskParity.VolatilityBelowMinimumInfeasible`, `RiskParity.TargetReturnOutOfRangeInfeasible` e `RiskParity.BoundsInfeasible` (lower·n > 1): `INFEASIBLE` com NaN.
  - `RiskParity.UnreachableVolatilityIsMaxIter`: Σ = diag(0.04, 0.09), long-only, s = 0.5 → `MAX_ITER`, NaN.
  - `RiskParity.EqualMeansRedundantTarget` (target = média comum → `OPTIMAL`, igual a Spinu) e `RiskParity.EqualMeansInconsistentTarget` (→ `INFEASIBLE`).
  - `RiskParity.ActiveBoundsKkt`: box `Capped` com um peso no teto, sinais de β e complementaridade.
  - `RiskParity.HomogeneityIdentity`: λ_orçamento = −λ_retorno·r − λ_vol·s − Σβᵢwᵢ em 1e-9 relativo.
  - `RiskParity.MultipliersAreDerivatives`: diferença central do f ótimo em r e em s (passo 1e-6 relativo, tolerância relativa 1e-5), em casos com mínimo local único.
  - `RiskParity.ScaleInvariance`: Σ → tΣ com s → √t·s e os mesmos pesos.
  - `AugmentedLagrangian.AgreesWithSpinuWhenBothApply`: só o orçamento como igualdade, partida perturbada, pesos concordam em 1e-8 (spec §9).
  - `AugmentedLagrangian.IterationCapsAreReported`.
  - `RiskParity.ActiveConstraintsUnsupportedForNow` (Task 2) sai.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** PASS em `gcc-release`.
- [ ] **Step 5: Commit** — `feat: RISK_PARITY with target return, target volatility and active bounds`

---

### Task 5: Bateria de fuzz, medição das tolerâncias e tetos de iteração

**Files:**
- Create: `tests/unit/fuzz/risk_parity_fuzz_test.cpp`
- Modify: `src/riskparity/*.cpp` (os fatores medidos), `CMakeLists.txt`, `docs/guide/validation.md`

**Interfaces:**
- Consumes: `generate`, `risk_parity_shapes`, `scan_reference`, `expect_risk_parity_kkt` (Task 0), `po::solve`.

- [ ] **Step 1: Testes**
  - `RiskParityFuzz.KktAndFeasibility`:
    - `risk_parity_shapes(400, 1)`;
    - todo `OPTIMAL` tem pesos somando 1, dentro dos bounds, `target_return` e vol(w) = s atendidos, KKT na faixa do cond, e resíduo igual a `risk_parity_value`;
    - todo `INFEASIBLE` é confirmado pela QP de mínima variância ou por s < σ_min.
  - `RiskParityFuzz.AgainstScanReference`:
    - n = 3 com `target_return` ou `target_volatility`, 200 casos;
    - o resíduo de PortOpt ≤ o da varredura·(1 + 1e-9);
    - um caso que falha só é aceito se o ponto da varredura for outro ponto de KKT, o que prova mínimo local;
    - a contagem desses casos entra em `validation.md`.
  - `RiskParityFuzz.NeverSilentlyWrong`:
    - qualquer `MAX_ITER` tem NaN em pesos, multiplicadores e resíduo;
    - nunca há `OPTIMAL` fora de KKT ou de viabilidade;
    - a taxa de `MAX_ITER` é medida e registrada, por forma.
  - `RiskParityFuzz.NoSampledFeasiblePointBeatsOptimum`: 10.000 pontos viáveis amostrados (semente fixa) em `None`, `LongOnly` e `Capped` sem restrições extras, nenhum com f menor que o ótimo·(1 − 1e-9).
  - `RiskParityFuzz.UnevenBudgetsAndConditioning` (Review Focus 4): bᵢ de 1e-6 a 1, cond(Σ) ∈ {1e2, 1e4, 1e6, 1e8}, escala diária e anual, orçamento só e com bounds `Capped`.
- [ ] **Step 2:** Rodar. Falhas viram análise, não afrouxamento de faixa.
- [ ] **Step 3:** Medir e fixar, com cada valor em `validation.md` junto da medição:
  - o fator de parada do Spinu;
  - o fator de KKT do Lagrangiano aumentado;
  - os tetos `kMaxOuter` e `kMaxInner` (10× o maior valor observado nos casos que convergem; partida 50 e 1000).
- [ ] **Step 4:** PASS em `gcc-release`.
- [ ] **Step 5: Commit** — `test: risk parity fuzz battery; measured tolerances`

---

### Task 6: Interface C e métricas

**Files:**
- Modify: `include/portopt/portopt.h`, `src/capi/spec_api.cpp`, `src/capi/metrics_api.cpp`, `src/capi/result_api.cpp` (só se o resíduo pedir), `tests/unit/capi/c_header_check.c`
- Test: `tests/unit/capi/spec_test.cpp`, `tests/unit/capi/solve_test.cpp`, `tests/unit/capi/metrics_test.cpp`

**Interfaces:**
- Produces:
  ```c
  /* Volatility target sigma of PO_OBJECTIVE_RISK_PARITY: sqrt(w' Sigma w) = sigma. Positive and finite. */
  PO_API po_status po_spec_set_target_volatility(po_spec* s, double sigma);
  /* Risk budget b of PO_OBJECTIVE_RISK_PARITY: n entries, each positive and finite, summing to 1 within
     n * DBL_EPSILON. n must be the asset count of the data given to po_solve (checked by po_solve). */
  PO_API po_status po_spec_set_risk_budget(po_spec* s, const double* b, size_t n);
  /* Risk contributions c_i = w_i (Sigma w)_i / w' Sigma w; out has n entries (n = the weight count), the cash
     entry is 0, the entries sum to 1; all NaN when the volatility is 0. Needs the covariance. */
  PO_API po_status po_risk_contributions(const po_data* d, const double* w, size_t n, double* out);
  ```
  O header documenta:
  - que todo objetivo é resolvido (a tabela ganha `RISK_PARITY`) e a linha da matriz;
  - a forma f = Σ(cᵢ − bᵢ)², o ponto de KKT local e o resíduo em `po_result_residual`;
  - `PO_CONSTRAINT_VOLATILITY` para RISK_PARITY com `target_volatility` (a derivada de f* em s; NaN sem o setter);
  - a regra nova de `MAX_ITER` (sem ponto, em todo objetivo) e o `INFEASIBLE` certificado;
  - o termo de volatilidade da fórmula de estacionaridade (hoje só do MAX_RETURN) e a ausência de caixa e de leverage.
  `c_header_check.c` ganha as três funções, e o valor de `PO_ABI_VERSION` continua 0.

- [ ] **Step 1: Testes**
  - `CapiSpec.RiskBudgetDomain`: n = 0, NaN, ≤ 0, soma 1.1 e soma 1 + n·ε·4 são `INVALID_ARGUMENT` e deixam a spec inalterada, e a soma 1 ± (n·ε) é aceita.
  - `CapiSpec.TargetVolatilityDomain`: 0, negativo, NaN e ±inf.
  - `CapiSolve.RiskParityEndToEnd`: bit a bit igual a `po::solve`, em quatro casos (orçamento só, com `target_return`, com `target_volatility`, com bounds `Capped`).
  - `CapiSolve.RiskParityRejectsCashAndLeverage`.
  - `CapiSolve.RiskParityResidualAndMultipliers`: `po_result_residual` bate com `risk_budget_residual`, multiplicadores NaN onde o setter não foi chamado.
  - `CapiSolve.MaxIterHasNoPoint` (`UnreachableVolatility`) e `CapiSolve.MaxIterHasNoPointForQpObjectives` (a regra uniforme).
  - `Metrics.RiskContributions` (à mão, com caixa, volatilidade 0) e `Metrics.RiskContributionsSizeMismatch`.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** `ctest` PASS em `gcc-release` (inclui `c_header_check` e `style`).
- [ ] **Step 5: Commit** — `feat(capi): RISK_PARITY, risk budget, target volatility, risk contributions`

---

### Task 7: Python e documentos

**Files:**
- Modify: `python/src/portopt/_native.py`, `_spec.py`, `_data.py`, `_result.py`
- Test: `python/tests/test_risk_parity.py`
- Modify: `docs/guide/math.md`, `docs/guide/validation.md`, `ARCHITECTURE.md`, `README.md`, `CHANGELOG.md`

**Interfaces:**
- Produces:
  - `Spec.target_volatility: float | None = None`.
  - `Spec.risk_budget: ArrayLike | None = None`, canonicalizado para tupla para que a igualdade e o hash continuem por valor (como os bounds).
  - `Data.risk_contributions(weights) -> NDArray[np.float64]`.
  - `Result.residual: float`.
  - `SCALAR_SETTERS` ganha `("target_volatility", "po_spec_set_target_volatility")`; `risk_budget` usa `po_spec_set_risk_budget`.

- [ ] **Step 1: Testes** (`test_risk_parity.py`)
  - `test_equal_risk_contribution_diagonal` (numpy, forma fechada);
  - `test_risk_budget_matches_contributions`;
  - `test_residual_zero_unconstrained`;
  - `test_target_return_residual_positive_and_met`;
  - `test_target_volatility_met`;
  - `test_risk_contributions_sum_to_one_and_match_numpy`;
  - `test_cash_and_leverage_raise_value_error`;
  - `test_unreachable_volatility_is_max_iter_with_nan`;
  - `test_spec_with_risk_budget_array_hashes`.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar.
- [ ] **Step 4:** Documentação:
  - `math.md`: seção RISK_PARITY com c, f, ∇f (derivação de J), Spinu, o Lagrangiano aumentado, a restauração e os multiplicadores, a identidade de homogeneidade, o caráter local e o `MAX_ITER`; a tabela "Which setters and data each objective takes" ganha a linha nova.
  - `validation.md`: a seção "Numerical decisions of RISK_PARITY" (medições das Tasks 2, 3 e 5), as referências por varredura e as faixas por cond.
  - `ARCHITECTURE.md`: `riskparity`, `solve`, o diagrama de camadas, a contagem de módulos e `risk_parity_reference.hpp`.
  - `README.md`: um exemplo de RISK_PARITY com `target_volatility`.
  - `CHANGELOG.md`.
- [ ] **Step 5:** Verificar:
  - pytest PASS;
  - mypy strict (`--config-file python/pyproject.toml`) limpo;
  - `check_style.py` limpo;
  - `ctest` PASS em `gcc-release`.
- [ ] **Step 6: Commit** — `feat(python): RISK_PARITY, risk budget and target volatility; docs for plan 4`
