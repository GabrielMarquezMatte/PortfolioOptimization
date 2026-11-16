# PortOpt — Plano 4: RISK_PARITY, target_volatility e risk_budget

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Entregar `RISK_PARITY` de ponta a ponta (núcleo → interface C → Python): Newton de Spinu para o caso sem restrições extras, uma sequência de QPs de Gauss-Newton sobre o poliedro (orçamento, `target_return`, bounds) resolvidos pelo Goldfarb-Idnani existente, com Lagrangiano aumentado só para `target_volatility` e polimento Newton-KKT, os setters `risk_budget` e `target_volatility`, o resíduo e a métrica `risk_contributions`, validados por KKT, formas fechadas e referências independentes por varredura.

**Architecture:**
- **Pasta nova `src/riskparity/`** (`spinu`, `objective`, `polyhedral_descent`, `augmented_lagrangian`, `newton_polish`, `risk_parity`). `src/problems/` já tem 8 módulos (limite do STYLEGUIDE), então `solve.*` e `requirements.*` saem dela para `src/solve/`. Ordem das dependências: `core` ← `qp` ← `problems` ← `riskparity` ← `solve` ← `capi`.
- **Caminho rápido:** sem `target_return`, sem `target_volatility` e com a solução de Spinu dentro dos bounds, o resultado é a de Spinu, `OPTIMAL`, todos os multiplicadores 0 e resíduo no piso do arredondamento.
- **Caminho geral:** viabilidade e σ_min pela mínima variância (`problems::solve_mean_variance`), depois, partindo da solução de Spinu, passos de Gauss-Newton com termo proximal sobre o poliedro Ω = {orçamento, `target_return`, bounds}, cada passo um QP do `qp::solve` com busca de Armijo (as linhas lineares valem exatamente em todo iterado e os multiplicadores saem do QP). Só vol(w) = s entra num Lagrangiano aumentado (regra de ρ do ALGENCAN). Um polimento Newton-KKT com a Hessiana exata leva a estacionaridade ao piso do arredondamento. Revisão por pesquisa: `docs/research/2026-11-16-risk-parity-and-numerics.md`.
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
  - Spinu: F(y) = ½ y'Σy − b' log y, y > 0, ∇F = Σy − b/y e ∇²F = Σ + diag(b/y²), estritamente convexa e coerciva. O mínimo satisfaz yᵢ(Σy)ᵢ = bᵢ, logo c = b, e w = y/1'y. A busca em linha é Armijo com y > 0 mantido por fração até a fronteira (F não é autoconcordante padrão com bᵢ < 1; F/min bᵢ é, mas a cota de iterações escala com 1/min bᵢ). O critério de parada é testado **depois** de dar o passo completo, ou por ‖Δy/y‖∞ ≤ c·ε (a medição fixa o c), porque o decremento mede F − F*, não o erro do ponto (erro relativo do peso ~1e-7/√bᵢ se o iterado devolvido for o do teste).
  - Igualdades, nesta ordem e só as presentes: h₁ = 1'w − 1, h₂ = μ'w − r (ambas lineares, no poliedro Ω com os bounds) e h₃ = √(w'Σw) − s (no Lagrangiano aumentado L = f − λh₃ + (ρ/2)h₃², atualização λ ← λ − ρh₃, ρ ×10 quando |h₃| não cai por 0.5, λ salvaguardado, tolerância interna decrescente: ABMS 2007, Alg. 3.1). Estacionaridade: ∇f(w) = Σₖ λₖ ∇hₖ + Σᵢ βᵢeᵢ, que é a convenção do header, com ∇h₃ = Σw/vol(w).
  - Gauss-Newton: J tem Jw = 0 (c é homogênea de grau 0), então J'J é singular e o QP leva o termo proximal τI, τ = 0.05·tr(Σ)/(2n). A Hessiana exata de f é indefinida em geral (problema não convexo: Bai-Scheinberg-Tütüncü 2016, Ex. 3.1 e 3.2) e só entra no polimento final.
- **Semântica do RISK_PARITY** (decisões deste plano, para conferir; as da spec estão em §3.3):
  - Sem caixa: `risk_free` nos dados é `INVALID_ARGUMENT`, `"risk_free does not apply to RISK_PARITY"`.
  - `leverage` é `INVALID_ARGUMENT`, `"leverage does not apply to RISK_PARITY"`.
  - `target_return` exige μ.
  - `risk_budget`: n entradas (conferido no solve, `"risk_budget has k entries; expected n"`), cada bᵢ > 0 finito, |Σbᵢ − 1| ≤ n·ε. A cota é a do erro de arredondamento da soma de n termos, não um chute. Padrão: 1/n.
  - `target_volatility`: s > 0 finito, só no RISK_PARITY.
  - `INFEASIBLE` só quando certificado: restrições lineares incompatíveis (a QP de mínima variância), `target_return` fora da faixa, ou s < σ_min.
    - s é comparado a σ_min com a constante 1e-12 relativa do MAX_RETURN, promovida a `problems::kMinimumVolatilityTolerance` (hoje `kMinimumTolerance` em `max_return.cpp`).
    - A constante vira `qp::kRelativeTolerance` (a tolerância de viabilidade do QP, que é o que a sustenta; ver `docs/research/…` item 5).
    - Com |s − σ_min| dentro dela, o único ponto viável é a carteira de mínima variância: `OPTIMAL`, com **λ_vol = −∞ e os demais multiplicadores NaN**, como o `MAX_RETURN` na fronteira vertical (aí a LICQ falha, em geral não existem multiplicadores de KKT e ∂f*/∂s = −∞; "multiplicadores de mínima norma" não existem).
    - Um s acima da volatilidade máxima de um conjunto viável limitado não é certificável de forma barata: o Lagrangiano aumentado estaciona com |h₃| > 0 e o status é `MAX_ITER`.
  - **`MAX_ITER` sem ponto, para todo objetivo:** pesos, multiplicadores e resíduo são NaN. Isso fecha a pendência do header ("documentado quando o RISK_PARITY for acrescentado") e uniformiza a regra.
  - `iterations` = iterações de Newton do Spinu + soma das iterações de `descend` + passos do polimento.
  - Multiplicadores: o de orçamento é finito (sem caixa); `target_return` e volatilidade são NaN quando o setter não foi chamado; leverage é NaN sempre.
  - `po_risk_contributions`: `out` tem n entradas (n = o número de pesos) e a do caixa é 0. Soma 1. É NaN quando a volatilidade é 0.
- **Matriz objetivo × setter**, acréscimo:

  | Objetivo | Exige | Aceita além de bounds |
  |---|---|---|
  | RISK_PARITY | covariance (mu com target_return) | target_return, target_volatility, risk_budget |

## Review Focus

1. **Conjunto viável de um ponto, dois casos diferentes.** n = 2 com `target_return` e μ distintos fixa w, mas a LICQ vale e os multiplicadores existem. s = σ_min fixa a carteira de mínima variância, a LICQ falha e λ_vol = −∞ (resto NaN). Os dois são `OPTIMAL`, sem laço infinito no Lagrangiano aumentado. Testes: `RiskParity.TwoAssetsWithTargetReturnIsThePoint` e `RiskParity.VolatilityEqualMinimumIsMinimumVariance` (Task 4).
2. **Médias iguais e `target_return` redundante ou inconsistente.** O redundante é `OPTIMAL` e igual à solução de Spinu. O inconsistente é `INFEASIBLE`. Testes: `RiskParity.EqualMeansRedundantTarget` e `RiskParity.EqualMeansInconsistentTarget` (Task 4).
3. **Ortante.** A spec §3.3 ("uma solução por ortante") está errada em geral: há no máximo 2^(n−1) soluções normalizadas, uma por par de ortantes opostos, e nenhuma quando a soma dos pesos do raio é 0 (Bai-Scheinberg-Tütüncü 2016, Lema 2.2). A solução de pesos positivos existe e é única para Σ definida positiva (Bruder-Roncalli 2012, §2.3.3). Com `lower < 0` e a solução de Spinu dentro do box, o resultado é a de pesos positivos; um `OPTIMAL` com bounds pode perder uma solução exata de outro ortante. Teste: `RiskParity.NegativeLowerBoundKeepsPositiveOrthant` (Task 2).
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

### Task 3: Objetivo, Jacobiano, Hessiana e a sequência de QPs sobre o poliedro

Contexto da revisão (`docs/research/2026-11-16-risk-parity-and-numerics.md`, item 1 e 2): a Task 3 original (L-BFGS projetado, `box_lbfgs`) foi **substituída**. O subproblema é um mínimos quadrados não linear sobre o poliedro Ω = {1'w = 1, μ'w = r (se houver), lo ≤ w ≤ hi}, resolvido por passos de Gauss-Newton com termo proximal, cada um um QP convexo do `qp::solve` (Goldfarb-Idnani), como o `riskParityPortfolio` faz com `quadprog`. As igualdades lineares e os bounds valem exatamente em todo iterado, e não há módulo de L-BFGS para validar.

**Files:**
- Create: `src/riskparity/objective.hpp`, `src/riskparity/objective.cpp`, `src/riskparity/polyhedral_descent.hpp`, `src/riskparity/polyhedral_descent.cpp`
- Test: `tests/unit/riskparity/objective_test.cpp`, `tests/unit/riskparity/polyhedral_descent_test.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: `Data`, `Spec`, `risk_budget_of` (Task 2), `problems::mean_variance_problem` (constrói as linhas de orçamento, retorno e bounds de `qp::Problem`), `qp::solve`.
- Produces:
  ```cpp
  namespace po::riskparity {
  // r(w) = c(w) - b, f = |r|^2. Fills value, gradient (2 J'r), the Jacobian J (n x n, row-major) and, when
  // hessian is non-empty, the exact Hessian of f: 2 J'J + 2 sum_i r_i hess(c_i). w'Sigma w must be positive.
  struct Linearization { double value; std::vector<double> gradient; std::vector<double> jacobian; std::vector<double> hessian; };
  Linearization linearize(const Data& data, std::span<const double> budget, std::span<const double> w, bool with_hessian);

  // Minimizes f(w) + extra(w) over the polyhedron of the spec (budget, target_return, bounds) from a feasible start.
  // extra is the augmented-Lagrangian term of the volatility row: value, gradient and the rank-one curvature
  // rho * grad h3 grad h3' (Gauss-Newton), or null when there is no target_volatility.
  struct VolatilityTerm { double multiplier; double penalty; double target; };   // L3 = -multiplier*h3 + penalty/2*h3^2
  struct DescentOptions { double proximal; double tolerance; int max_iterations; };
  struct DescentResult
  {
      std::vector<double> weights;
      std::vector<double> equality_multipliers;      // of the QP rows at the fixed point: budget, target_return (those present)
      std::vector<double> bound_multipliers;         // beta, convention of the header
      double stationarity;                           // relative KKT residual at weights
      int iterations;
      bool converged;
  };
  DescentResult descend(const Spec& spec, const Data& data, std::span<const double> budget,
                        std::optional<VolatilityTerm> volatility, std::span<const double> start,
                        const DescentOptions& options);
  }
  ```
  Passo: com Q = 2J'J + τI (+ ρ ∇h₃∇h₃' quando há volatilidade), τ = `options.proximal`·tr(Σ)/(2n) (a escolha do pacote `riskParityPortfolio`, 0.05), resolver o QP min ½ŵ'Qŵ + (∇L − Qw)'ŵ sobre Ω com `qp::solve`, direção d = ŵ − w, Armijo em L ao longo de w + γd (sempre viável pela convexidade de Ω). No ponto fixo (d = 0), ∇L(w) = N'u, e os multiplicadores do QP são os da estacionaridade original (com λ₃ = multiplier − penalty·h₃). Parada por KKT relativo (não pela variação de w ou de f). O QP infeasible não ocorre porque a partida é viável (Task 4 só chama depois da checagem de mínima variância).

- [ ] **Step 1: Testes**
  - `RiskParityObjective.GradientMatchesIndependentReference`: contra `risk_parity_gradient` (`fuzz/risk_parity/risk_parity_reference.hpp`) em 1e-12 relativo, n ∈ {2, 7, 30}.
  - `RiskParityObjective.ValueMatchesMetric`: igual a `metrics::risk_budget_residual`.
  - `RiskParityObjective.JacobianByCentralDifferences` e `RiskParityObjective.HessianByCentralDifferencesOfGradient`: passo 1e-6, relativo 1e-6 (n ∈ {3, 10}); `Jw = 0` (homogeneidade) em 1e-12.
  - `PolyhedralDescent.AtSpinuPointStaysPut`: partindo do ponto de Spinu (f = 0) com box folgado, devolve o mesmo ponto, `converged`, multiplicadores 0.
  - `PolyhedralDescent.TargetReturnLineAgainstScan`: n = 3, `target_return`, partida viável arbitrária, resíduo ≤ varredura·(1 + 1e-9) (ou outro KKT, caso mínimo local) e `expect_risk_parity_kkt`.
  - `PolyhedralDescent.ActiveBoundsKkt`: box `Capped` com um peso no teto, sinais e complementaridade de β.
  - `PolyhedralDescent.MonotoneDescent`: f + extra nunca cresce de uma iteração para a seguinte (`iterations` e o valor final ≤ o inicial).
  - `PolyhedralDescent.IterationCapIsReported`: `max_iterations = 1` dá `converged = false`.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** PASS em `gcc-release`. Registrar em `validation.md` o número de iterações por n e o valor de `proximal` usado.
- [ ] **Step 5: Commit** — `feat: risk parity objective derivatives and polyhedral Gauss-Newton descent`

---

### Task 4: Volatilidade-alvo por Lagrangiano aumentado, polimento, viabilidade e status

Contexto da revisão (`docs/research/…`, itens 1 e 5). O Lagrangiano aumentado entra **só** para h₃ = √(w'Σw) − s (as demais igualdades e os bounds estão em Ω). Sem `target_volatility` não há laço externo: é `descend` puro.

**Files:**
- Create: `src/riskparity/augmented_lagrangian.hpp`, `src/riskparity/augmented_lagrangian.cpp`, `src/riskparity/newton_polish.hpp`, `src/riskparity/newton_polish.cpp`
- Modify: `src/riskparity/risk_parity.cpp`, `src/problems/max_return.cpp` e `src/problems/constraints.hpp` (a constante compartilhada), `tests/unit/user_kkt.hpp` ou o `expect_risk_parity_kkt` (o caso s = σ_min), `CMakeLists.txt`
- Test: `tests/unit/riskparity/augmented_lagrangian_test.cpp`, `tests/unit/riskparity/newton_polish_test.cpp`, `tests/unit/riskparity/risk_parity_test.cpp`

**Interfaces:**
- Consumes: `spinu`, `linearize`, `descend` (Tasks 2–3), `problems::solve_mean_variance`, `linalg::qr`, `linalg::cholesky`.
- Produces:
  ```cpp
  namespace po::problems {
  // Relative tolerance of the QP feasibility (qp kRelativeTolerance, 1e-12), shared with MAX_RETURN; replaces
  // kMinimumTolerance in max_return.cpp.
  inline constexpr double kMinimumVolatilityTolerance = qp::kRelativeTolerance;
  }
  namespace po::riskparity {
  struct AugmentedResult
  {
      std::vector<double> weights;
      std::vector<double> equality_multipliers;   // budget, target_return, volatility (those present), header convention
      std::vector<double> bound_multipliers;
      int outer_iterations;
      int inner_iterations;
      bool converged;
  };
  // ABMS 2007 Algorithm 3.1: rho grows by 10 when |h3| does not fall by a factor 0.5, lambda is safeguarded in
  // [-lambda_max, lambda_max], the inner tolerance decreases (epsilon_{k+1} = max(epsilon_final, 0.1 epsilon_k)).
  AugmentedResult augmented_lagrangian(const Spec& spec, const Data& data, std::span<const double> budget,
                                       std::span<const double> start);
  // Newton on the KKT system with the exact Hessian of the Lagrangian, active rows as equalities; accepts the step
  // only when the feasibility and the KKT residual decrease. Returns the polished point and multipliers.
  struct PolishResult { std::vector<double> weights; std::vector<double> multipliers; bool improved; };
  PolishResult newton_polish(const Spec& spec, const Data& data, std::span<const double> budget,
                             std::span<const double> weights, std::span<const double> active_rows_info);
  }
  ```
  (a assinatura exata de `newton_polish` é do implementador; ela recebe o conjunto ativo do último QP.)

  Política em `solve_risk_parity` (a ordem é o contrato):
  1. Spinu; sem convergência → `MAX_ITER`.
  2. Caminho rápido (Task 2).
  3. A QP de mínima variância com os bounds e o `target_return` do spec (`problems::solve_mean_variance` sobre uma cópia com `objective = MinVariance`): `Infeasible` → `INFEASIBLE`. Com `target_volatility`:
     - s < σ_min·(1 − kMinimumVolatilityTolerance) → `INFEASIBLE`.
     - |s − σ_min| ≤ kMinimumVolatilityTolerance·s → o ponto de mínima variância, `OPTIMAL`, **λ_vol = −∞ e os demais multiplicadores NaN** (nesse ponto a LICQ falha e em geral não há multiplicadores de KKT, e ∂f*/∂s = −∞; é a mesma convenção do `MAX_RETURN` na fronteira vertical, `math.md`). Sem laço externo.
  4. Sem `target_volatility`: `descend` a partir de Spinu, projetado em Ω quando a partida de Spinu viola um bound (um QP de projeção sobre Ω). Com `target_volatility`: `augmented_lagrangian`.
  5. Polimento: `newton_polish` com o conjunto ativo do último QP; aceita só se a viabilidade e o resíduo KKT caem.
  6. Verificação final de KKT e viabilidade, relativa ao ruído do problema. Se falhar, `MAX_ITER`.
  7. Um β de sinal errado tira o bound do conjunto ativo e continua.
  8. `MAX_ITER` devolve pesos, multiplicadores e resíduo NaN.
  `iterations` = Newton do Spinu + iterações de `descend` (soma de todas as externas) + passos do polimento.

- [ ] **Step 1: Testes** (todo `OPTIMAL` com multiplicadores finitos chama `expect_risk_parity_kkt`; o caso s = σ_min tem o ramo próprio, abaixo)
  - `RiskParity.TwoAssetsWithTargetReturnIsThePoint`: n = 2, μ distintos. A LICQ vale e o ponto é isolado mas regular: w vem do sistema 2×2, os multiplicadores existem e são únicos, o resíduo é igual a `risk_parity_value`.
  - `RiskParity.TargetReturnLineAgainstScan`, `RiskParity.TargetVolatilityEllipseAgainstScan`: n = 3, na varredura de referência.
  - `RiskParity.VolatilityEqualMinimumIsMinimumVariance`: s = σ_min do `MIN_VARIANCE` → w igual ao do `MIN_VARIANCE`, `OPTIMAL`, λ_vol = −∞ e o resto NaN; s = σ_min·(1 + 1e-9) → resultado regular com multiplicadores finitos e λ_vol muito negativo.
  - `RiskParity.VolatilityBelowMinimumInfeasible`, `RiskParity.TargetReturnOutOfRangeInfeasible` e `RiskParity.BoundsInfeasible` (lower·n > 1): `INFEASIBLE` com NaN.
  - `RiskParity.UnreachableVolatilityIsMaxIter`: Σ = diag(0.04, 0.09), long-only, s = 0.5 → `MAX_ITER`, NaN.
  - `RiskParity.EqualMeansRedundantTarget` (target = média comum → `OPTIMAL`, igual a Spinu) e `RiskParity.EqualMeansInconsistentTarget` (→ `INFEASIBLE`).
  - `RiskParity.ActiveBoundsKkt`, `RiskParity.HomogeneityIdentity` (λ_orçamento = −λ_retorno·r − λ_vol·s − Σβᵢwᵢ em 1e-9 relativo), `RiskParity.ScaleInvariance` (Σ → tΣ, s → √t·s, mesmos pesos).
  - `RiskParity.MultipliersAreDerivatives`: diferença central do f ótimo em r e em s (passo 1e-6 relativo, tolerância 1e-5), em casos com mínimo local único.
  - `AugmentedLagrangian.SafeguardedMultiplierStaysBounded`: com s = σ_min·(1 + 1e-6) λ não explode e ρ cresce de forma finita.
  - `AugmentedLagrangian.IterationCapsAreReported`.
  - `NewtonPolish.ReachesRoundingFloor`: partindo de um ponto KKT a 1e-6, a estacionaridade cai para a faixa do cond em 1 a 2 passos; um passo que piora é recusado (`improved = false`).
  - `expect_risk_parity_kkt`: acrescentar os casos de rejeição apontados na revisão da Task 0 (λ_orçamento perturbado, β de sinal errado, β NaN num bound infinito) e uma variante com `target_return`.
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
  - o termo proximal e a tolerância de `descend` (Task 3) e o fator de KKT final;
  - os tetos de iteração externas e internas do `augmented_lagrangian` e de `descend` (10× o maior valor observado nos casos que convergem; partida 50 e 1000);
  - a aceitação do `newton_polish`;
  - **s ≈ σ_min por duas rotas:** medir em quantos casos o σ_min de `MIN_VARIANCE` difere em mais de `kMinimumVolatilityTolerance` do σ_min da QP interna (o `validation.md` do plano 3 registra 17 de 320 no `MAX_RETURN`). Um usuário que passa o s vindo de `MIN_VARIANCE` recebe `INFEASIBLE` nesses casos; registrar a taxa e, se for relevante, decidir com o autor.
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
  - `math.md`: seção RISK_PARITY com c, f, ∇f (derivação de J), Spinu (o critério de parada pós-passo e a autoconcordância de F/min bᵢ), o Gauss-Newton proximal sobre o poliedro, o Lagrangiano aumentado de vol(w) = s, o polimento, os multiplicadores, a identidade de homogeneidade, o caso s = σ_min (λ_vol = −∞), o caráter local e o `MAX_ITER`;
  - **corrigir a spec** §3.3 (a frase "uma solução por ortante") e §4.2 (a linha do RISK_PARITY com restrições: sequência de QPs sobre o poliedro mais Lagrangiano aumentado só para a volatilidade);
  - citar nas docs as fontes verificadas da pesquisa: Tibshirani 1996 §6 (a família do leverage no Goldfarb-Idnani), Niedermayer & Niedermayer 2007 Prop. 1/Cor. 1 (o caminho do `MAX_RETURN`), a prova da homogeneização com leverage (κ ≥ 1/(L‖v‖∞), o que exclui `UNBOUNDED`), a diferença |Σ| × |Rᵀ||R| na guarda do κ, e as referências do plano 4; a tabela "Which setters and data each objective takes" ganha a linha nova.
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

---

## Desvios do plano (registro)

O que mudou durante a execução, em relação ao texto das Tasks 2 a 7 acima. Onde este registro e o texto divergem, vale
o registro; as medições estão em `docs/guide/validation.md` e a pesquisa em `docs/research/2026-11-16-risk-parity-and-numerics.md`.

1. **Tasks 3 e 4 refeitas a partir da pesquisa.** O L-BFGS-B projetado (`box_lbfgs`) foi trocado por funções do objetivo
   (f, ∇f, J, Hessiana exata) e uma sequência de QPs de Gauss-Newton proximal sobre o poliedro {orçamento, target_return,
   bounds} com o `qp::solve` existente. O Lagrangiano aumentado ficou só para vol(w) = s (regra de ρ de ABMS 2007, Alg. 3.1,
   com a atribuição do plano a N&W §17.4 corrigida), e a Task 4 ganhou o polimento por Newton no sistema KKT.
2. **Peso proximal τ.** A fórmula do plano, τ = 0,05·tr(Σ)/(2n), mistura unidades (J'J é homogênea de grau 0 em Σ):
   com escala 1 e cond(Σ) = 1e8 convergiram 22 de 300 execuções, e com escala 1e-4 o `qp::solve` lançou exceção em 16 de
   300. Foi substituída por τ₀ = 0,05·tr(2J'J)/n, sem escala, com atualização de Levenberg-Marquardt (τ /= 3 após passo
   cheio, τ *= 4 após encurtar, piso 100·n·ε·max (2J'J)ᵢᵢ) e Cholesky próprio de Q (τ *= 10 em vez de exceção).
3. **Produto Σw compensado.** A soma em `long double` não vale em MSVC (`long double` = `double`), então Σw e os produtos
   do KKT usam o produto escalar compensado Dot2 (Ogita-Rump-Oishi 2005, `std::fma` e TwoSum) em `double`. O oráculo dos
   testes também deixou de depender de precisão estendida (as partes que ainda a usam são puladas onde ela não existe).
4. **`-ffp-contract=off`.** O GCC 16 com `-mfma` fundia a·b+c e quebrava a compensação; a opção entrou em `portopt_core`
   para GCC e Clang.
5. **Piso de arredondamento derivado no verificador de KKT.** A banda de estacionaridade do plano é inalcançável quando o
   resíduo tem o seu próprio piso (H·dw com H ≈ 2J'J, que cresce como cond^2,3). O verificador usa um piso derivado
   (κ = 2) da perturbação de um arredondamento de cada peso e de cada termo, e a decisão `OPTIMAL` do solver é a mesma
   medida (fator 16, parada interna com fator 2), mais a linha do orçamento resolvida (exposição bruta ≤ ~4,5e3).
   Nenhuma banda foi afrouxada.
6. **Sem carregar (λ, ρ) para o segundo esforço.** O plano previa continuar o multiplicador; a medição mostrou que
   carregar só λ ou o par transforma 9 ou 10 dos 72 `OPTIMAL` em `MAX_ITER` (todos n = 50), então cada esforço reinicia o método.
7. **s = σ_min.** Com |s − σ_min| ≤ 1e-12·s (`qp::kRelativeTolerance`, promovida de `kMinimumTolerance`), o único ponto viável
   é a carteira de mínima variância: `OPTIMAL` com λ_vol = −∞ e os demais multiplicadores NaN (a LICQ falha), como o
   `MAX_RETURN` na fronteira vertical.
8. **`MAX_ITER` sem ponto, uniforme.** Pesos, multiplicadores e resíduo NaN em todo objetivo; inclui target_return sem
   bounds (sem minimizador) e a linha do orçamento não resolvível em `double`.
9. **Limitação da precisão dos pesos.** Com vol-alvo e poucos bounds ativos em cond(Σ) ≈ 1e8 os pesos só são determinados até
   cerca de 1e-4 relativo (medido até 3,7e-5; o valor do objetivo fica dentro da banda). Um problema só com o orçamento e
   bounds (Capped, orçamento igual, n = 20, cond 1e8, semente 31696) é `MAX_ITER`. Documentado, não corrigido. Custo: n = 50,
   cond 1e8, sem box e vol-alvo leva cerca de 10 s por solve; n ≤ 20, no máximo cerca de 0,15 s.
10. **Tasks 6 e 7 numa só entrega.** Sem mudança de interface. Os testes de `MAX_ITER` de objetivos QP não alcançam o limite
    de iterações em problemas de teste, então a regra uniforme é verificada em `to_result`, onde ela é feita
    (`CapiSolve.MaxIterHasNoPointForQpObjectives`). Os testes de métricas do C ficaram em `CapiMetrics.*`.
