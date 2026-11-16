# PortOpt — Plano 3: correção do MAX_RETURN, leverage, MAX_SHARPE e MAX_DIVERSIFICATION

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fechar o achado aberto do plano 2 (MAX_RETURN devolvendo `OPTIMAL` errado em casos mal condicionados, e o QP aceitando um ótimo que viola bound), montar uma bateria de fuzz com referências independentes, e entregar o leverage Σ|w| ≤ L em todos os objetivos de QP, `MAX_SHARPE` e `MAX_DIVERSIFICATION` por homogeneização, com o status `UNBOUNDED` e as métricas `sharpe_ratio` e `diversification_ratio`.

**Architecture:**
- **Leverage:** é uma família implícita de linhas no Goldfarb-Idnani, `a'x − ‖x‖₁ ≥ b`. O solver gera sob demanda a linha mais violada, `(a − sign(x))'x ≥ b`, e a Hessiana continua definida positiva.
- **MAX_SHARPE e MAX_DIVERSIFICATION:** resolvidos num construtor homogeneizado próprio, `min y'Σy` s.a. `v'y = 1`, com `v = μ̃` ou `v = σ`. As linhas homogêneas são `yᵢ − lo·1'y ≥ 0`, `hi·1'y − yᵢ ≥ 0` e `1'y ≥ 0`, e a família homogênea do leverage é `L·1'y − ‖y‖₁ ≥ 0`.
- **Volta aos pesos do usuário:** w = y/1'y, e os multiplicadores são mapeados por mudança de escala exata (fórmula na Task 3).
- **MAX_RETURN:** o caminho paramétrico passa a reconhecer linhas geradas no conjunto ativo.

**Tech Stack:** o dos planos 1–2.

**Spec:** `docs/superpowers/specs/2026-11-15-portopt-design.md`. As seções que importam aqui são §3.2 leverage, §3.3 MAX_SHARPE e MAX_DIVERSIFICATION (caixa com peso 0, `INFEASIBLE`, `UNBOUNDED`), §4.2 homogeneização e família implícita, e §5.2 e §5.3.

## Roteiro dos planos

1. Fundação + `MIN_VARIANCE` (concluído).
2. Objetivos de QP (concluído; um achado aberto, Task 0 deste plano).
3. **Correção do MAX_RETURN, leverage, homogeneização** (este).
4. `RISK_PARITY`: Spinu, Lagrangiano aumentado, `target_volatility`, `risk_budget`; métrica `risk_contributions`.
5. `cpp/portopt.hpp`, benchmarks, wheels/PyPI, release.

## Global Constraints

- Todas as dos planos 1 e 2 continuam valendo (`docs/superpowers/plans/2026-11-15-portopt-01-fundacao.md` e `…-02-objetivos-qp.md`, seções Global Constraints). Em especial:
  - a interface C só cresce por acréscimo;
  - a positividade definida é decidida uma vez (fator de `Data`, √c·L);
  - a convenção de multiplicadores do header vale: cada multiplicador é ∂f*/∂(lado direito); uma restrição "≤" tem multiplicador ≤ 0;
  - commits **sem** trailer `Co-Authored-By`.
- **Lição do plano 2: tolerância é medida, nunca chutada.**
  - Toda tolerância nova de decisão no código (constante, segmento, posto, ilimitado) é relativa ao ruído do problema e fica justificada em `docs/guide/validation.md` com a medição que a sustenta.
  - Todo teste roda nas condições da bateria de fuzz da Task 0: cond(Σ) ∈ {1e2, 1e4, 1e6, 1e8}, escala diária e anual, médias iguais, rf acima de todos os μ, bounds apertados.
- **Tolerâncias de teste por condicionamento** (KKT e concordância com a referência, relativas):

  | cond(Σ) | tolerância |
  |---|---|
  | ≤ 1e4 | 1e-10 |
  | 1e6 | 1e-9 |
  | 1e8 | 1e-7 |

  Um teste que não atinge a sua faixa é reportado com os resíduos. Nunca se afrouxa uma faixa.
- **Leverage:**
  - Σ|wᵢ| ≤ L, só sobre os ativos de risco; o caixa fica fora.
  - Setter `leverage(L)` com L > 0 e finito.
  - Sem caixa, L < 1 é inviável (Σ|w| ≥ |1'w| = 1) e o QP devolve `INFEASIBLE`.
  - Vale para MIN_VARIANCE, MAX_RETURN, MAX_SHARPE, UTILITY, IMPORTANCE, KELLY, TRACKING_ERROR e MAX_DIVERSIFICATION; não vale para RISK_PARITY.
  - Multiplicador `PO_CONSTRAINT_LEVERAGE` = λ_lev ≤ 0, com ∇f ∋ λ_lev·g, g ∈ ∂‖w‖₁. É NaN sem leverage e 0 com leverage inativo.
- **Formas de minimização** (vão para o header e `math.md`), com μ̃ = μ − rf·1, rf = 0 sem caixa, e σ = √diag Σ:
  - `MAX_SHARPE`: f = −μ̃'w / √(w'Σw).
  - `MAX_DIVERSIFICATION`: f = −σ'w / √(w'Σw).
  - Com caixa, nos dois objetivos o peso do caixa é fixo em 0. O termo de orçamento sai como β₀ (convenção ∇f = … − β₀·1) e o multiplicador de orçamento é NaN.
- **Matriz objetivo × setter**, acréscimos ao plano 2:

  | Objetivo | Exige | Aceita além de bounds |
  |---|---|---|
  | MAX_SHARPE | mu, covariance | leverage |
  | MAX_DIVERSIFICATION | covariance | leverage |

  A coluna `leverage` é aceita em todos os objetivos acima e rejeitada em nenhum dos atuais. A mensagem segue o padrão `"leverage does not apply to <OBJECTIVE>"` (só RISK_PARITY, no plano 4).
- **Status `UNBOUNDED`:** quando o supremo existe mas não é atingido (a linha `1'y ≥ 0` fica ativa no ótimo homogeneizado). Pesos e multiplicadores saem NaN.
- **`INFEASIBLE` em MAX_SHARPE e MAX_DIVERSIFICATION:**
  - quando nenhuma carteira viável tem v'w > 0 (o QP homogeneizado é inviável);
  - com caixa, quando os bounds do caixa excluem 0.

## Review Focus

1. **Leverage com pesos exatamente zero.** sign(0) é qualquer valor em [−1, 1]. O solver gera s = +1 nesses índices, e o helper de KKT aceita qualquer gᵢ ∈ [−1, 1] onde wᵢ = 0. Teste: `Leverage.ZeroWeightsSubgradient` (Task 1).
2. **MAX_SHARPE com rf acima de todos os μ**, com ou sem bounds:
   - long-only, ou com qualquer cone em que v'w ≤ 0: `INFEASIBLE`;
   - venda a descoberto sem bounds e 1'Σ⁻¹μ̃ ≤ 0: `UNBOUNDED`.

   Testes: `MaxSharpe.RiskFreeAboveAllMeansLongOnly` e `MaxSharpe.UnboundedWithShorts` (Task 3).
3. **Leverage ativo no MAX_RETURN.** O conjunto ativo ao longo do caminho inclui linhas geradas, e o sinal delas muda quando um peso cruza zero. Teste: `MaxReturn.LeverageFuzz` (Task 2).
4. **Ótimo do QP que viola bound em t grande.** O achado do plano 2 vira teste: `GoldfarbIdnani.FinalFeasibilityAtLargeScale` (Task 0).
5. **Python com `leverage` e bounds em array**: a igualdade e o hash do `Spec` continuam valendo (Task 5).

---

### Task 0: Correção do MAX_RETURN, viabilidade final do QP e bateria de fuzz

**Files:**
- Modify: `src/problems/parametric_path.cpp` (`weight_noise`), `src/problems/max_return.cpp` (verificação antes de `at_limit`), `src/qp/goldfarb_idnani.cpp` (re-checagem final)
- Create: `tests/unit/fuzz/random_problems.hpp` (gerador), `tests/unit/fuzz/frontier_reference.hpp` (referência), `tests/unit/fuzz/max_return_fuzz_test.cpp`
- Modify: `tests/unit/qp/goldfarb_idnani_test.cpp`, `CMakeLists.txt`, `docs/guide/validation.md`

**Interfaces:**
- Consumes: plano 2 (`Path`, `PathPoint`, `weights_constant`, `at_limit`, `qp::solve`, `user_kkt.hpp`, `random_covariance.hpp`).
- Produces:
  ```cpp
  namespace po::testing {
  struct ProblemShape
  {
      std::size_t n;
      double condition;        // cond(Sigma)
      double scale;            // 1 (annual) or 1e-4 (daily) on Sigma; sqrt-scaled on mu
      enum class Means { Spread, Equal, NearlyEqual, BelowRiskFree } means;
      enum class Box { None, LongOnly, Tight, Fixed } box;   // Tight: hi = 2/n; Fixed: some lo == hi
      bool cash;
      std::uint64_t seed;
  };
  struct GeneratedProblem { po::Spec spec; po::Data data; };
  GeneratedProblem generate(const ProblemShape& shape);
  std::vector<ProblemShape> fuzz_shapes(std::size_t count, std::uint64_t seed);   // covers every combination
  // MAX_RETURN reference by another route: the largest target return r with MIN_VARIANCE(target_return = r)
  // feasible and of volatility <= s, by bisection on r to 1e-13 relative; returns NaN when infeasible.
  double max_return_reference(const po::Spec& spec, const po::Data& data, double cap);
  }
  ```
  `random_covariance(n, seed, condition)` já existe; `generate` multiplica por `scale` e escolhe μ pelo padrão. `NearlyEqual` é espalhamento de 1e-8 a 1e-5 relativo, e `BelowRiskFree` é rf acima de todo μ.

Correções (rulings do plano 2):
1. `weight_noise` usa `t·‖μ̃‖∞`, não `t·|Σ⁻¹μ̃|max`: o ruído da equação de estacionaridade, que não carrega cond(Σ).
2. Antes de `at_limit` (teto inativo), verificar no ponto: viabilidade primal de toda linha (relativa a max(1, |cᵢ|, |Cᵢ||x|)) e KKT do usuário com λ_vol = 0. Se falhar, a decisão "constante" estava errada: continua o `climb`, ou `MAX_ITER` ao fim das dobras.
3. Goldfarb-Idnani: ao não achar candidata, re-checar toda desigualdade com a tolerância calculada no x **final** (não no max|x| das iterações). Se alguma for violada, ela volta a ser candidata e o algoritmo continua.

- [ ] **Step 1: Testes**
  - `GoldfarbIdnani.FinalFeasibilityAtLargeScale`: o caso do re-review do plano 2 (G = Σ, g = −t μ̃ com t = 1e12, long-only, n = 30, uma variância 1e-10) → `Optimal` sem nenhuma linha violada além da tolerância no x final; `expect_kkt`.
  - `MaxReturnFuzz.AgainstFrontierReference`: `fuzz_shapes(400, 1)` com `box` ∈ {None, LongOnly, Tight}, caps em 5 frações entre σ_min e max σᵢ. Para cada caso, o status é `OPTIMAL` (ou `INFEASIBLE` exatamente quando a referência é NaN); o retorno concorda com `max_return_reference` na faixa de tolerância do cond; o `expect_user_kkt` passa na mesma faixa; e vol ≤ s.
  - `MaxReturnFuzz.NeverSilentlyWrong`: os casos `NearlyEqual` com uma variância de 1e-10 a 1e-6 (os do re-review). O resultado é `OPTIMAL` correto pela referência, ou `MAX_ITER`. Nunca `OPTIMAL` fora da referência, e nunca um peso fora dos bounds.
  - Os testes `MaxReturn.*` existentes passam sem mudar valores.
- [ ] **Step 2:** Rodar, ver falhar (o `NeverSilentlyWrong` e o `FinalFeasibilityAtLargeScale` falham hoje).
- [ ] **Step 3:** Implementar as três correções e o gerador.
- [ ] **Step 4:** PASS em `gcc-release` e `gcc-debug`. Rodar os scripts `scratchpad/rereview/*.py` e `scratchpad/final/*.py` do plano 2 contra a DLL nova: todos limpos. Registrar os números em `validation.md`.
- [ ] **Step 5: Commit** — `fix: MAX_RETURN constant-weight noise and QP final feasibility; fuzz harness`

---

### Task 1: Família implícita do leverage no QP e nos objetivos de QP

**Files:**
- Modify: `src/qp/problem.hpp`, `src/qp/goldfarb_idnani.cpp`, `src/qp/active_set.hpp`
- Modify: `src/core/spec.hpp` (`std::optional<double> leverage`), `src/core/result.hpp` (comentário: termo λ_lev·g)
- Modify: `src/problems/constraints.hpp`, `src/problems/constraints.cpp` (a família e o multiplicador), `src/problems/requirements.cpp`
- Modify: `tests/unit/kkt.hpp`, `tests/unit/user_kkt.hpp`
- Create: `tests/unit/fuzz/sign_pattern_reference.hpp`, `tests/unit/fuzz/leverage_fuzz_test.cpp`
- Test: `tests/unit/qp/goldfarb_idnani_test.cpp`, `tests/unit/problems/leverage_test.cpp`

**Interfaces:**
- Consumes: Task 0.
- Produces:
  ```cpp
  namespace po::qp {
  // a'x - ||x||_1 >= b, the implicit family of rows (a - s)'x >= b over s in {-1, 1}^n.
  struct GrossExposure { std::vector<double> offset; double rhs = 0.0; };
  // Problem: std::optional<GrossExposure> gross_exposure;
  // Solution (Optimal only): the active generated rows and their multipliers, so that
  //   G x + g = E'y + C'z + sum_j zeta_j (a - s_j),  zeta_j >= 0.
  // std::vector<double> gross_exposure_signs;        // k x n, entries -1 or +1
  // std::vector<double> gross_exposure_multipliers;  // k
  }
  // problems::ConstraintRows gains: bool gross_exposure (the spec set leverage); to_result fills
  // Result::leverage_multiplier = -sum_j zeta_j / multiplier_scale (NaN without leverage, 0 when inactive).
  namespace po::testing {
  // Brute force for n <= 8: for every sign pattern s, the QP with s_i x_i >= 0 and s'x <= L added; the best
  // Optimal value. An independent route to the leverage-constrained optimum.
  std::optional<double> sign_pattern_optimum(const po::qp::Problem& problem_without_family, double leverage);
  }
  ```
  - **Geração no GI:** a violação da família é b − a'x + ‖x‖₁. Ela concorre com a desigualdade mais violada pela mesma regra de tolerância, calculada sobre a linha gerada. A linha gerada é (a − s)'x ≥ b, com sᵢ = sign(xᵢ) e s = +1 onde xᵢ = 0.
  - **Re-geração:** uma linha gerada que sai do conjunto ativo é descartada, e pode voltar a ser gerada depois.
  - **Linhas geradas:** são numeradas depois das explícitas para fins de conjunto ativo.
  - **Identidade do conjunto ativo** (usada no caminho da Task 2): índices explícitos mais os vetores de sinal das linhas geradas ativas.
  - **Objetivos de QP:** para MIN_VARIANCE, UTILITY, IMPORTANCE, KELLY, TRACKING_ERROR e MAX_RETURN, `add_constraint_rows` põe `gross_exposure = {0, −L}` quando a spec tem leverage.

- [ ] **Step 1: Testes**
  - `GoldfarbIdnani.GrossExposureMatchesSignPatterns`: sementes 1..100, n em 2..8, G e g aleatórios, 1'x = 1, L ∈ [1, 3]. O valor é igual a `sign_pattern_optimum` em 1e-10 relativo, `expect_kkt` passa com a família, e Σ|x| ≤ L(1 + 1e-12).
  - `GoldfarbIdnani.GrossExposureInactive`: L grande. O resultado é bit a bit igual ao problema sem família, sem linhas geradas.
  - `Leverage.MinVarianceAndUtilityKkt`: sem bounds, L = 1.3. `expect_user_kkt` com o termo λ_lev·g: onde wᵢ ≠ 0, o resíduo é λ_lev·sign(wᵢ); onde wᵢ = 0, o resíduo tem módulo ≤ |λ_lev|. λ_lev ≤ 0.
  - `Leverage.ZeroWeightsSubgradient`: um caso construído com wᵢ = 0 no ótimo e leverage ativo.
  - `Leverage.BelowOneWithoutCashIsInfeasible` (L = 0.9), `Leverage.CashExcluded` (com caixa, L conta só os de risco), `Leverage.MultiplierIsDerivativeOfOptimum` ((f*(L+h) − f*(L−h))/2h = λ_lev em 1e-6), `Leverage.MultiplierNaNWithoutLeverage`.
  - `LeverageFuzz.AgainstSignPatterns`: `fuzz_shapes` com n ≤ 8, os seis objetivos de QP (MAX_RETURN entra na Task 2), concordância na faixa do cond.
  - `Requirements.LeverageAcceptedByQpObjectives`.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** PASS em ambos os presets, Task 0 incluída.
- [ ] **Step 5: Commit** — `feat: gross exposure as an implicit constraint family; leverage for the QP objectives`

---

### Task 2: Leverage no caminho do MAX_RETURN

**Files:**
- Modify: `src/problems/parametric_path.hpp`, `src/problems/parametric_path.cpp`, `src/problems/max_return.cpp`
- Test: `tests/unit/problems/max_return_test.cpp`, `tests/unit/fuzz/max_return_fuzz_test.cpp`

**Interfaces:**
- Consumes: Tasks 0–1 (`gross_exposure_signs`, identidade do conjunto ativo).
- Produces: `same_segment`, `weights_constant`, `multiplier_slope`, `vertex_exit`, `on_segment` e `slope_solution` tratam as linhas geradas como linhas ativas comuns, identificadas pelo vetor de sinal. Mudança de sinal de uma linha gerada é mudança de segmento. O multiplicador de leverage do usuário é λ_lev = −Σζⱼ/t*, ou o limite da inclinação no teto inativo.

- [ ] **Step 1: Testes**
  - `MaxReturn.LeverageActiveAlongPath`: sem bounds, L = 1.5, s alto o bastante para o leverage travar antes da vol. O resultado é `OPTIMAL` com λ_vol = 0, λ_lev < 0 e KKT.
  - `MaxReturn.LeverageAndVolatilityBothActive`: os dois multiplicadores < 0, KKT e diferenças finitas em s e em L.
  - `MaxReturnFuzz.LeverageFuzz`: os casos de `fuzz_shapes` com leverage. A referência é `max_return_reference` com a MIN_VARIANCE também com o leverage, na faixa do cond.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** PASS.
- [ ] **Step 5: Commit** — `feat: MAX_RETURN path with the gross exposure family`

---

### Task 3: MAX_SHARPE e MAX_DIVERSIFICATION por homogeneização

**Files:**
- Create: `src/problems/homogeneous.hpp`, `src/problems/homogeneous.cpp` (o oitavo módulo de `src/problems/`)
- Modify: `src/problems/solve.cpp`, `src/problems/requirements.cpp`, `CMakeLists.txt`
- Create: `tests/unit/fuzz/ratio_reference.hpp`
- Test: `tests/unit/problems/max_sharpe_test.cpp`, `tests/unit/fuzz/ratio_fuzz_test.cpp`

**Interfaces:**
- Consumes: Tasks 0–2; `Data::covariance_factor()`.
- Produces:
  ```cpp
  namespace po::problems {
  // min y' Sigma y  s.t.  v'y = 1,  y_i - lo_i 1'y >= 0,  hi_i 1'y - y_i >= 0 (finite bounds),  1'y >= 0,
  // and with leverage L 1'y - ||y||_1 >= 0 (gross_exposure {L 1, 0}); G = 2 Sigma, factor sqrt(2) L.
  // v = mu~ (MAX_SHARPE) or sigma (MAX_DIVERSIFICATION).
  Result solve_max_sharpe(const Spec& spec, const Data& data);
  Result solve_max_diversification(const Spec& spec, const Data& data);
  }
  namespace po::testing {
  // The ratio maximum by another route: golden-section search on the target return r of
  // MIN_VARIANCE(target_return = r, same bounds and leverage) over the efficient branch, maximizing
  // (r - rf) / vol(r); for MAX_DIVERSIFICATION the same with mu := sigma and rf := 0. NaN when no r gives a
  // positive ratio.
  double max_ratio_reference(const po::Spec& spec, const po::Data& data, std::span<const double> v);
  }
  ```
  Regras:
  - **Bounds e caixa:** os bounds homogeneizados usam os bounds dos ativos de risco. Com caixa, o peso do caixa é fixo em 0. Se os bounds do caixa excluem 0, o status é `INFEASIBLE` antes de resolver.
  - **Status:**
    - QP `Infeasible` → `INFEASIBLE`;
    - linha `1'y ≥ 0` no conjunto ativo do ótimo → `UNBOUNDED`, pesos e multiplicadores NaN;
    - senão `OPTIMAL`, com w = y/κ, κ = 1'y.
  - **Multiplicadores no w do usuário (exatos):** com σ_y = √(y'Σy) e c = κ/(2σ_y³):
    - linha inferior i: βᵢ += c·zᵢ;
    - linha superior j: βⱼ −= c·zⱼ;
    - λ_lev = −c·Σζ;
    - termo de orçamento λ = c·(η − Σ_inf zᵢ·loᵢ + Σ_sup zⱼ·hiⱼ + L·Σζ), onde η é o multiplicador de `1'y ≥ 0`, nulo quando inativo.

    Sem caixa, λ é o `budget_multiplier`. Com caixa, o budget é NaN e β₀ = −λ. A derivação, a partir de ∇f(w) = −κv/σ_y + κΣy/σ_y³ e 2Σy = 2σ_y²·v + C'z + η·1 + Σζⱼ(L·1 − sⱼ), vai em `math.md`.

- [ ] **Step 1: Testes** (todo `OPTIMAL` chama `expect_user_kkt` com ∇f do teste: −v/σ_w + (v'w)Σw/σ_w³)
  - `MaxSharpe.UnconstrainedIsTangency`: sem bounds e com 1'Σ⁻¹μ̃ > 0, w = Σ⁻¹μ̃ / 1'Σ⁻¹μ̃ em 1e-12.
  - `MaxSharpe.CashWeightZero`: com rf, o caixa tem peso 0 exato, o budget é NaN e β₀ é finito com KKT. Com `bounds(0.05, 0.3)` mais rf, o status é `INFEASIBLE`.
  - `MaxSharpe.RiskFreeAboveAllMeansLongOnly` → `INFEASIBLE`. `MaxSharpe.UnboundedWithShorts` (sem bounds, 1'Σ⁻¹μ̃ ≤ 0) → `UNBOUNDED`, pesos NaN.
  - `MaxSharpe.LongOnlyAndLeverageKkt`.
  - `MaxSharpe.NoSharperFeasiblePortfolio`: 10.000 carteiras viáveis aleatórias (semente fixa) com Sharpe ≤ ótimo·(1 + 1e-12).
  - `MaxDiversification.UnconstrainedClosedForm`: w ∝ Σ⁻¹σ. `MaxDiversification.LongOnlyKkt`, `MaxDiversification.CashWeightZero`.
  - `RatioFuzz.AgainstFrontierReference`: `fuzz_shapes(300, 3)` nos dois objetivos, com e sem leverage. Concordância do valor da razão com `max_ratio_reference` na faixa do cond, e status igual ao da referência (NaN ↔ `INFEASIBLE`/`UNBOUNDED`).
  - `Requirements.MaxSharpeRequiresMu`, `Requirements.MaxDiversificationNeedsNoMu`, e `Requirements.UnimplementedObjectivesStayUnsupported` (agora só RISK_PARITY).
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** PASS.
- [ ] **Step 5: Commit** — `feat: MAX_SHARPE and MAX_DIVERSIFICATION by homogenization; UNBOUNDED`

---

### Task 4: Interface C e métricas

**Files:**
- Modify: `include/portopt/portopt.h`, `src/capi/spec_api.cpp`, `src/capi/metrics_api.cpp`, `src/metrics/portfolio.hpp`, `src/metrics/portfolio.cpp`
- Test: `tests/unit/capi/spec_test.cpp`, `tests/unit/capi/solve_test.cpp`, `tests/unit/capi/metrics_test.cpp`, `tests/unit/metrics/portfolio_test.cpp`

**Interfaces:**
- Produces:
  ```c
  PO_API po_status po_spec_set_leverage(po_spec* s, double max_gross);   /* max_gross > 0, finite */
  /* (return - rf) / volatility over the risky weights (rf = 0 without cash); NaN when the volatility is 0. */
  PO_API po_status po_sharpe_ratio(const po_data* d, const double* w, size_t n, double* out);
  /* sigma'w / volatility over the risky weights; needs the covariance; NaN when the volatility is 0. */
  PO_API po_status po_diversification_ratio(const po_data* d, const double* w, size_t n, double* out);
  ```
  ```cpp
  namespace po::metrics {
  double sharpe_ratio(const Data& data, std::span<const double> w);
  double diversification_ratio(const Data& data, std::span<const double> w);
  }
  ```
  O header documenta:
  - os objetivos resolvidos: todos menos RISK_PARITY;
  - as formas f de MAX_SHARPE e MAX_DIVERSIFICATION;
  - o caixa em peso 0, com β₀ e budget NaN;
  - `PO_SOLVE_UNBOUNDED`, com pesos e multiplicadores NaN;
  - `PO_CONSTRAINT_LEVERAGE` (λ_lev ≤ 0, NaN sem leverage, 0 inativo), com o termo λ_lev·g, g ∈ ∂‖w‖₁, na fórmula de estacionaridade;
  - a coluna leverage na matriz.

- [ ] **Step 1: Testes**
  - `CapiSpec.LeverageDomain`.
  - `CapiSolve.SharpeAndDiversificationEndToEnd`: bit a bit igual a `po::solve`, com e sem caixa e com leverage.
  - `CapiSolve.UnboundedIsResultNotError`.
  - `CapiSolve.LeverageMultiplierOnlyWithLeverage`.
  - `Metrics.SharpeAndDiversification` (à mão; caixa com peso; vol 0 → NaN).
  - `Metrics.SharpeAtOptimumIsMaximal`: a métrica no ótimo do MAX_SHARPE é igual a −f*.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** `ctest` PASS (inclui `c_header_check` e `style`).
- [ ] **Step 5: Commit** — `feat(capi): leverage, MAX_SHARPE, MAX_DIVERSIFICATION, ratio metrics`

---

### Task 5: Python e documentos

**Files:**
- Modify: `python/src/portopt/_native.py`, `_spec.py` (campo `leverage`), `_data.py` (`sharpe_ratio`, `diversification_ratio`), `_result.py` (`leverage_multiplier`; status `UNBOUNDED`)
- Test: `python/tests/test_ratio_objectives.py`
- Modify: `docs/guide/math.md`, `docs/guide/validation.md`, `ARCHITECTURE.md`, `README.md`, `CHANGELOG.md`

**Interfaces:**
- Produces: o campo `Spec.leverage: float | None = None`; os métodos `Data.sharpe_ratio(weights)` e `Data.diversification_ratio(weights)`; o campo `Result.leverage_multiplier: float`.

- [ ] **Step 1: Testes**:
  - `test_tangency_closed_form` (numpy);
  - `test_sharpe_cash_weight_zero`;
  - `test_sharpe_unbounded_status`;
  - `test_diversification_closed_form`;
  - `test_leverage_caps_gross_exposure`;
  - `test_leverage_multiplier_nan_without_leverage`;
  - `test_spec_with_leverage_and_array_bounds_hashes`;
  - `test_ratio_metrics_match_numpy`.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar.
- [ ] **Step 4:** Documentação:
  - `math.md`: a família implícita, a homogeneização com a derivação dos multiplicadores, `UNBOUNDED`, o caixa em peso 0 e o leverage no caminho do MAX_RETURN.
  - `validation.md`: a bateria de fuzz, as referências independentes (fronteira por bissecção, golden-section, padrões de sinal), as faixas de tolerância por cond e a medição que justifica cada tolerância de decisão.
  - `ARCHITECTURE.md`: `homogeneous`, a família no QP, `tests/unit/fuzz`.
  - `README.md`: um exemplo de MAX_SHARPE com leverage.
  - `CHANGELOG.md`.
- [ ] **Step 5:** Verificar:
  - pytest PASS;
  - mypy strict (`--config-file python/pyproject.toml`) limpo;
  - `check_style.py` limpo;
  - `ctest` PASS nos dois presets.
- [ ] **Step 6: Commit** — `feat(python): leverage and ratio objectives; docs for plan 3`
