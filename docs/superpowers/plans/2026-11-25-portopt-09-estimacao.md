# PortOpt — Plano 9: utilidades de estimação

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Momentos amostrais, Bayes-Stein e Black-Litterman como utilidades da API C, ao lado do core, com os bindings Python, C# e C++, o pivô na mensagem de Σ não PD e a conversão do posterior para a entrada da otimização.

**Architecture:** Um módulo novo, `src/estimation/`, depende só de `linalg` e `core` (`Data`, `Failure`). Σ chega por um `po_data` já validado, e as visões chegam por um handle novo, `po_views`. As entradas C ficam em `src/capi/estimation/`. Python e C# escondem os handles atrás de funções puras sobre arrays.

**Tech Stack:** C++23, CMake 3.28, GoogleTest, MinGW-GCC (`gcc-release`), MSVC (`native-release`), Python 3.11–3.14 (ctypes, numpy, pytest, mypy, uv), .NET 10 (xUnit, `LibraryImport`).

**Spec:** `docs/superpowers/specs/2026-11-25-portopt-estimation-design.md` ("spec §N").

## Global Constraints

- Valem as dos planos 1 a 8: `STYLEGUIDE.md` (sem `else`; no máximo 60 linhas, 5 parâmetros e 3 níveis de aninhamento por função, **a interface C e o Python incluídos**; no máximo 8 módulos por pasta), commits **sem** trailer `Co-Authored-By` (CONTRIBUTING.md), `-ffp-contract=off`, sem `long double` nem `-ffast-math`, `docs/superpowers/` em português, documentos do repositório em inglês.
- `src/estimation/` só usa `+ - * /` e `sqrt`. Nada de libm.
- A interface C só cresce: nenhuma assinatura existente muda, e `PO_ABI_VERSION` fica em 1. `PO_VERSION` no CMake fica em 1.0.0, e as entradas vão para `Unreleased` e `PublicAPI.Unshipped.txt` (spec §2).
- Erros: `PO_INVALID_ARGUMENT` (ou `PO_INVALID_HANDLE` para handle nulo), com o motivo em `po_last_error()`. Nenhum código novo. Em uma falha, nenhuma saída é escrita.
- O core (`solve`, `problems`, `qp`, `riskparity`, `metrics`) não inclui nada de `estimation`.
- Testes locais: `cmake --build --preset gcc-release && ctest --preset gcc-release`, e antes do último commit de cada task C++ também `native-release` (MSVC). Não rodar `gcc-debug`.
- Python: `PORTOPT_LIB=<lib construída> uv run --project python --python 3.14 pytest python/tests -q` e o `mypy` do CI (`.github/workflows/ci.yml`).
- C#: `PORTOPT_LIB=<lib construída> dotnet test --solution csharp/PortOpt.slnx -c Release`.
- Proibido para deixar o CI verde: `ctest -E`, `GTEST_SKIP`, `NOLINT`, afrouxar tolerância de teste existente.

## Review Focus

1. **Escala diária × anualizada no Bayes-Stein.** O DataPanel vai chamar com retornos diários. `w_hat` tem de ser invariante a multiplicar os retornos por uma constante e tem de mudar quando `mu` e `S` são multiplicados por 252. (Task 3: `BayesStein.WeightIgnoresTheUnit`, `BayesStein.WeightDependsOnTheFrequency`.)
2. **Visões com confiança total (Ω = 0) e uma visão quase redundante.** Com Ω = 0, `P mu = Q` dentro de `1e-15`. P com linhas quase paralelas, mas de posto cheio, não pode dar NaN em silêncio: ou a conta sai finita, ou o Cholesky de `A` recusa com o pivô. (Task 5: `BlackLitterman.FullConfidenceReproducesTheView`, `BlackLitterman.NearlyRedundantViewsAreFiniteOrRejected`.)
3. **`po_data` com rf passado a uma utilidade.** Tem de ser recusado, não ignorado. (Task 3: `EstimationInputs.RiskFreeDataIsRejected`.)
4. **Pesos de mercado com soma de ponto flutuante** (por exemplo `1/3` três vezes, ou capitalizações normalizadas). Devem passar; uma soma de 0.99 deve falhar. (Task 4: `Equilibrium.FloatingPointSumIsAccepted`, `Equilibrium.WeightsOffByOnePercentAreRejected`.)
5. **Python com `p` de uma dimensão** (`[1, -1]` em vez de `[[1, -1]]`), o erro natural de quem tem uma visão só. Tem de dar `ValueError` com mensagem clara, não um erro nativo nem uma visão errada. (Task 9: `test_views_reject_one_dimensional_p`.)

---

### Task 1: O pivô de Cholesky na mensagem de Σ não positiva definida

**Files:**
- Modify: `src/linalg/cholesky.hpp`, `src/linalg/cholesky.cpp`
- Modify: `src/core/data.hpp` (expõe `kSymmetryTolerance`), `src/core/data.cpp:15,190-197`
- Test: `tests/unit/linalg/cholesky_test.cpp`, `tests/unit/core/data_test.cpp:56,96,108`, `tests/unit/capi/data_test.cpp:77,83`, `python/tests/test_min_variance.py:89` (continua passando, porque o `match` é por prefixo)

**Interfaces:**
- Produces: `std::expected<std::vector<double>, std::size_t> po::linalg::cholesky_or_pivot(std::span<const double> a, std::size_t n)`, que devolve o fator ou o índice (base zero) do primeiro pivô recusado. `cholesky` continua com a mesma assinatura e o mesmo resultado, implementado sobre ela.
- Produces: `inline constexpr double po::kSymmetryTolerance = 1e-10;` em `core/data.hpp`, no lugar da constante anônima de `data.cpp`.
- Produces: o formato de mensagem `"<nome> is not positive definite: pivot <j> is not positive"`.

- [ ] **Step 1: Testes que falham.** Em `cholesky_test.cpp`: `Cholesky.ReportsTheFirstRejectedPivot`. Uma matriz 4×4 positiva definida nas duas primeiras linhas e colunas, com `a[2][2]` tal que o pivô 2 é negativo (por exemplo `[[4,2,0,0],[2,2,0,0],[0,0,-1,0],[0,0,0,1]]`), dá `error() == 2`. Também um caso com o pivô 1 recusado (cobre a segunda coluna do par de `factor_in_place`) e um caso PD em que `cholesky_or_pivot` e `cholesky` devolvem fatores bit a bit iguais. Em `data_test.cpp` (core e capi), as mensagens esperadas viram `"covariance is not positive definite: pivot <j> is not positive"`, com o `j` que cada caso produz.
- [ ] **Step 2: Rodar e ver falhar.** `ctest --preset gcc-release -R "Cholesky|Data"`. Esperado: falha de compilação (`cholesky_or_pivot` não existe).
- [ ] **Step 3: Implementar.** `finish_column` e `factor_in_place` passam a devolver o índice recusado (`std::optional<std::size_t>`, vazio no sucesso). `cholesky` vira `cholesky_or_pivot(a, n)` convertida para `optional`. `Data::create` usa `cholesky_or_pivot` e monta a mensagem com o índice.
- [ ] **Step 4: Rodar e ver passar.** `ctest --preset gcc-release` inteiro (os chamadores do QP, do Spinu e da descida não mudam), depois `native-release`.
- [ ] **Step 5: Commit.** `feat(core): the covariance message names the rejected pivot`

---

### Task 2: Momentos amostrais

**Files:**
- Create: `src/estimation/moments.hpp`, `src/estimation/moments.cpp`
- Create: `tests/unit/estimation/moments_test.cpp`
- Modify: `CMakeLists.txt` (`portopt_core` ganha `src/estimation/moments.cpp`; `portopt_tests` ganha o teste)

**Interfaces:**
- Produces: `struct po::estimation::Moments { std::vector<double> mu; std::vector<double> covariance; };`
- Produces: `std::expected<Moments, Failure> po::estimation::sample_moments(std::span<const double> returns, std::size_t t, std::size_t n);`

- [ ] **Step 1: Testes que falham.**
  - `Moments.SmallPanelByHand`: `returns = {0.01, 0.02, 0.03, -0.01, -0.01, 0.02}` (t = 3, n = 2) dá `mu = {0.01, 0.01}` e `covariance = {0.0004, -0.0003, -0.0003, 0.0003}`, cada um com `EXPECT_NEAR(..., 1e-18)`.
  - `Moments.CovarianceIsExactlySymmetric`: com um painel 50×4 de `random_covariance.hpp` (ou literais), `cov[i*n+j] == cov[j*n+i]` bit a bit.
  - `Moments.OneObservationIsRejected`: t = 1 dá a mensagem `"returns need at least two observations"`.
  - `Moments.NonFiniteIsRejected`: um NaN em (1, 0) dá `"returns is not finite at (1, 0)"`.
  - `Moments.ZeroAssetsIsRejected`: n = 0 dá `"asset count must be positive"`.
- [ ] **Step 2: Rodar e ver falhar.** `ctest --preset gcc-release -R Moments`. Esperado: falha de compilação.
- [ ] **Step 3: Implementar.** São duas passadas: primeiro a média da coluna, somando em ordem de `t`; depois o triângulo superior dos produtos dos desvios, em ordem de `t`, dividido por `t - 1` e espelhado. O overflow de `t * n` é conferido como em `Data::set_returns` (`"returns size is too large"`).
- [ ] **Step 4: Rodar e ver passar.** `ctest --preset gcc-release -R Moments`, e depois o `ctest` inteiro (a regra de estilo de módulos por pasta inclusive).
- [ ] **Step 5: Commit.** `feat(estimation): sample moments of a return panel`

---

### Task 3: Entradas das utilidades e Bayes-Stein

**Files:**
- Create: `src/estimation/inputs.hpp`, `src/estimation/inputs.cpp`
- Create: `src/estimation/bayes_stein.hpp`, `src/estimation/bayes_stein.cpp`
- Create: `tests/unit/estimation/bayes_stein_test.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Produces: `std::optional<Failure> po::estimation::check_inputs(const Data& data, bool needs_mu);`. Recusa `data` sem covariância (`"data has no covariance"`), sem mu quando `needs_mu` (`"data has no mu"`) e com rf (`"data must not have a risk-free rate"`).
- Produces: `struct po::estimation::BayesStein { std::vector<double> mu; double grand_mean; double weight; };`
- Produces: `std::expected<BayesStein, Failure> po::estimation::bayes_stein(const Data& data, std::size_t t);`

- [ ] **Step 1: Testes que falham.** O fixture é `S = {0.0004, 0.0001, 0.0001, 0.0009}` (n = 2), com t = 60.
  - `BayesStein.ConstantMeanIsItsOwnTarget`: `mu_hat = {0.001, 0.001}` dá `grand_mean` ≈ 0.001 (`1e-18`), `weight` ≈ 1 (`1e-12`) e `mu` ≈ `mu_hat`.
  - `BayesStein.MatchesTheClosedForm`: `mu_hat = {0.001, 0.003}`. O valor esperado sai das fórmulas de spec §6.2 escritas no teste com a inversa 2×2 explícita (a referência por outro caminho).
  - `BayesStein.WeightIgnoresTheUnit`: multiplicar `mu_hat` por 100 e S por 10⁴ mantém `weight` dentro de `1e-12` relativo.
  - `BayesStein.WeightDependsOnTheFrequency`: multiplicar `mu_hat` e S por 252, com o mesmo t, muda `weight` em mais de `1e-3` relativo.
  - `BayesStein.TooFewObservationsIsRejected`: t = 4 com n = 2 dá `"bayes-stein needs t > n + 2; got t = 4 for n = 2"`.
  - `EstimationInputs.RiskFreeDataIsRejected`: `Data` com `set_risk_free(0.01)` dá `"data must not have a risk-free rate"`. O mesmo para sem mu e sem covariância.
- [ ] **Step 2: Rodar e ver falhar.** `ctest --preset gcc-release -R "BayesStein|EstimationInputs"`.
- [ ] **Step 3: Implementar.** As resoluções usam `data.covariance_factor()` com `solve_lower` e `solve_lower_transposed`, e o fator de escala `(t - n - 2) / (t - 1)` aplicado ao resultado (spec §6.2). Sem fatoração nova.
- [ ] **Step 4: Rodar e ver passar.** O `ctest` inteiro.
- [ ] **Step 5: Commit.** `feat(estimation): Bayes-Stein shrinkage of the mean (Jorion 1986)`

---

### Task 4: Visões, equilíbrio e as duas formas de Ω

**Files:**
- Create: `src/estimation/views.hpp`, `src/estimation/views.cpp`
- Create: `src/estimation/black_litterman.hpp`, `src/estimation/black_litterman.cpp` (por ora só o equilíbrio e Ω)
- Create: `tests/unit/estimation/views_test.cpp`, `tests/unit/estimation/black_litterman_test.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: `check_inputs` (Task 3), `kSymmetryTolerance` (Task 1), `linalg::householder_qr` e `has_full_column_rank`.
- Produces: a classe `po::estimation::Views`:
  - `static std::expected<Views, Failure> create(std::size_t k, std::size_t n, std::span<const double> p, std::span<const double> q);`
  - `std::expected<void, Failure> set_omega(std::span<const double> omega);`
  - os acessores `view_count()`, `asset_count()`, `p()`, `q()`, `has_omega()` e `omega()`.
- Produces: `std::expected<std::vector<double>, Failure> equilibrium(const Data& data, std::span<const double> w_mkt, double delta);`
- Produces: `std::expected<std::vector<double>, Failure> omega_he_litterman(const Data& data, const Views& views, double tau);`
- Produces: `std::expected<std::vector<double>, Failure> omega_idzorek(const Data& data, const Views& views, double tau, std::span<const double> confidence);`
- Produces: `inline constexpr double kMarketWeightTolerance = 1e-8;`

O fixture comum é Σ = `{0.04, 0.01, 0.01, 0.09}`, com `P = {1, -1}` (k = 1) e τ = 0.05.

- [ ] **Step 1: Testes que falham.**
  - `Equilibrium.ByHand`: `w = {0.6, 0.4}` e δ = 2.5 dão `Pi = {0.07, 0.105}` (`1e-16`).
  - `Equilibrium.FloatingPointSumIsAccepted`: n = 3, Σ diagonal e `w = {1.0/3, 1.0/3, 1.0/3}` são aceitos.
  - `Equilibrium.WeightsOffByOnePercentAreRejected`: `{0.6, 0.39}` dá `"w_mkt must sum to 1 within 1e-8"`.
  - `Equilibrium.DeltaMustBePositive`: δ = 0 e δ = NaN dão `"delta must be positive and finite"`.
  - `Views.RedundantViewsAreRejected`: `P = {1, -1, 2, -2}` (k = 2, n = 2) dá `"views are redundant: P does not have full row rank"`.
  - `Views.MoreViewsThanAssetsAreRejected`: k = 3 com n = 2 dá `"k must not exceed n"`.
  - `Views.EmptyViewsAreAccepted`: k = 0, com `p` e `q` vazios.
  - `Views.OmegaMustBeSymmetricWithNonNegativeDiagonal`: as mensagens são `"omega is not symmetric at (0, 1)"` e `"omega has a negative diagonal entry at 1"`. Uma falha mantém o Ω anterior.
  - `Omega.HeLitterman`: dá `{0.0055}` (`1e-18`).
  - `Omega.Idzorek`: `c = 0.5` dá `{0.0055}`, `c = 1` dá `{0}` exato e `c = 0.25` dá `{0.0165}`. `c = 0` e `c = 1.5` dão `"confidence must be in (0, 1] at 0"`.
  - `Omega.AssetCountsMustAgree`: `views` com n = 3 e `data` com n = 2 dão `"views have 3 assets; data has 2"`.
  - `Omega.TauMustBePositive`: τ = 0 dá `"tau must be positive and finite"`.
- [ ] **Step 2: Rodar e ver falhar.** `ctest --preset gcc-release -R "Equilibrium|Views|Omega"`.
- [ ] **Step 3: Implementar.** O posto vem do Householder QR de `P'` (`n x k`): transponha P para uma matriz temporária e chame `householder_qr(pt, n, k)` e `has_full_column_rank`. `p_k Sigma p_k'` usa Σ de `data.covariance()`.
- [ ] **Step 4: Rodar e ver passar.** O `ctest` inteiro.
- [ ] **Step 5: Commit.** `feat(estimation): Black-Litterman views, equilibrium and omega`

---

### Task 5: O posterior do Black-Litterman

**Files:**
- Modify: `src/estimation/black_litterman.hpp`, `src/estimation/black_litterman.cpp`
- Test: `tests/unit/estimation/black_litterman_test.cpp`

**Interfaces:**
- Produces: `struct po::estimation::Posterior { std::vector<double> mu; std::vector<double> covariance; };` (`covariance` vazio quando não é pedida)
- Produces: `std::expected<Posterior, Failure> posterior(const Data& data, const Views& views, double tau, bool with_covariance);`

Fixture: o de Task 4, com `Pi = {0.07, 0.105}` no `mu` de `data` e `Q = {0.02}`.

- [ ] **Step 1: Testes que falham.**
  - `BlackLitterman.HeLittermanByHand`: Ω = `{0.0055}` dá `mu = {0.0775, 0.085}` (`1e-16`). `Sigma_post` sai da forma escrita no teste: `1.05 Sigma - B B' / 0.011`, com `B = {0.0015, -0.004}`.
  - `BlackLitterman.NoViewsReturnsThePrior`: k = 0 dá `mu == Pi` bit a bit e `Sigma_post = 1.05 Sigma`.
  - `BlackLitterman.FullConfidenceReproducesTheView`: Ω = `{0}` dá `|mu[0] - mu[1] - 0.02| <= 1e-15`. Também com n = 4, k = 2 e Ω = 0.
  - `BlackLitterman.PosteriorCovarianceMatchesTheInverseForm`: com n = 3, k = 2 e Ω diagonal positiva, `Sigma_post` fica igual a `Sigma + ((tau Sigma)^-1 + P' Omega^-1 P)^-1` (inversas por Gauss-Jordan no teste), com `1e-14` relativo. A saída é exatamente simétrica.
  - `BlackLitterman.OmegaIsRequired`: k = 1 sem `set_omega` dá `"omega is not set"`.
  - `BlackLitterman.NearlyRedundantViewsAreFiniteOrRejected`: `P = {1, -1, 1, -1 + 1e-9}` com Ω = 0. Ou todas as saídas são finitas, ou o erro é `"P tau Sigma P' + omega is not positive definite: pivot 1 is not positive"`. Nunca NaN.
  - `BlackLitterman.WithoutCovarianceLeavesItEmpty`: `with_covariance = false` dá `covariance.empty()` e o mesmo `mu`.
- [ ] **Step 2: Rodar e ver falhar.** `ctest --preset gcc-release -R BlackLitterman`.
- [ ] **Step 3: Implementar.** A conta segue spec §6.6: `B = tau Sigma P'`, `A = P B + Omega`, o fator `L` por `cholesky_or_pivot` (a mensagem da Task 1, com o nome `P tau Sigma P' + omega`), `mu = Pi + B (L L')^-1 (Q - P Pi)` e, com covariância, `X = L^-1 B'` (k × n) e `Sigma_post[i][j] = (1 + tau) Sigma[i][j] - sum_r X[r][i] X[r][j]`, só para `j >= i`, espelhado. Para manter as funções em até 60 linhas, cada termo vai numa função própria. O arquivo fica abaixo de cerca de 300 linhas. Se passar disso, divida em `black_litterman` (equilíbrio e Ω) e `posterior`.
- [ ] **Step 4: Rodar e ver passar.** O `ctest` inteiro, depois `native-release`.
- [ ] **Step 5: Commit.** `feat(estimation): Black-Litterman posterior mean and covariance`

---

### Task 6: A interface C

**Files:**
- Modify: `include/portopt/portopt.h` (seção "Estimation utilities" depois das métricas; os comentários de contrato copiam spec §6, a nota de frequência do Bayes-Stein e a de excesso do rf inclusive)
- Modify: `include/portopt/portopt.abi` (`typedef struct po_views po_views` e as 9 linhas `function`, na ordem do header)
- Modify: `src/capi/handles.hpp` (`struct po_views { po::estimation::Views views; };`)
- Create: `src/capi/estimation/views_api.cpp`, `src/capi/estimation/estimation_api.cpp`
- Create: `tests/unit/capi/estimation/estimation_test.cpp`
- Modify: `CMakeLists.txt` (`portopt_capi` e `portopt_tests`)

**Interfaces:**
- Consumes: tudo das Tasks 2 a 5.
- Produces: as 9 funções de spec §6, com as assinaturas exatas de lá: `po_sample_moments`, `po_bayes_stein`, `po_bl_equilibrium`, `po_views_create`, `po_views_set_omega`, `po_views_free`, `po_bl_omega_he_litterman`, `po_bl_omega_idzorek` e `po_bl_posterior`.

- [ ] **Step 1: Testes que falham.** O padrão é o de `tests/unit/capi/metrics_test.cpp`, com um `ViewsPtr` local (deleter `po_views_free`) no arquivo de teste.
  - `CapiEstimation.EachFunctionMatchesTheCore`: uma chamada de cada função sobre os fixtures das Tasks 2 a 5, com os mesmos valores.
  - `CapiEstimation.NullHandles`: `d` ou `v` nulos dão `PO_INVALID_HANDLE`, e a mensagem é `"d must not be null"` ou `"v must not be null"`.
  - `CapiEstimation.NullArguments`: cada ponteiro obrigatório nulo dá `PO_INVALID_ARGUMENT`. `out_cov` nulo no posterior dá `PO_OK`.
  - `CapiEstimation.OutputUntouchedOnFailure`: as saídas preenchidas com `-7.0` continuam `-7.0` depois de uma falha de domínio, em cada função.
  - `CapiEstimation.FailedCreateLeavesNull`: `po_views_create` com P redundante deixa `*out == nullptr`. `po_views_free(nullptr)` não faz nada.
  - `CapiEstimation.PivotMessage`: uma Σ com o pivô 2 não positivo, passada ao `po_data_create`, dá a mensagem com `pivot 2`.
- [ ] **Step 2: Rodar e ver falhar.** `ctest --preset gcc-release -R CapiEstimation`.
- [ ] **Step 3: Implementar.** Com `guard`, `null_handle`, `null_argument` e `from_failure` de `capi/handles.hpp`. As saídas são copiadas só depois do sucesso da função do core.
- [ ] **Step 4: Rodar e ver passar.** O `ctest` inteiro: `abi_contract`, `style` e `portopt_c_header_check` inclusive. Depois `native-release`. Rode também o `clang-tidy` do CI sobre os arquivos novos (`git ls-files 'src/*.cpp' | xargs clang-tidy-20 -p build/tidy`, `.github/workflows/ci.yml:182`), se estiver disponível. Se não estiver, anote que o CI confere.
- [ ] **Step 5: Commit.** `feat(capi): the estimation utilities and po_views`

---

### Task 7: Os bits do determinismo

**Files:**
- Create: `tests/unit/estimation/determinism_test.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: as funções do core das Tasks 2 a 5.

- [ ] **Step 1: O teste, com valores vazios.** Siga o cabeçalho e o formato de `tests/unit/solve/determinism_test.cpp`. O caso fixo:
  - um painel literal de 8 × 3 em `sample_moments`;
  - `bayes_stein` sobre ele, com t = 8 e n = 3 (`t > n + 2`);
  - `equilibrium` com `w = {0.5, 0.3, 0.2}` e δ = 2.5;
  - `omega_idzorek` com `P = {1, -1, 0, 0, 1, -1}`, τ = 0.05 e `c = {0.6, 0.9}`;
  - `posterior`.

  Cada saída é comparada bit a bit (`std::bit_cast<std::uint64_t>`) com hex-floats registrados no arquivo.
- [ ] **Step 2: Rodar, capturar e registrar.** Rode em `gcc-release`, copie os valores impressos na falha para o arquivo e rode de novo. Esperado: PASS.
- [ ] **Step 3: Conferir em MSVC.** `ctest --preset native-release -R EstimationDeterminism`. Esperado: PASS sem mudar valor nenhum. Se houver divergência, é defeito: pare e leve ao autor, sem afrouxar.
- [ ] **Step 4: Commit.** `test(estimation): the committed bits of the estimation utilities`

---

### Task 8: O wrapper C++

**Files:**
- Modify: `cpp/include/portopt/portopt.hpp`
- Test: `cpp/tests/wrapper_test.cpp`

**Interfaces:**
- Consumes: as 9 funções C (Task 6).
- Produces, em `namespace po`:
  - `class Views` (move-only, como `Data`), com `static Expected<Views> create(std::size_t assets, std::span<const double> p, std::span<const double> q)`, em que k é `q.size()` e `p.size()` tem de ser `k * assets`; com `Expected<void> set_omega(std::span<const double> omega)`; e com `const po_views* get() const noexcept`;
  - `struct Moments { std::vector<double> mu; std::vector<double> cov; };`, `struct BayesSteinEstimate { std::vector<double> mu; double grand_mean; double weight; };` e `struct Posterior { std::vector<double> mu; std::vector<double> cov; };`;
  - `Expected<Moments> sample_moments(std::span<const double> returns, std::size_t assets)`;
  - `Expected<BayesSteinEstimate> bayes_stein(const Data& data, std::size_t periods)`;
  - `Expected<std::vector<double>> bl_equilibrium(const Data& data, std::span<const double> w_mkt, double delta)`;
  - `Expected<std::vector<double>> bl_omega_he_litterman(const Data& data, const Views& views, double tau)`;
  - `Expected<std::vector<double>> bl_omega_idzorek(const Data& data, const Views& views, double tau, std::span<const double> confidence)`;
  - `Expected<Posterior> bl_posterior(const Data& data, const Views& views, double tau)`.

- [ ] **Step 1: Testes que falham.** `Wrapper.EstimationRoundTrip` vai do painel até um `solve` de `MAX_SHARPE` sobre o posterior, com pesos que somam 1. `Wrapper.ViewsSizeMismatch` verifica que `p.size() != q.size() * assets` dá `Error` antes da chamada nativa.
- [ ] **Step 2: Rodar e ver falhar.** `ctest --preset gcc-release -R Wrapper`.
- [ ] **Step 3: Implementar.** Siga `Data`, `detail::check` e `detail::check_size`.
- [ ] **Step 4: Rodar e ver passar.** `ctest --preset gcc-release` e `python cpp/scripts/check_abi_coverage.py`. Esperado: `ok: the wrapper uses N of N PO_API functions`.
- [ ] **Step 5: Commit.** `feat(cpp): estimation utilities in the C++ wrapper`

---

### Task 9: O binding Python

**Files:**
- Modify: `python/src/portopt/_native.py` (as 9 entradas de `PROTOTYPES`)
- Create: `python/src/portopt/_estimation.py`
- Modify: `python/src/portopt/__init__.py` (exporta e acrescenta a `__all__`)
- Create: `python/tests/test_estimation.py`

**Interfaces:**
- Consumes: as 9 funções C; `Data` de `_data.py` (o `pointer`); `_native.create_handle`, `check` e `as_doubles`/`doubles_pointer` de `_spec.py`.
- Produces:
  - `sample_moments(returns: ArrayLike) -> tuple[NDArray, NDArray]`;
  - `class BayesStein(NamedTuple): mu: NDArray; grand_mean: float; weight: float`;
  - `bayes_stein(mu_hat: ArrayLike, cov: ArrayLike, t: int) -> BayesStein`;
  - `bl_equilibrium(cov: ArrayLike, w_mkt: ArrayLike, delta: float) -> NDArray`;
  - `@dataclass(frozen=True) class Views: p: NDArray; q: NDArray; omega: NDArray | None = None`, que valida as formas em `__post_init__` (p 2-D `(k, n)`; q 1-D com k entradas; omega `(k, k)` ou None só com `k = 0`);
  - `bl_omega_he_litterman(cov: ArrayLike, p: ArrayLike, tau: float) -> NDArray` e `bl_omega_idzorek(cov: ArrayLike, p: ArrayLike, tau: float, confidence: ArrayLike) -> NDArray`;
  - `class BlackLittermanPosterior(NamedTuple): mu: NDArray; cov: NDArray`, com `def to_data(self, risk_free: float | None = None) -> Data`;
  - `bl_posterior(cov: ArrayLike, pi: ArrayLike, tau: float, views: Views) -> BlackLittermanPosterior`.

  As funções de Ω recebem `p` solto porque Ω ainda não existe quando elas rodam. Elas criam um `po_views` com `q = 0`.

- [ ] **Step 1: Testes que falham.** As referências são numpy, escritas no teste com `np.linalg.solve`.
  - `test_sample_moments_match_numpy`: `np.mean(axis=0)` e `np.cov(rowvar=False, ddof=1)`, com `rtol=1e-14`.
  - `test_bayes_stein_matches_numpy`: as fórmulas de spec §6.2 em numpy.
  - `test_relative_view_matches_numpy`: n = 3, P = `[[1, -1, 0]]`, Ω de He-Litterman e o posterior contra numpy. É a referência pedida pelo DataPanel.
  - `test_full_confidence_reproduces_the_view`: Idzorek com `c = 1` dá `P @ mu ≈ Q` (`atol=1e-15`).
  - `test_no_views_returns_the_prior`: `Views(np.empty((0, 3)), np.empty(0))` dá `mu == pi`.
  - `test_to_data_adds_the_rate`: `to_data(0.02)` dá um `Data` com `weight_count == n + 1`, e `portfolio_return` do vetor com só o primeiro ativo é `mu[0] + 0.02`. `to_data()` dá `weight_count == n`.
  - `test_views_reject_one_dimensional_p`: `Views(np.array([1.0, -1.0]), np.array([0.02]), ...)` dá `ValueError` com `"p must be 2-D"`.
  - `test_redundant_views_raise`: dá `ValueError`, com a mensagem nativa `"views are redundant"`.
- [ ] **Step 2: Rodar e ver falhar.** `PORTOPT_LIB=... uv run --project python --python 3.14 pytest python/tests/test_estimation.py -q`.
- [ ] **Step 3: Implementar.** Use os padrões de `_data.py`: formas conferidas antes da chamada nativa, `Data(mu, cov)` para o `po_data` e `create_handle(lib.po_views_create, lib.po_views_free, k, n, p, q)`. Cada função tem no máximo 5 parâmetros, e não há `else` nem `elif`.
- [ ] **Step 4: Rodar e ver passar.** O pytest inteiro em 3.14 e em 3.11, o `mypy` do CI e o `ctest -R style`.
- [ ] **Step 5: Commit.** `feat(python): estimation utilities`

---

### Task 10: O binding C#

**Files:**
- Modify: `csharp/src/PortOpt/Interop/Native/NativeMethods.cs` (as 9 `LibraryImport`)
- Modify: `csharp/src/PortOpt/Interop/Handles.cs` (`ViewsHandle`)
- Modify: `csharp/src/PortOpt/Market/MarketData.cs` (`internal DataHandle Handle => _handle;`, se ainda não houver um acesso interno)
- Create: `csharp/src/PortOpt/Estimation/Estimation.cs`, `csharp/src/PortOpt/Estimation/Records.cs`
- Modify: `csharp/src/PortOpt/PublicAPI.Unshipped.txt`
- Create: `csharp/tests/PortOpt.Tests/EstimationTests.cs`

**Interfaces:**
- Consumes: as 9 funções C; `MarketData` e `MarketInputs`; `Calls.Check`.
- Produces, em `namespace PortOpt`:
  - `public static class Estimation`, com:
    - `Moments SampleMoments(double[,] returns)`;
    - `BayesSteinEstimate BayesStein(double[] muHat, double[,] covariance, int periods)`;
    - `double[] Equilibrium(double[,] covariance, double[] marketWeights, double delta)`;
    - `double[,] OmegaHeLitterman(double[,] covariance, double[,] p, double tau)`;
    - `double[,] OmegaIdzorek(double[,] covariance, double[,] p, double tau, double[] confidence)`;
    - `BlackLittermanPosterior Posterior(double[,] covariance, double[] pi, double tau, Views views)`;
  - `public sealed record Moments(double[] Mu, double[,] Covariance);`
  - `public sealed record BayesSteinEstimate(double[] Mu, double GrandMean, double Weight);`
  - `public sealed record Views(double[,] P, double[] Q, double[,]? Omega);`
  - `public sealed record BlackLittermanPosterior(double[] Mu, double[,] Covariance)`, com `public MarketInputs ToMarketInputs(double? riskFree = null)`, que copia os arrays (o `MarketInputs` não pode partilhar o array do record).

- [ ] **Step 1: Testes que falham.** Os mesmos casos e valores de Task 9, com tolerância `1e-14`. As referências numéricas são as literais do Python, copiadas.
  - `Posterior_reproduces_a_full_confidence_view`.
  - `ToMarketInputs_adds_the_rate_and_sets_cash`: com rf, `Mu[i] == posterior.Mu[i] + rf` e `RiskFree == rf`; sem rf, `RiskFree` é null.
  - `Redundant_views_throw_PortOptException`: `Status == Status.InvalidArgument`.
  - `Views_with_a_wrong_shape_throw_ArgumentException`.
  - Os testes de contrato e de cobertura existentes (`ContractTests`, `CoverageTests`) passam sem mudança.
- [ ] **Step 2: Rodar e ver falhar.** `PORTOPT_LIB=... dotnet test --solution csharp/PortOpt.slnx -c Release`.
- [ ] **Step 3: Implementar.** Siga `MarketData` (`unsafe`, `fixed`, `Calls.Check`). Cada método cria um `MarketData` e, no BL, um `ViewsHandle`, e descarta os dois com `using`. As entradas de `PublicAPI.Unshipped.txt` saem do analisador (o build falha até que estejam lá).
- [ ] **Step 4: Rodar e ver passar.** `dotnet test` inteiro e `csharp/scripts/check_package.sh`, se o ambiente permitir (o consumer com Native AOT). Se não permitir, anote que o CI confere.
- [ ] **Step 5: Commit.** `feat(csharp): estimation utilities and the posterior to MarketInputs`

---

### Task 11: Documentação

**Files:**
- Modify: `docs/guide/math.md` (seção "Estimation utilities", com as fórmulas de spec §6, a nota de frequência, a de excesso do rf para `Pi` e `Q` e as referências: Black e Litterman, 1992; He e Litterman, 1999; Idzorek, 2005; Walters, 2014; Jorion, 1986)
- Modify: `docs/guide/validation.md` (as regras de domínio e as tolerâncias: `1e-10` de simetria, `1e-8` de soma dos pesos de mercado, a regra de pivô e o posto pelo QR)
- Modify: `ARCHITECTURE.md` (a camada `src/estimation`, `src/capi/estimation`, a tabela de responsabilidades e a frase de abertura com as utilidades)
- Modify: `CHANGELOG.md` (`## [Unreleased]` → `### Added`, e `### Changed` para a mensagem do pivô)
- Modify: `python/README.md`, `csharp/src/PortOpt/README.md` e `cpp/README.md` (um exemplo: do painel ao posterior com uma visão relativa e daí ao `solve`)
- Modify: `docs/superpowers/specs/2026-11-15-portopt-design.md:44` (uma nota de que o não-objetivo de estimar μ ou Σ foi revisto pelo spec de 2026-11-25)

- [ ] **Step 1: Escrever.** A documentação do repositório é em inglês, e a nota no spec antigo é em português.
- [ ] **Step 2: Conferir.** Rode cada exemplo dos READMEs: o Python com `PORTOPT_LIB`, e o C# e o C++ compilando um consumer mínimo ou um teste temporário que você não commita. Também `ctest --preset gcc-release -R style`.
- [ ] **Step 3: Commit.** `docs: estimation utilities`
