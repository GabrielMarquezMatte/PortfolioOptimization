# PortOpt — Plano 1: fundação e MIN_VARIANCE de ponta a ponta

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Um objetivo (`MIN_VARIANCE`, com bounds, `target_return` e caixa) resolvido de ponta a ponta — núcleo C++23 → interface C → Python — e validado por KKT e formas fechadas, com o contrato (enums completos, status, regras de erro) já na forma final.

**Architecture:** Núcleo C++23 interno em `src/`, exposto só por `include/portopt/portopt.h` (handles opacos, setters, `po_status`, sem structs públicas). Um QP denso de Goldfarb-Idnani próprio resolve tudo; `src/problems/` traduz spec + dados para o QP e a solução de volta para as variáveis do usuário. Python via ctypes.

**Tech Stack:** C++23 (GCC ≥ 15, MSVC VS 18, Clang), CMake ≥ 3.28 + Ninja + `CMakePresets.json`, FetchContent, GoogleTest v1.17.0, Python ≥ 3.11 + numpy + pytest + mypy, uv.

**Spec:** `docs/superpowers/specs/2026-11-15-portopt-design.md`

## Roteiro dos planos

Cada plano entrega software testável e é escrito quando o anterior termina.

1. **Fundação + `MIN_VARIANCE` de ponta a ponta** (este).
2. Objetivos de QP: `UTILITY`, `IMPORTANCE`, `KELLY`, `TRACKING_ERROR` (dados de retornos), `MAX_RETURN` (busca em δ); métricas `utility` e `tracking_error`.
3. Homogeneização e leverage: `MAX_SHARPE`, `MAX_DIVERSIFICATION`, status `UNBOUNDED`, família implícita do leverage em todos os objetivos de QP; métricas `sharpe_ratio` e `diversification_ratio`.
4. `RISK_PARITY`: Newton de Spinu, Lagrangiano aumentado, resíduo; métrica `risk_contributions`.
5. `cpp/portopt.hpp`, benchmarks, wheels/PyPI, release, preparação da 1.0.

## Global Constraints

- A interface C é o único contrato estável. Nenhuma struct pública; configuração só por setters sobre handles opacos. Configuração nova é função nova, nunca parâmetro novo nem `_ex`.
- Prefixo `po_` para funções e tipos, `PO_` para macros; namespace C++ `po`. Tipos na fronteira: `int32_t`, `size_t`, `double`, ponteiros.
- Enums da interface C são `typedef int32_t` + `#define`; valores só crescem por acréscimo.
- Toda função falível devolve `po_status`; mensagem em `po_last_error()`, thread-local, `""` depois de sucesso.
- Nenhuma exceção atravessa a interface C (`std::bad_alloc` → `PO_OUT_OF_MEMORY`, resto → `PO_INTERNAL_ERROR`), por um único `guard`.
- Sem estado global mutável. Handles independentes; spec e dados copiados no `po_solve`; `po_result` imutável.
- Nenhuma gambiarra numérica (variância `eps` do caixa, grid search, bounds ignorados).
- Só `double`. Nenhuma dependência além de GoogleTest (testes).
- Anti-objetivo: não fazer mais do que o necessário. Nada fora do que este plano lista.
- O código do `quadprog` (GPL) não é copiado nem traduzido; o Goldfarb-Idnani vem do paper (Goldfarb & Idnani, *Math. Programming* 27, 1983).
- Commits sem trailer `Co-Authored-By`.
- Estilo: `STYLEGUIDE.md` do VolEngine, integralmente, com os nomes deste projeto (criado na Task 1). Em resumo: nunca `else`; `if` sempre com chaves e corpo em linha própria; ternário só com nomes/literais/membros; aninhamento ≤ 3; função ≤ 60 linhas; ≤ 5 parâmetros; ≤ 8 módulos por pasta; Allman, 4 espaços; tipos `PascalCase`, funções e variáveis `snake_case`, membros privados com sufixo `_`, constantes `kPascalCase`; comentários só para restrição externa ou invariante; medições em `docs/`; nunca inversa explícita (Cholesky/solve); covariâncias simetrizadas.
- Código, identificadores, mensagens de erro e documentos do repositório em inglês; documentos em `docs/superpowers/` em português.
- Warnings são erros (`-Wall -Wextra -Wpedantic -Werror`, `-Wnrvo` no GCC 14+; `/W4 /WX` no MSVC).

## Review Focus

1. **Σ mal condicionada** (autovalor relativo 1e-10, mas definida positiva): `po_data_create` aceita e o solve dá `OPTIMAL` com KKT dentro da tolerância escalada. Teste na Task 4 (`MinVariance.IllConditionedCovariance`).
2. **Igualdade redundante e consistente** (todos os μ iguais e `target_return` igual a esse μ, sem caixa): `OPTIMAL`, não `INFEASIBLE`. Testes na Task 3 (`GoldfarbIdnani.RedundantConsistentEquality`) e na Task 4 (`MinVariance.TargetReturnEqualToCommonMean`).
3. **Bounds que fixam um ativo** (lo == hi): o ativo sai exatamente no valor e os outros otimizam. Teste na Task 4 (`MinVariance.FixedAsset`).
4. **Arrays numpy não contíguos, Fortran-order ou inteiros**: o wrapper converte com `np.ascontiguousarray(..., dtype=np.float64)`. Teste na Task 6 (`test_data_accepts_fortran_and_int_arrays`).
5. **Um único ativo**: `w = [1]`; `target_return` diferente de μ₁ dá `INFEASIBLE`. Teste na Task 4 (`MinVariance.SingleAsset`).

---

### Task 1: Esqueleto do repositório, styleguide, header C e erro

**Files:**
- Create: `.gitignore`, `.clang-format`, `.clang-tidy`, `.editorconfig` (copiados do VolEngine)
- Create: `STYLEGUIDE.md` (o do VolEngine com `ve`/`VE`/`volengine` → `po`/`PO`/`portopt`; a seção "Layers" reescrita para as camadas da Task 7; exemplos de domínio GARCH trocados por exemplos de carteira onde aparecem)
- Create: `LICENSE` (MIT, "Copyright (c) 2026 Gabriel Matte")
- Create: `CMakeLists.txt`, `CMakePresets.json`
- Create: `include/portopt/portopt.h`
- Create: `src/capi/error.hpp`, `src/capi/error.cpp`, `src/capi/version_api.cpp`
- Create: `tests/style/check_style.py` (copiado do VolEngine; `IGNORED_DIRECTORIES = ("build", "python/.venv")`)
- Test: `tests/unit/capi/abi_test.cpp`, `tests/unit/capi/c_header_check.c`

**Interfaces:**
- Produces: alvos CMake `portopt_core` (STATIC, PIC, `target_include_directories(PUBLIC src include)`), `portopt` (SHARED, define `PO_BUILDING`), `portopt_tests` (GoogleTest, `tests/unit` no include path), teste ctest `style`. Opções `PO_BUILD_TESTS` (ON quando top-level), `PO_BUILD_BENCHMARKS` (OFF, sem alvo ainda). `project(PortOpt VERSION ${PO_VERSION})` com `PO_VERSION` padrão `0.1.0`.
- Produces (header):
  ```c
  #define PO_API  /* __declspec(dllexport/dllimport) no Windows por PO_BUILDING; visibility("default") fora */
  #define PO_ABI_VERSION 0
  PO_API int32_t po_abi_version(void);
  PO_API const char* po_version(void);          /* "0.1.0", de PO_VERSION */
  typedef int32_t po_status;
  #define PO_OK 0
  #define PO_INVALID_ARGUMENT 1
  #define PO_INVALID_HANDLE 2
  #define PO_UNSUPPORTED 3
  #define PO_OUT_OF_MEMORY 4
  #define PO_INTERNAL_ERROR 5
  PO_API const char* po_last_error(void);
  ```
- Produces (C++): `po::capi::set_error(po_status, std::string_view) -> po_status`, `clear_error()`, `last_error() -> const char*`, `template <typename F> po_status guard(F&&)` — mesma forma de `VolEngine/src/capi/error.hpp`, sem o `no_standard_errors`.

- [ ] **Step 1:** `git init`, branch `main`. Copiar `.gitignore`, `.clang-format`, `.clang-tidy`, `.editorconfig` e `tests/style/check_style.py` do VolEngine.
- [ ] **Step 2: Escrever os testes**

  `abi_test.cpp`:
  ```cpp
  TEST(Abi, VersionMatchesHeader) { EXPECT_EQ(po_abi_version(), PO_ABI_VERSION); }
  TEST(Abi, VersionStringIsProjectVersion) { EXPECT_STREQ(po_version(), "0.1.0"); }
  TEST(Abi, LastErrorStartsEmpty) { EXPECT_STREQ(po_last_error(), ""); }
  TEST(Guard, BadAllocBecomesOutOfMemory) {
      EXPECT_EQ(po::capi::guard([]() -> po_status { throw std::bad_alloc(); }), PO_OUT_OF_MEMORY);
      EXPECT_STREQ(po_last_error(), "out of memory");
  }
  TEST(Guard, OtherExceptionBecomesInternalError) {
      EXPECT_EQ(po::capi::guard([]() -> po_status { throw std::runtime_error("boom"); }), PO_INTERNAL_ERROR);
      EXPECT_STREQ(po_last_error(), "boom");
  }
  TEST(Guard, SuccessClearsLastError) {
      po::capi::set_error(PO_INVALID_ARGUMENT, "x");
      EXPECT_EQ(po::capi::guard([] { return PO_OK; }), PO_OK);
      EXPECT_STREQ(po_last_error(), "");
  }
  ```
  `c_header_check.c`: `#include <portopt/portopt.h>` e uma função que chama `po_abi_version()`; compilado como C11 num alvo OBJECT `portopt_c_header_check`, como no VolEngine.
- [ ] **Step 3:** Escrever `CMakeLists.txt` (forma do VolEngine sem Boost, Highway e thread pool; mesmas flags de warning e o mesmo bloco de link estático MinGW/Linux, com `portopt.def`) e `CMakePresets.json` (os presets `gcc-debug`, `gcc-release`, `native-release`, `msvc-release` do VolEngine; `gcc-bench` e `msvc-profile` ficam para o plano 5). `po_version()` lê `PO_VERSION` por `target_compile_definitions(portopt PRIVATE PO_VERSION_STRING="${PROJECT_VERSION}")`.
- [ ] **Step 4:** Escrever o header, `error.*` e `version_api.cpp` (cada `extern "C"` dentro de `guard`, exceto `po_abi_version`, `po_version` e `po_last_error`, que não falham).
- [ ] **Step 5:** Escrever `STYLEGUIDE.md` e `LICENSE`.
- [ ] **Step 6: Verificar**

  Run: `cmake --preset gcc-release && cmake --build --preset gcc-release && ctest --preset gcc-release`
  Expected: todos os testes `Abi.*` e `Guard.*` PASS, `style` PASS.
- [ ] **Step 7: Commit** — `git add -A && git commit -m "build: repository skeleton, C header, error boundary"`

---

### Task 2: Cholesky e matrizes de teste

**Files:**
- Create: `src/linalg/cholesky.hpp`, `src/linalg/cholesky.cpp`
- Create: `tests/unit/random_covariance.hpp`
- Test: `tests/unit/linalg/cholesky_test.cpp`

**Interfaces:**
- Produces:
  ```cpp
  namespace po::linalg {
  // Lower factor L (row-major n x n, upper triangle zero) with A = L L'; nullopt when A is not positive definite.
  std::optional<std::vector<double>> cholesky(std::span<const double> a, std::size_t n);
  void solve_lower(std::span<const double> l, std::size_t n, std::span<double> b);        // L x = b, in place
  void solve_lower_transposed(std::span<const double> l, std::size_t n, std::span<double> b); // L' x = b
  }
  namespace po::testing {
  // Seeded SPD matrix: Q diag(eigenvalues) Q', Q from Householder of a seeded normal matrix; row-major.
  std::vector<double> random_covariance(std::size_t n, std::uint64_t seed, double condition = 100.0);
  std::vector<double> random_vector(std::size_t n, std::uint64_t seed, double low, double high);
  }
  ```
  `random_covariance` usa `std::mt19937_64` e autovalores log-espaçados de 1 a `condition`, escalados por 0.04 (vol ≈ 20% na unidade do teste). Só `<random>` com distribuição normal própria (Box-Muller) para ser idêntico entre compiladores.

- [ ] **Step 1: Testes**
  - `Cholesky.ReconstructsMatrix`: `random_covariance(5, 1)`, `L L'` igual a A em 1e-14 relativo.
  - `Cholesky.RejectsIndefinite`: `{1, 2, 2, 1}` → `nullopt`.
  - `Cholesky.SolvesSystem`: com A de `random_covariance(6, 2)` e b de `random_vector`, `solve_lower` seguido de `solve_lower_transposed` dá x com ‖Ax − b‖∞ ≤ 1e-12·‖b‖∞.
  - `RandomCovariance.IsDeterministicAndSymmetric`: duas chamadas com a mesma semente dão vetores idênticos; A simétrica exatamente.
- [ ] **Step 2:** Rodar, ver falhar (não compila).
- [ ] **Step 3:** Implementar.
- [ ] **Step 4:** `ctest --preset gcc-release` PASS.
- [ ] **Step 5: Commit** — `feat(linalg): Cholesky factor and triangular solves`

---

### Task 3: QP de Goldfarb-Idnani

**Files:**
- Create: `src/qp/problem.hpp`, `src/qp/goldfarb_idnani.hpp`, `src/qp/goldfarb_idnani.cpp`
- Create: `tests/unit/kkt.hpp`
- Test: `tests/unit/qp/goldfarb_idnani_test.cpp`

**Interfaces:**
- Consumes: `po::linalg::cholesky`, `solve_lower*` (Task 2).
- Produces:
  ```cpp
  namespace po::qp {
  // min 0.5 x'Gx + g'x  s.t.  E x = e (rows of E),  C x >= c (rows of C).
  struct Problem
  {
      std::size_t n = 0;
      std::vector<double> hessian;          // G, n x n row-major, positive definite
      std::vector<double> linear;           // g, n
      std::vector<double> equality_rows;    // E, m_eq x n row-major
      std::vector<double> equality_rhs;     // e
      std::vector<double> inequality_rows;  // C, m_in x n row-major
      std::vector<double> inequality_rhs;   // c
  };
  enum class Status { Optimal, Infeasible, MaxIterations };
  // Stationarity: G x + g = E' y + C' z, with z >= 0 and z_i (C_i x - c_i) = 0.
  struct Solution
  {
      Status status = Status::Infeasible;
      std::vector<double> x;                // NaN when not Optimal
      std::vector<double> equality_multipliers;    // y
      std::vector<double> inequality_multipliers;  // z
      int iterations = 0;
  };
  // Hessian not positive definite is a precondition violation: throws std::invalid_argument.
  Solution solve(const Problem& problem);
  }
  namespace po::testing {
  struct KktResiduals { double stationarity; double primal; double dual; double complementarity; };
  KktResiduals kkt_residuals(const po::qp::Problem& problem, const po::qp::Solution& solution);
  // Each residual relative to max(1, |G|max * |x|max, |g|max) (stationarity) or max(1, |rhs|max) (the others).
  void expect_kkt(const po::qp::Problem& problem, const po::qp::Solution& solution, double tolerance = 1e-10);
  }
  ```

  Algoritmo (dual active-set do paper; o implementador segue o paper, com estas escolhas fixadas):
  1. Fatorar G = LL'; x = −G⁻¹g; J = L⁻ᵀ; conjunto ativo vazio; R vazio (triangular superior q×q).
  2. Escolher a restrição a adicionar: primeiro cada igualdade ainda não ativa, na ordem (tratada como `n'x ≥ b` se `n'x − b < 0`, como `−n'x ≥ −b` caso contrário; violação nula → marcar ativa sem passo); depois a desigualdade mais violada com `C_i x − c_i < −tol`. Nenhuma → `Optimal`.
  3. d = J'n⁺; z = J₂d₂; r = R⁻¹d₁.
  4. t₁ = min uⱼ/rⱼ sobre desigualdades ativas com rⱼ > 0 (igualdades nunca saem); t₂ = −s(x)/(z'n⁺) se ‖z‖ > tol, senão ∞. t = min(t₁, t₂).
  5. t = ∞ → `Infeasible`. t₂ = ∞ → passo só dual, remover a restrição k, voltar ao 3. Senão x += t z, u += t(−r, 1); t = t₂ → adicionar (atualizar J e R por Givens) e voltar ao 2; t = t₁ → remover k (Givens) e voltar ao 3.
  6. Limite de 10·(n + m_eq + m_in) + 10 iterações → `MaxIterations`.
  `tol` = 1e-12 relativo a max(1, |rhs da restrição|, ‖linha‖·‖x‖∞). Restrição com linha nula e rhs satisfeito é ignorada; com rhs violado, `Infeasible`.

- [ ] **Step 1: Testes** (`G` e `g` de `random_covariance` / `random_vector`; cada teste `Optimal` chama `expect_kkt`)
  - `GoldfarbIdnani.Unconstrained`: x = −G⁻¹g em 1e-12 relativo.
  - `GoldfarbIdnani.EqualityOnly`: n = 4, E = 1', e = 1, g = 0 → x = G⁻¹1 / 1'G⁻¹1 em 1e-12; y > 0.
  - `GoldfarbIdnani.ActiveBounds`: n = 5, x ≥ 0 com g escolhido para que o irrestrito tenha 2 negativos → esses saem 0, z correspondentes > 0, KKT.
  - `GoldfarbIdnani.InfeasibleBounds`: n = 3, 1'x = 1, x ≤ 0.3 → `Infeasible`, x todo NaN.
  - `GoldfarbIdnani.RedundantConsistentEquality`: E = [1'; 2·1'], e = [1, 2] → `Optimal`, KKT.
  - `GoldfarbIdnani.InconsistentEqualities`: E = [1'; 1'], e = [1, 2] → `Infeasible`.
  - `GoldfarbIdnani.EqualBoundsFixVariable`: x₀ ≥ 0.2 e −x₀ ≥ −0.2 → x₀ = 0.2 em 1e-14.
  - `GoldfarbIdnani.RejectsIndefiniteHessian`: `EXPECT_THROW(solve(p), std::invalid_argument)`.
  - `GoldfarbIdnani.RandomProblemsSatisfyKkt`: sementes 1..50, n em 2..30, 1'x = 1 e bounds aleatórios lo ∈ [−0.5, 1/n], hi ∈ [1/n, 1] (o ponto 1/n é sempre viável) → todos `Optimal` com KKT.
- [ ] **Step 2:** Rodar, ver falhar.
- [ ] **Step 3:** Implementar `solve` e o helper KKT.
- [ ] **Step 4:** PASS em `gcc-release` e `gcc-debug`.
- [ ] **Step 5: Commit** — `feat(qp): dense Goldfarb-Idnani dual active-set solver`

---

### Task 4: Spec, dados e `MIN_VARIANCE` no núcleo

**Files:**
- Create: `src/core/objective.hpp`, `src/core/spec.hpp`, `src/core/data.hpp`, `src/core/data.cpp`, `src/core/result.hpp`, `src/core/failure.hpp`
- Create: `src/problems/constraints.hpp`, `src/problems/constraints.cpp`, `src/problems/min_variance.hpp`, `src/problems/min_variance.cpp`, `src/problems/solve.hpp`, `src/problems/solve.cpp`
- Test: `tests/unit/core/data_test.cpp`, `tests/unit/problems/min_variance_test.cpp`, `tests/unit/problems/solve_test.cpp`

**Interfaces:**
- Consumes: `po::qp::solve`, `Problem`, `Solution` (Task 3); `po::linalg::cholesky` (Task 2).
- Produces:
  ```cpp
  namespace po {
  enum class Objective : std::int32_t { MinVariance = 0, MaxReturn, MaxSharpe, Utility, Importance, Kelly,
                                        TrackingError, MaxDiversification, RiskParity };
  struct Spec
  {
      Objective objective = Objective::MinVariance;
      double lower = -std::numeric_limits<double>::infinity();   // scalar bounds, cash included
      double upper = std::numeric_limits<double>::infinity();
      std::vector<double> asset_lower;  // non-empty replaces the scalars; size checked at solve
      std::vector<double> asset_upper;
      std::optional<double> target_return;
  };
  enum class FailureKind { InvalidArgument, Unsupported };
  struct Failure { FailureKind kind; std::string message; };
  class Data
  {
  public:
      // Copies; covariance symmetrized. Failure when non-finite, asymmetric beyond tolerance or not PD.
      static std::expected<Data, Failure> create(std::size_t n, std::span<const double> mu,
                                                 std::span<const double> covariance);  // empty span = absent
      std::size_t asset_count() const;            // N
      std::size_t weight_count() const;           // N, or N + 1 with cash
      bool has_mu() const; bool has_covariance() const;
      std::span<const double> mu() const;
      std::span<const double> covariance() const; // N x N row-major
      std::optional<double> risk_free() const;
      void set_risk_free(double rate);            // caller checked finiteness
  };
  enum class SolveStatus : std::int32_t { Optimal = 0, Infeasible = 1, Unbounded = 2, MaxIterations = 3 };
  struct Result
  {
      SolveStatus status = SolveStatus::Infeasible;
      std::vector<double> weights;                 // weight_count(); NaN when no point
      double budget_multiplier = NaN;              // NaN when the constraint is absent
      double target_return_multiplier = NaN;
      double volatility_multiplier = NaN;
      double leverage_multiplier = NaN;
      std::vector<double> bound_multipliers;       // weight_count(); lower minus upper
      double residual = 0.0;
      int iterations = 0;
  };
  std::expected<Result, Failure> solve(const Spec& spec, const Data& data);   // problems/solve.hpp
  }
  ```
  Tolerância de simetria: |aᵢⱼ − aⱼᵢ| ≤ 1e-10·max(|aᵢᵢ|, |aⱼⱼ|). Mensagens nomeiam o problema ("covariance is not positive definite", "covariance is not symmetric at (1, 2)").

  **Convenção dos multiplicadores (vale para todos os planos; vai para o header na Task 5):** para o objetivo f em forma de minimização exatamente como na spec (`MIN_VARIANCE`: f = w'Σw, sem ½; objetivos de max: f = −objetivo), ∇f(w) = λ_budget·1 + λ_ret·∇(retorno) + Σᵢ βᵢ eᵢ − β₀·1, com β = bound_multipliers (β ≥ 0 no bound inferior ativo, ≤ 0 no superior) e β₀ a entrada do caixa (o caixa é w₀ = 1 − 1'w, então seu gradiente em w é −1).

  Tradução de `MIN_VARIANCE` para o QP (`min_variance.cpp`, com as linhas de bounds/caixa/retorno em `constraints.cpp`, reutilizadas pelos planos 2–4):
  - G = 2Σ, g = 0.
  - Sem caixa: igualdade 1'w = 1; com `target_return` r, igualdade μ'w = r.
  - Com caixa: sem orçamento; `target_return` vira (μ − rf)'w = r − rf; bounds do caixa lo ≤ 1 − 1'w ≤ hi viram −1'w ≥ lo − 1 e 1'w ≥ 1 − hi (só os finitos).
  - Bounds finitos dos ativos de risco viram linhas eᵢ'w ≥ lo e −eᵢ'w ≥ −hi.
  - O resultado traz w₀ = 1 − 1'w por último; multiplicadores remapeados para a convenção acima.

  `solve` valida, na ordem, e devolve `Failure{InvalidArgument}` com mensagem: objetivo sem os dados que exige (`MIN_VARIANCE` exige covariância; `target_return` exige μ); `asset_lower`/`asset_upper` com tamanho ≠ `weight_count()`; setter fora da matriz da seção 3.3 da spec para o objetivo. Objetivos ainda não implementados devolvem `Failure{Unsupported, "objective not implemented yet"}`.

- [ ] **Step 1: Testes** (Σ e μ de `random_covariance` / `random_vector`; os casos `Optimal` checam KKT do QP construído e a convenção de multiplicadores nas variáveis do usuário por uma função de teste `expect_user_kkt(spec, data, result)` em `min_variance_test.cpp`)
  - `Data.RejectsNonFinite`, `Data.RejectsAsymmetric`, `Data.RejectsIndefinite`, `Data.AcceptsAbsentMuAndCovariance`, `Data.SymmetrizesWithinTolerance`.
  - `MinVariance.UnconstrainedMatchesClosedForm`: n = 6 → Σ⁻¹1 / 1'Σ⁻¹1 em 1e-12 relativo; `budget_multiplier` = 2/(1'Σ⁻¹1) em 1e-12.
  - `MinVariance.FrontierIsParabola`: sem bounds, 5 alvos r; variância = (A r² − 2B r + C)/(AC − B²) com A = 1'Σ⁻¹1, B = 1'Σ⁻¹μ, C = μ'Σ⁻¹μ, em 1e-10 relativo.
  - `MinVariance.LongOnlyBoundsKkt`: `lower = 0`, `upper = 1`, n = 20, sementes 1..20 → `Optimal`, todos os pesos ≥ −1e-14, KKT.
  - `MinVariance.AssetBounds`: vetores por ativo, um deles ativo no superior → peso igual ao bound e β < 0.
  - `MinVariance.InfeasibleBounds`: n = 3, `upper = 0.3` → `Infeasible`, pesos NaN.
  - `MinVariance.CashWithoutTargetIsAllCash`: rf = 0.02 → pesos de risco 0, caixa 1, `budget_multiplier` NaN.
  - `MinVariance.CashFrontierIsCapitalMarketLine`: rf = 0.02, alvo r → vol = |r − rf| / √((μ − rf)'Σ⁻¹(μ − rf)) em 1e-10.
  - `MinVariance.CashBoundsApply`: rf, `lower = 0`, `upper = 1`, alvo acima de rf → caixa ∈ [0, 1], KKT.
  - `MinVariance.AssetBoundsSizeIncludesCash`: com rf, `asset_lower` de tamanho N → `Failure{InvalidArgument}`; N + 1 → ok.
  - `MinVariance.IllConditionedCovariance`: `random_covariance(10, 7, 1e10)` → `Optimal` e KKT na tolerância 1e-8.
  - `MinVariance.TargetReturnEqualToCommonMean`: μ = 0.1·1, alvo 0.1 → `Optimal`, igual ao irrestrito.
  - `MinVariance.FixedAsset`: `asset_lower[0] = asset_upper[0] = 0.25` → w₀ = 0.25 em 1e-14, KKT.
  - `MinVariance.SingleAsset`: n = 1 → w = {1}; com alvo ≠ μ₁ → `Infeasible`.
  - `Solve.MissingCovarianceIsInvalid`, `Solve.TargetReturnWithoutMuIsInvalid`, `Solve.UnimplementedObjectiveIsUnsupported` (cada um dos outros 8 objetivos).
- [ ] **Step 2:** Rodar, ver falhar.
- [ ] **Step 3:** Implementar.
- [ ] **Step 4:** PASS em `gcc-release` e `gcc-debug`.
- [ ] **Step 5: Commit** — `feat: spec, data and MIN_VARIANCE solve`

---

### Task 5: Interface C de spec, dados, solve, resultado e métricas

**Files:**
- Modify: `include/portopt/portopt.h`
- Create: `src/capi/handles.hpp`, `src/capi/spec_api.cpp`, `src/capi/data_api.cpp`, `src/capi/solve_api.cpp`, `src/capi/result_api.cpp`, `src/capi/metrics_api.cpp`
- Create: `src/metrics/portfolio.hpp`, `src/metrics/portfolio.cpp`
- Create: `tests/unit/capi_scope.hpp` (deleters RAII para os três handles, como no VolEngine)
- Test: `tests/unit/capi/spec_test.cpp`, `tests/unit/capi/data_test.cpp`, `tests/unit/capi/solve_test.cpp`, `tests/unit/capi/metrics_test.cpp`, `tests/unit/metrics/portfolio_test.cpp`

**Interfaces:**
- Consumes: `po::Spec`, `po::Data`, `po::Result`, `po::solve` (Task 4).
- Produces (header; texto dos comentários do header é do implementador, conteúdo fixado aqui):
  ```c
  typedef struct po_spec po_spec;
  typedef struct po_data po_data;
  typedef struct po_result po_result;

  #define PO_OBJECTIVE_MIN_VARIANCE 0
  #define PO_OBJECTIVE_MAX_RETURN 1
  #define PO_OBJECTIVE_MAX_SHARPE 2
  #define PO_OBJECTIVE_UTILITY 3
  #define PO_OBJECTIVE_IMPORTANCE 4
  #define PO_OBJECTIVE_KELLY 5
  #define PO_OBJECTIVE_TRACKING_ERROR 6
  #define PO_OBJECTIVE_MAX_DIVERSIFICATION 7
  #define PO_OBJECTIVE_RISK_PARITY 8

  #define PO_SOLVE_OPTIMAL 0
  #define PO_SOLVE_INFEASIBLE 1
  #define PO_SOLVE_UNBOUNDED 2
  #define PO_SOLVE_MAX_ITER 3

  #define PO_CONSTRAINT_BUDGET 0
  #define PO_CONSTRAINT_TARGET_RETURN 1
  #define PO_CONSTRAINT_VOLATILITY 2
  #define PO_CONSTRAINT_LEVERAGE 3

  PO_API po_status po_spec_create(int32_t objective, po_spec** out);
  PO_API po_status po_spec_set_bounds(po_spec* s, double lower, double upper);
  PO_API po_status po_spec_set_asset_bounds(po_spec* s, const double* lower, const double* upper, size_t n);
  PO_API po_status po_spec_set_target_return(po_spec* s, double r);
  PO_API void po_spec_free(po_spec* s);

  PO_API po_status po_data_create(size_t n, const double* mu, const double* cov, po_data** out);
  PO_API po_status po_data_set_risk_free(po_data* d, double rate);
  PO_API po_status po_data_asset_count(const po_data* d, size_t* out);
  PO_API po_status po_data_weight_count(const po_data* d, size_t* out);
  PO_API void po_data_free(po_data* d);

  PO_API po_status po_solve(const po_spec* s, const po_data* d, po_result** out);

  PO_API po_status po_result_status(const po_result* r, int32_t* out);
  PO_API po_status po_result_weight_count(const po_result* r, size_t* out);
  PO_API po_status po_result_weights(const po_result* r, double* out, size_t n);
  PO_API po_status po_result_iterations(const po_result* r, int32_t* out);
  PO_API po_status po_result_multiplier(const po_result* r, int32_t constraint, double* out);
  PO_API po_status po_result_bound_multipliers(const po_result* r, double* out, size_t n);
  PO_API po_status po_result_residual(const po_result* r, double* out);
  PO_API void po_result_free(po_result* r);

  PO_API po_status po_portfolio_return(const po_data* d, const double* w, size_t n, double* out);
  PO_API po_status po_portfolio_volatility(const po_data* d, const double* w, size_t n, double* out);
  ```
  O header documenta: a convenção de multiplicadores da Task 4; que NaN significa restrição ausente; que `*out` é `NULL` depois de falha de `*_create`/`po_solve`; que `n` dos acessores de vetor deve ser exatamente a contagem (senão `PO_INVALID_ARGUMENT`); que `po_solve` é seguro em paralelo com handles somente leitura; que objetivos fora de `MIN_VARIANCE` devolvem `PO_UNSUPPORTED` nesta versão.
- Produces (C++): `po::metrics::portfolio_return(const Data&, std::span<const double> w) -> double`, `portfolio_volatility(...) -> double`; com caixa, retorno = Σwᵢμᵢ + w₀·rf e vol só dos ativos de risco.
- Regras: ponteiro de handle nulo → `PO_INVALID_HANDLE`; ponteiro de saída nulo → `PO_INVALID_ARGUMENT`; setter com valor fora do domínio (NaN; `lower > upper`; `lower = +∞`; `upper = −∞`) → `PO_INVALID_ARGUMENT` e spec inalterada; `po_spec_set_bounds` limpa `asset_*` e `po_spec_set_asset_bounds` substitui os escalares (a última chamada vale); `n = 0` em `po_data_create` → `PO_INVALID_ARGUMENT`; `Failure` → `PO_INVALID_ARGUMENT` ou `PO_UNSUPPORTED` com a mensagem.

- [ ] **Step 1: Testes**
  - `CapiSpec.RejectsUnknownObjective` (−1 e 9), `CapiSpec.InvalidBoundsLeaveSpecUnchanged`, `CapiSpec.LastBoundsCallWins`.
  - `CapiData.CopiesInputs` (alterar o array depois de criar não muda o solve), `CapiData.NullCovarianceAllowedAtCreate`, `CapiData.RejectsNaN`, `CapiData.WeightCountWithCash` (N → N + 1).
  - `CapiSolve.MinVarianceEndToEnd`: n = 4, mesma resposta que `po::solve` bit a bit; status `PO_SOLVE_OPTIMAL`; `PO_CONSTRAINT_TARGET_RETURN` NaN sem alvo.
  - `CapiSolve.InfeasibleIsResultNotError`: devolve `PO_OK`, status `PO_SOLVE_INFEASIBLE`, pesos NaN.
  - `CapiSolve.UnsupportedObjective`: `PO_OBJECTIVE_MAX_SHARPE` → `PO_UNSUPPORTED`, `*out == NULL`, `po_last_error()` não vazio.
  - `CapiSolve.ResultOutlivesSpecAndData`: liberar spec e dados e ler o resultado.
  - `CapiSolve.WrongVectorSizeIsInvalid`: `po_result_weights` com n errado.
  - `CapiSolve.NullHandles`: cada função com handle nulo → `PO_INVALID_HANDLE`; cada `*_free(NULL)` não faz nada.
  - `Metrics.ReturnAndVolatility`: pesos iguais, comparar com w'μ e √(w'Σw); com caixa, retorno inclui w₀·rf.
  - `CapiMetrics.WrongWeightCountIsInvalid`, `CapiMetrics.MissingMuIsInvalid`.
- [ ] **Step 2:** Rodar, ver falhar.
- [ ] **Step 3:** Implementar.
- [ ] **Step 4:** `ctest --preset gcc-release` PASS (inclui `c_header_check` e `style`).
- [ ] **Step 5: Commit** — `feat(capi): spec, data, solve, result and portfolio metrics`

---

### Task 6: Pacote Python

**Files:**
- Create: `python/pyproject.toml` (do VolEngine: nome `portopt`, versão `0.1.0`, descrição "Portfolio optimization over the PortOpt C API", grupo `dev` = pytest, mypy, hatchling; sem grupo `bench`), `python/hatch_build.py` (do VolEngine, nomes trocados)
- Create: `python/src/portopt/__init__.py`, `_native.py`, `_spec.py`, `_data.py`, `_result.py`, `py.typed`
- Test: `python/tests/conftest.py`, `python/tests/test_loader.py`, `python/tests/test_min_variance.py`

**Interfaces:**
- Consumes: o header da Task 5.
- Produces:
  ```python
  class Objective(enum.IntEnum): MIN_VARIANCE = 0; MAX_RETURN = 1; ...; RISK_PARITY = 8
  class SolveStatus(enum.IntEnum): OPTIMAL = 0; INFEASIBLE = 1; UNBOUNDED = 2; MAX_ITER = 3

  class Spec:
      def __init__(self, objective: Objective, *, bounds: tuple[float | ArrayLike, float | ArrayLike] | None = None,
                   target_return: float | None = None) -> None
  class Data:
      def __init__(self, mu: ArrayLike | None, cov: ArrayLike | None, *, risk_free: float | None = None) -> None
      asset_count: int; weight_count: int
      def portfolio_return(self, weights: ArrayLike) -> float
      def volatility(self, weights: ArrayLike) -> float
  @dataclass(frozen=True)
  class Result:
      status: SolveStatus; weights: NDArray[np.float64]; iterations: int
      budget_multiplier: float; target_return_multiplier: float; bound_multipliers: NDArray[np.float64]
  def solve(spec: Spec, data: Data) -> Result
  class PortOptError(Exception)            # PO_INTERNAL_ERROR, PO_OUT_OF_MEMORY, PO_INVALID_HANDLE
  ```
  - `_native.py`: descoberta por `PORTOPT_LIB`, senão `src/portopt/_lib` (`libportopt.dll`, `portopt.dll`, `libportopt.so`, `libportopt.dylib`); checagem `po_abi_version() == 0`; tabela `PROTOTYPES`; posse dos handles (cada classe libera o seu no `__del__`/finalizer).
  - `INVALID_ARGUMENT` → `ValueError(po_last_error())`; `UNSUPPORTED` → `NotImplementedError(...)`.
  - `bounds` com escalares chama `po_spec_set_bounds`; com arrays, `po_spec_set_asset_bounds` (escalar e array misturados: o escalar é expandido).
  - `Result` copia tudo do `po_result` e libera o handle na hora (resultado é pequeno).
  - Arrays convertidos com `np.ascontiguousarray(x, dtype=np.float64)`; `cov` 2-D quadrada, `mu` 1-D, senão `ValueError` antes da chamada C.

- [ ] **Step 1: Testes** (`conftest.py` pula tudo com mensagem clara se a biblioteca não for encontrada)
  - `test_loader.py::test_abi_version_matches`, `::test_missing_library_message_names_env_var`.
  - `test_min_variance.py::test_unconstrained_matches_closed_form` (numpy `solve`, 1e-12), `::test_long_only_weights_nonnegative_and_sum_to_one`, `::test_frontier_with_list_comprehension` (5 alvos, variância crescente nos dois lados do mínimo), `::test_cash_weight_is_last`, `::test_infeasible_returns_status_not_exception`, `::test_invalid_covariance_raises_value_error`, `::test_unsupported_objective_raises_not_implemented`, `::test_data_accepts_fortran_and_int_arrays` (`np.asfortranarray` e `dtype=int` dão o mesmo resultado que float C-order), `::test_metrics_match_numpy`.
- [ ] **Step 2:** `uv sync --project python` e rodar: `PORTOPT_LIB=$PWD/build/gcc-release/libportopt.dll uv run --project python pytest python/tests -q` → FAIL.
- [ ] **Step 3:** Implementar.
- [ ] **Step 4:** Mesmo comando → PASS; `uv run --project python mypy --strict python/src/portopt` sem erros; `python tests/style/check_style.py` PASS.
- [ ] **Step 5: Commit** — `feat(python): portopt package over the C API`

---

### Task 7: Documentos do projeto e CI

**Files:**
- Create: `README.md`, `ARCHITECTURE.md`, `CONTRIBUTING.md`, `CHANGELOG.md`, `CITATION.cff`, `docs/guide/math.md`, `docs/guide/validation.md`
- Create: `.github/workflows/ci.yml`

**Interfaces:**
- Consumes: tudo acima.

- [ ] **Step 1:** `ARCHITECTURE.md` na forma do VolEngine, descrevendo só o que existe: camadas `core` ← `linalg`, `qp` ← `problems` ← `metrics`, `capi` (`linalg` não depende de nada; `qp` depende só de `linalg`); tabela de responsabilidades por arquivo; regras de estabilidade da spec §5.1 (tabela do VolEngine adaptada); anatomia de um `po_solve` (validação → tradução para QP → Goldfarb-Idnani → remapeamento para as variáveis do usuário → publicação do handle).
- [ ] **Step 2:** `docs/guide/math.md`: notação da spec §3.1, `MIN_VARIANCE` com e sem caixa, as linhas de restrição de `constraints.cpp` e a convenção de multiplicadores. Os planos 2–4 acrescentam suas seções.
- [ ] **Step 3:** `docs/guide/validation.md`: a política (KKT, formas fechadas, consistência, invariantes; sem oráculo externo, porque as equações são do autor) e as tolerâncias usadas: KKT 1e-10 relativo (1e-8 com condição 1e10), formas fechadas 1e-12, fronteira 1e-10, bounds 1e-14.
- [ ] **Step 4:** `README.md` (o que é, exemplo Python de `MIN_VARIANCE` com fronteira, build), `CONTRIBUTING.md` (o do VolEngine sem benchmarks, oráculos e rugarch; com "do not copy or translate quadprog code (GPL)"; variável `PORTOPT_LIB`), `CHANGELOG.md` com `## [Unreleased]` listando o plano 1, `CITATION.cff` (título "PortOpt", abstract "Portfolio optimization in C++23 behind a stable C API, with a Python package over it.", versão 0.1.0, MIT, keywords portfolio-optimization, quadratic-programming, finance).
- [ ] **Step 5:** `ci.yml` a partir do VolEngine: matriz windows-msvc, linux-gcc (GCC 15 do PPA), macos-clang com `native-release`; build, ctest, pytest com `PORTOPT_LIB`, mypy, clang-tidy sobre `src/` no Linux.
- [ ] **Step 6: Verificar** — `ctest --preset gcc-release` PASS, pytest PASS, e `clang-tidy -p build/gcc-release` sobre cada `.cpp` de `src/` sem avisos (se o clang-tidy estiver instalado localmente; senão o CI cobre).
- [ ] **Step 7: Commit** — `docs: architecture, math, validation, CI`
