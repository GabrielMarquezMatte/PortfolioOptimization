# PortOpt — Plano 8: bits portáveis e o cruzamento do MAX_RETURN

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** CI verde nas três plataformas com os mesmos bits de solução e de instância de teste em todas, e o `MAX_RETURN` sem caso pinado, antes da v0.1.0 final.

**Architecture:** Primeiro o cruzamento do `MAX_RETURN` passa a sair de uma tangente exata da face (null-space, referencial do shift) com certificado KKT num ponto, validado contra as duas versões (Windows e Linux) da Σ do caso hoje pinado. Depois `src/` perde a libm transcendental (`hypot`, `log` → `po::portable`), sob A/B de performance. Por fim os geradores de teste passam a usar só `po::portable` e `engine()` cru, e a documentação é remedida.

**Tech Stack:** C++23, CMake 3.28, GoogleTest, MinGW-GCC 16 (`gcc-release`, `gcc-bench`), MSVC (`native-release`), Docker (`quay.io/pypa/manylinux_2_28_x86_64` com `gcc-toolset-15`), Python 3.14 (`tests/benchmarks/quick.py`), GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-11-18-portopt-portable-and-max-return-design.md`; pesquisa `docs/research/2026-11-18-max-return-path-and-portable-generators.md` ("pesquisa §N").

## Global Constraints

- Valem as dos planos 1–7: `STYLEGUIDE.md`, commits **sem** `Co-Authored-By`, `-ffp-contract=off` em todo alvo de produção e de teste, sem `long double` nem `-ffast-math`, `docs/superpowers/` em português, documentos do repositório em inglês.
- `src/` só usa `+ − × ÷`, `sqrt` e `fma` da biblioteca padrão de matemática ao fim do plano (`grep -rnE "std::(hypot|log|exp|pow|sin|cos|log10|log2|cbrt|tan|atan)" src` vazio).
- `tests/` não usa `std::uniform_*_distribution`, `std::normal_distribution` nem `std::log/exp/cos/sin/pow` para gerar dados (fora `std::log10` em nomes de propriedade).
- Proibido para deixar o CI verde: `ctest -E`, `GTEST_SKIP`, `#ifdef` de SO em teste, banda de `expect_user_kkt`/`expect_risk_parity_kkt` afrouxada, novo pin por nome.
- Uma cópia só da matemática portável: `src/core/portable_math.{hpp,cpp}`, namespace `po::portable`.
- O porte do fdlibm mantém no topo de `portable_math.cpp` o aviso da Sun: "Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved. Developed at SunSoft, a Sun Microsystems, Inc. business. Permission to use, copy, modify, and distribute this software is freely granted, provided that this notice is preserved."
- Testes locais: `ctest --preset gcc-release` e `ctest --preset native-release` (MSVC). Não rodar `gcc-debug`.
- Performance (condição do autor): A/B de `docs/performance/measuring.md` (`quick.py run` duas vezes antes; A e B na mesma sessão). Caso pior que `noise` → **parar e levar ao autor** com a tabela.
- Toda mudança de bits de solução: `portopt_profile --max-n 100 --json` antes e depois, `quick.py compare-solutions A B`, resumo no corpo do commit.
- Tags (`v0.1.0-rc.2`, `v0.1.0`) e qualquer push de tag são do autor.

## Review Focus

1. **Argumentos extremos de `po::portable::hypot`** (subnormais, 1e-300 com 1e300, zeros, um zero): resultado igual a `sqrt(a²+b²)` exato arredondado, sem underflow/overflow intermediário — os Givens do QP recebem isso em problemas mal escalados. (Task 3: `PortableMath.HypotScalesExtremeArguments`.)
2. **`po::portable::log` em subnormais e em 1** (o Spinu leva `y_i` a valores minúsculos): `log(denorm_min)` dentro de 1 ulp de −744.44007192138126…, `log(1) == 0` exato. (Task 3: `PortableMath.LogHandlesSubnormalsAndOne`.)
3. **Face sem família de leverage** (long-only comum, só budget e bounds): a tangente vale também ali; um MAX_RETURN de n = 3 long-only com cap no meio do segmento sai igual ao de antes dentro da banda. (Task 2: `ParametricPath.TangentCrossingMatchesTheSolveWithoutLeverage`.)
4. **Tangente que sai da face antes do cap** (o segmento termina num breakpoint antes de `vol = s`): o certificado recusa e a bissecção segue, nunca uma resposta que viola uma linha. (Task 2: `ParametricPath.TangentCrossingRejectsAPointOffItsFace`.)
5. **MSVC × MinGW na mesma máquina**: os testes de bits (Tasks 3, 4, 5) passam nos dois presets antes do push, para que uma divergência apareça localmente e não só no CI.

---

### Task 1: Fixtures das duas Σ do caso pinado — teste que falha

**Files:**
- Create: `tests/unit/fuzz/known_shape_fixture.hpp`
- Modify: `tests/unit/fuzz/max_return_fuzz_test.cpp` (dois testes novos, no fim)

**Interfaces:**
- Produces:
  ```cpp
  namespace po::testing {
  struct ShapeFixture {
      std::size_t n;
      std::vector<double> mu;          // n
      std::vector<double> covariance;  // n x n, row-major
      double risk_free;
      double leverage;
      std::vector<double> fractions;   // cap fractions of the fuzz
      std::vector<double> caps;        // the caps the fuzz computed on that platform, same order
  };
  ShapeFixture known_shape_windows();  // shape seed 3963545680699380072, as MinGW-GCC on Windows generates it
  ShapeFixture known_shape_linux();    // the same shape, as GCC 15 + glibc generates it
  GeneratedProblem fixture_problem(const ShapeFixture& fixture);  // Box None, cash, spec.leverage = fixture.leverage
  }
  ```
  Todos os doubles em hex-float (`%a`), exatos.

- [ ] **Step 1:** Escrever no scratchpad (não no repositório) `dump_shape.cpp`, que inclui `tests/unit/fuzz/random_problems.hpp`, acha em `with_leverage(fuzz_shapes(kLeverageShapes, 3))` o shape de seed `3963545680699380072`, gera o problema e imprime em `%a` n, μ, Σ, `risk_free`, `leverage` e os caps das frações de `kCapFractions` calculados como `LeverageFuzz` os calcula (copiar a função de caps de `max_return_fuzz_test.cpp`). Compilar com `-std=c++23 -O2 -ffp-contract=off` contra `portopt_core`.
- [ ] **Step 2:** Rodar no Windows (MinGW, `gcc-release`) e no Linux: `docker run --rm -v "$PWD:/work" -w /work quay.io/pypa/manylinux_2_28_x86_64 sh -c ". /opt/rh/gcc-toolset-15/enable && <cmake configure+build de portopt_core com PO_BUILD_TESTS=OFF> && g++ ... dump_shape.cpp ... && ./dump_shape"`. Conferir no Linux: 561 das 576 entradas de Σ diferem do Windows e μ é igual (`validation.md` "Other platforms"). Se não bater, **parar** e mostrar ao autor.
- [ ] **Step 3:** Gravar as duas saídas em `known_shape_fixture.hpp` (comentário no topo: de onde vieram, o comando do Step 2, e que são congeladas — não mudam quando os geradores mudarem).
- [ ] **Step 4: Testes** em `max_return_fuzz_test.cpp`, reusando `checked_against_reference(problem, cap, condition_tolerance(1e8))` sem nenhum relaxamento:
  ```cpp
  TEST(MaxReturnFuzz, KnownShapeUnderTheWindowsCovariance)  // every fraction: OPTIMAL, full band
  TEST(MaxReturnFuzz, KnownShapeUnderTheLinuxCovariance)    // fractions 0.3, 0.5, 0.7, 0.9: OPTIMAL, full band
  ```
  Cada um: para cada `(fraction, cap)` da fixture, `SCOPED_TRACE(fraction)`, `std::optional<Result> r = checked_against_reference(fixture_problem(f), cap, tol)`; `ASSERT_TRUE(r)`; `EXPECT_EQ(r->status, SolveStatus::Optimal)`.
- [ ] **Step 5:** `ctest --preset gcc-release -R KnownShape` → os dois falham hoje: Windows na fração 0,7 (complementaridade 1.006e-7 > 1e-7), Linux `MAX_ITER` em 0,3–0,9. Anotar as falhas exatas para o commit.
- [ ] **Step 6: Commit** — `test: the known MAX_RETURN shape under its Windows and Linux covariances` (corpo: as falhas do Step 5).

---

### Task 2: Cruzamento por tangente exata com certificado

**Files:**
- Modify: `src/problems/parametric_path.hpp`, `src/problems/parametric_path.cpp`, `src/problems/max_return.cpp`
- Modify: `tests/unit/fuzz/max_return_fuzz_test.cpp` (remover `KnownResidual`, `kKnownResidual`, `pinned_against_reference`; `leverage_case` chama só `checked_against_reference`)
- Create: `tests/unit/problems/parametric_path_test.cpp` (e registrá-lo em `CMakeLists.txt`, lista de `portopt_tests`)

**Interfaces:**
- Consumes: Task 1 (os dois testes `KnownShape*`).
- Produces (em `parametric_path.hpp`):
  ```cpp
  // A solution of Q(t) at the crossing vol(w(t)) = cap, multipliers in the convention of PathPoint::solution.
  struct Crossing { double t; qp::Solution solution; };
  // From the face and shift of point: the exact tangent dw on that face (null space of the active normals over
  // the free coordinates, QR of L'Z), anchored at a face solve at point.t; the root of ||L'(w_a + h dw)|| = cap by
  // the cancellation-free formula; certified by the KKT conditions of Q(t*) at one point (research §2). nullopt
  // when the certificate fails or the reduced Hessian has a pivot below its noise.
  std::optional<Crossing> tangent_crossing(const Path& path, const PathPoint& point);
  ```
  `same_segment(earlier, later)` muda de semântica: mesma face (linhas explícitas ativas iguais; família ativa nos dois ou em nenhum; mesmo suporte com sinais e mesmo conjunto de zeros D), sem comparar `gross_exposure_signs`.

- [ ] **Step 1: Testes que falham** em `parametric_path_test.cpp`, sobre `make_path`/`solve_at` de `po::problems`:
  - `ParametricPath.SameSegmentIgnoresTheSignVectorRepresentation`: dois `PathPoint` com os mesmos `active_inequalities` e `x` (suporte {0: +, 2: −}, zero em 1), mas `gross_exposure_signs` `{1,1,−1}` num e `{1,−1,−1}` no outro → `same_segment == true`; trocar o sinal de `x[2]` → `false`.
  - `ParametricPath.TangentCrossingMatchesTheSolveWithoutLeverage` (Review Focus 3): n = 3 long-only, μ e Σ literais (escolha do builder, Σ bem condicionada), cap no meio do primeiro segmento; `point = solve_at(path, t)` com t dentro dele; `tangent_crossing(path, point)` tem valor; os pesos a 1e-12 dos do `po::solve` de antes desta task (gravar os de antes como literais).
  - `ParametricPath.TangentCrossingRejectsAPointOffItsFace` (Review Focus 4): mesmo problema, cap acima do fim do segmento de `point` (o próximo breakpoint vem antes de `vol = s`) → `nullopt`.
- [ ] **Step 2:** `ctest --preset gcc-release -R "ParametricPath|KnownShape"` → falham (símbolo ausente / semântica antiga).
- [ ] **Step 3:** Implementar em `parametric_path.cpp` (pesquisa §5a itens 1–4): identidade de face em `same_segment`; o solve de face (coordenadas livres = complemento de D; N das linhas ativas, budget/cash exatos, linhas da família como −σ_S; `linalg::householder_qr` de N; Z; QR de LᵀZ com `path.factor`; `dw = Z H_r⁻¹ Zᵀ v`; `w_a` em `point.t`; multiplicadores por `linalg::least_squares` em `N u = Σw − t v`, no referencial de `point.shift`); a raiz pela fórmula de `crossing_fraction` com `start = Lᵀw_a`, `direction = Lᵀdw`; o certificado: primal com `row_tolerance`/`family_holds`/sinais do suporte, dual explícito `u_i ≥ −tol_d`, dual da família na forma invariante (`−ζ_v + 2t*a ≤ s_f,i ρ_i ≤ ζ_v` em D, pesquisa §2 item 3), `tol_d = kNoiseFactor·n·ε·(t*|v|max + |Σ|max|w|max)`. Devolver os multiplicadores no referencial de Q(t) (`to_unshifted`).
- [ ] **Step 4:** Em `max_return.cpp`: em `bisect`, quando `same_segment(bracket.lower, bracket.upper)`, `tangent_crossing(path, bracket.lower)`; com valor → `answer(data, path, c.solution, c.t, -path.cap / c.t)`; sem valor → segue a bissecção pelo ponto médio. Em doubles adjacentes: `tangent_crossing` de `lower` (de `upper` se `lower.t == 0`), depois do outro; nenhum → `broke_down`. Remover `kCrossingSanity`, `kAdjacentEnd`, `confirms`, `at_adjacent`, `at_crossing` e, se ficarem sem uso, `crossing_fraction`/`on_segment` do header. Remover o pin do fuzz.
- [ ] **Step 5:** `ctest --preset gcc-release` e `ctest --preset native-release` inteiros verdes, incluindo `KnownShape*` e `MaxReturnFuzz.LeverageFuzz` sem pin.
- [ ] **Step 6:** `compare-solutions` (antes = `HEAD` da Task 1): só casos `max_return/*` mudam, cada um dentro da banda, nenhum status muda. A/B `quick.py` filtrando `max_return` (regra das Global Constraints).
- [ ] **Step 7:** `validation.md`: a seção "The known residual" vira história curta (o que era, e que a tangente o resolve, com `vol/s − 1` e o erro do retorno medidos agora no caso); o comentário de `LeverageFuzz` perde "one is the known residual".
- [ ] **Step 8: Commit** — `fix(max_return): exact tangent crossing with a one-point KKT certificate` (corpo: compare-solutions e A/B).

---

### Task 3: `po::portable` e a libm fora de `src/`

**Files:**
- Create: `src/core/portable_math.hpp`, `src/core/portable_math.cpp`, `tests/unit/core/portable_math_test.cpp`
- Modify: `CMakeLists.txt` (fonte em `portopt_core`, teste em `portopt_tests`), `src/qp/active_set.hpp:132`, `src/riskparity/spinu.cpp:72`

**Interfaces:**
- Produces:
  ```cpp
  namespace po::portable {
  double exp(double x);              // fdlibm 5.3 e_exp.c; finite x; < 1 ulp
  double log(double x);              // fdlibm 5.3 e_log.c; finite x > 0, subnormals included; < 1 ulp
  double hypot(double a, double b);  // power-of-two scaling, sqrt, unscaling; no intermediate over/underflow
  }
  ```
  Palavras por `std::bit_cast<std::uint64_t>`; constantes em hex literal; `assert` fora do domínio; `#pragma fp_contract(off)` sob `_MSC_VER` no `.cpp`.

- [ ] **Step 1: Testes que falham** em `portable_math_test.cpp`:
  - `PortableMath.ExpHasTheCommittedBits`, `LogHasTheCommittedBits`, `HypotHasTheCommittedBits`: ~10 entradas cada, comparadas por `std::bit_cast<std::uint64_t>` com valores gravados. Os valores gravados são a saída da implementação no MinGW, e cada um tem de estar a ≤ 1 ulp do valor de referência em 60 dígitos (script `mpmath` no scratchpad; gravar no comentário do teste o comando e a maior distância em ulp encontrada).
  - `PortableMath.HypotScalesExtremeArguments` (Review Focus 1): `hypot(0x3p-1074, 0x4p-1074) == 0x5p-1074` (subnormais, exato), `hypot(0x1p-1074, 0x1p-1074) == 0x1p-1074` (√2 denorm_min arredonda para denorm_min), `hypot(1e-300, 1e300) == 1e300`, `hypot(0, 0) == 0`, `hypot(0, -3) == 3`, `hypot(3, 4) == 5`, `hypot(0x1p1000, 0x1p1000) == 0x1.6a09e667f3bcdp1000`.
  - `PortableMath.LogHandlesSubnormalsAndOne` (Review Focus 2): `log(1.0) == 0.0`; `log(std::numeric_limits<double>::denorm_min())` a ≤ 1 ulp de −744.44007192138126231410729844608163411300…; `log(0x1p-1022)` idem de −708.39641853226410621…
  - `PortableMath.FmaIsCorrectlyRounded`: `std::fma(1 + 0x1p-52, 1 - 0x1p-52, -1.0) == -0x1p-104`, `std::fma(0x1.0000001p0, 0x1.0000001p0, -0x1.0000002p0) == 0x1p-56`, mais dois casos de arredondamento difícil (escolha do builder). Se falhar numa plataforma: spec §6.4 (`fma` próprio por Dekker/Veltkamp em `po::portable`), não tolerância.
- [ ] **Step 2:** `ctest --preset gcc-release -R PortableMath` → falha (ausente).
- [ ] **Step 3:** Implementar `portable_math.{hpp,cpp}`. Rodar os testes em `gcc-release` **e** `native-release`: mesmos bits (Review Focus 5).
- [ ] **Step 4:** Antes da troca: `cmake --build --preset gcc-bench`, `quick.py run` duas vezes (descartar a primeira), depois `quick_a.json`; `portopt_profile --max-n 100 --json solutions_a.json`. Docker Linux (`manylinux_2_28_x86_64`, `gcc-toolset-15`, `-DCMAKE_BUILD_TYPE=Release -DPO_BUILD_BENCHMARKS=ON -DPO_BUILD_TESTS=OFF`): o mesmo par `quick.py run --exe build/linux-bench/portopt_benchmarks` (duas vezes, a segunda é a A).
- [ ] **Step 5:** Trocar `std::hypot` → `po::portable::hypot` em `make_givens` e `std::log` → `po::portable::log` em `spinu.cpp:72`. `grep -rnE "std::(hypot|log|exp|pow|sin|cos)" src` vazio.
- [ ] **Step 6:** B nos dois ambientes; `quick.py compare` A × B. Pior que `noise` em qualquer caso → **parar e levar ao autor** com as duas tabelas (opções do spec §6.3).
- [ ] **Step 7:** `compare-solutions solutions_a.json solutions_b.json`: nenhum status muda; pesos e objetivos dentro da banda; casos com iterações diferentes listados no commit. Status mudou → parar e levar ao autor.
- [ ] **Step 8:** `ctest --preset gcc-release` e `ctest --preset native-release` verdes. Registrar as tabelas A/B em `docs/performance/` (arquivo novo `portable-math.md`, inglês, curto: máquina, comandos, medianas, veredito; e que o macOS não foi medido, spec §6.3).
- [ ] **Step 9: Commit** — `perf: portable exp, log and hypot; no libm transcendental in src/` (corpo: compare-solutions e veredito do A/B).

---

### Task 4: Soluções-ouro (guarda de bits da biblioteca)

**Files:**
- Create: `tests/unit/solve/determinism_test.cpp` (registrar em `CMakeLists.txt`)

**Interfaces:**
- Consumes: a biblioteca final das Tasks 2–3.

- [ ] **Step 1:** Um `TEST(Determinism, <Objective>HasTheCommittedBits)` por objetivo (os nove de `po::Objective`), n = 4 (5 com cash onde o objetivo aceita), μ, Σ (e retornos/benchmark para `TRACKING_ERROR`) literais **não gerados**, cada um com uma restrição que o exercite (bounds, leverage, target return, target volatility ou risk budget, escolha do builder, ao menos um com leverage ativa e um MAX_RETURN com cap ativo). Compara status, `iterations`, pesos, os quatro multiplicadores e `residual` por `std::bit_cast<std::uint64_t>` com hex-floats gravados. Comentário no topo: os valores são regravados no mesmo commit de qualquer mudança legítima de bits, com o `compare-solutions` no corpo.
- [ ] **Step 2:** Gravar os valores do MinGW; `ctest --preset native-release -R Determinism` verde com os mesmos valores (se não: spec §6.4, causa na raiz).
- [ ] **Step 3: Commit** — `test: committed solution bits for every objective`

---

### Task 5: Geradores de teste portáveis

**Files:**
- Modify: `tests/unit/random_covariance.hpp`, `tests/unit/fuzz/random_problems.hpp`, `tests/unit/fuzz/max_return_fuzz_test.cpp:120`, `tests/unit/linalg/qr_test.cpp:67`, `tests/unit/problems/tracking_error_test.cpp:354`, `tests/unit/fuzz/risk_parity/risk_parity_fuzz_test.cpp:202`, `tests/unit/fuzz/risk_parity/risk_parity_reference.hpp:369-370,552`, `tests/unit/riskparity/spinu_test.cpp:58`, `tests/unit/problems/max_sharpe_test.cpp:221,254,286`, `tests/unit/riskparity/objective_test.cpp:303`
- Create: `tests/unit/core/generators_test.cpp` (registrar)

**Interfaces:**
- Consumes: `po::portable::exp`, `log` (Task 3).
- Produces (em `random_covariance.hpp`, namespace `po::testing`):
  ```cpp
  double exp10(double x);                                     // po::portable::exp(x * 0x1.26bb1bbb55516p+1)
  double standard_normal(std::mt19937_64& engine);            // Marsaglia polar: only log and sqrt
  double uniform_between(std::mt19937_64& engine, double low, double high);  // low + (high - low) * uniform_open
  std::uint64_t uniform_index(std::mt19937_64& engine, std::uint64_t count); // rejection on engine(), [0, count)
  std::int64_t uniform_integer(std::mt19937_64& engine, std::int64_t low, std::int64_t high);  // [low, high]
  ```

- [ ] **Step 1: Teste que falha** `Generators.HaveTheCommittedBits` em `generators_test.cpp`: os 8 primeiros `standard_normal` da semente 42, `exp10(-7.5)`, `uniform_index(engine, 17)` × 5, `uniform_integer(engine, -8, 8)` × 5, e `random_covariance(4, 7, 1e8)` inteira (16 entradas), por bits, com os valores gravados do MinGW. Ver falhar (funções ausentes).
- [ ] **Step 2:** Implementar as funções; trocar `std::exp/log` dos autovalores em `random_covariance` por `po::portable`; `std::pow(10, x)` → `exp10(x)` nos arquivos listados; as distribuições de `<random>` → `uniform_between`/`uniform_integer`. A elipse de `risk_parity_reference.hpp`: parametrização racional (`((1−τ²)/(1+τ²), 2τ/(1+τ²))`, τ em pontos que cubram o mesmo arco) ou sin/cos portáteis, escolha do builder, com a resolução do scan igual ou melhor que hoje.
- [ ] **Step 3:** `Generators.HaveTheCommittedBits` verde em `gcc-release` e `native-release` (mesmos bits).
- [ ] **Step 4:** `ctest --preset gcc-release` e `native-release` inteiros. **Triagem, um teste por vez, uma linha por teste no commit** (spec §6.5): (a) contagem ou valor medido gravado → remedir e atualizar; (b) referência de resolução limitada → melhorar a referência; (c) defeito do solver → corrigir se local e pequeno (com teste próprio e `compare-solutions`), senão **parar e levar ao autor** com o caso (shape, seed, o que falha, por quanto). Nada de banda afrouxada ou pin.
- [ ] **Step 5:** `grep -rnE "std::(uniform_(int|real)_distribution|normal_distribution)|std::(log|exp|cos|sin|pow)\(" tests/unit` só mostra `std::log10` em nomes de propriedade.
- [ ] **Step 6: Commit** — `test: portable generators; every test instance has the same bits everywhere` (corpo: a triagem).

---

### Task 6: `validation.md`, CHANGELOG e docs remedidos

**Files:**
- Modify: `docs/guide/validation.md`, `CHANGELOG.md`, `ARCHITECTURE.md` (se citar libm/`hypot`), comentários de teste com contagens (`LeverageFuzz`, `AgainstFrontierReference` etc.)

- [ ] **Step 1:** Remedir em `gcc-release` cada número de `validation.md` que vem de instância gerada (contagens por status, piores concordâncias e resíduos por cond, tempos de execução) e atualizar texto e comentários dos testes. A seção "Other platforms" passa a dizer: os geradores e a biblioteca dão os mesmos bits em Windows (MinGW, MSVC), Linux e macOS, conferido por `Generators.HaveTheCommittedBits` e `Determinism.*`.
- [ ] **Step 2:** CHANGELOG, seção `[0.1.0]` (não lançada; a data fica para o dia da tag): `Changed` — bits de solução do QP (hypot portátil) e do `MAX_RETURN` no cruzamento, com o resumo do `compare-solutions` das Tasks 2 e 3; `Fixed` — o caso do `MAX_RETURN` que virava `MAX_ITER` sob perturbação de 3,6e-13 de Σ; remover de "Known limitations" o item das plataformas e a frase "Portable generators … are the next plan".
- [ ] **Step 3:** `python python/scripts/release_version.py check` sai 0; `ctest --preset gcc-release` verde (o teste `style` inclusive).
- [ ] **Step 4: Commit** — `docs: validation remeasured with the portable generators`

---

### Task 7: CI verde nas três plataformas e entrega

- [ ] **Step 1:** `git push origin master`; `gh run watch` de `ci.yml`, `cpp.yml` e `wheels.yml` até o fim.
- [ ] **Step 2:** Os três verdes, os três jobs `Build & Test` com todos os passos `success` (nenhum `skipped` além dos de outra plataforma). Vermelho → diagnosticar com `gh run view --log-failed`; divergência de bits → spec §6.4 (causa na raiz, commit novo); causa não achada → **parar e levar ao autor**.
- [ ] **Step 3: PARAR e entregar ao autor:** (1) `git tag v0.1.0-rc.2 && git push origin v0.1.0-rc.2`; (2) depois do `release.yml` verde, num venv limpo `uv pip install --pre pyportopt==0.1.0rc2` e o exemplo do README; (3) a data de `## [0.1.0]` no CHANGELOG = dia da tag, commit, push, `git tag v0.1.0 && git push origin v0.1.0`.
