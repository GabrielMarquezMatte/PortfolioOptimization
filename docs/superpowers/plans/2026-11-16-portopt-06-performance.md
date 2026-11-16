# PortOpt — Plano 6: performance, no método do VolEngine

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Medir onde o tempo da biblioteca vai (suíte Google Benchmark e linha de base) e aplicar as otimizações que **não alteram nenhum bit** do resultado (risco 0), na ordem do ganho esperado, uma por commit com A/B, antes de fechar a v1.

**Architecture:** O método do VolEngine porta quase inteiro: Google Benchmark via FetchContent, suíte rápida com `quick.py` (run/compare), A/B na mesma sessão, preset `gcc-bench` genérico (os binários distribuídos são genéricos), preset `msvc-profile` para o uProf. Muda o conteúdo: instâncias de portfólio, contadores (iterações, QPs resolvidos, alocações) e a comparação de **soluções** (hex-float, contagem de saídas idênticas bit a bit) no lugar de log-verossimilhança. A pesquisa levantou hipóteses de hot spot por leitura do código; **nenhuma está medida na biblioteca**, então a Task 1 as confirma ou refuta antes de qualquer otimização.

**Tech Stack:** C++23, Google Benchmark v1.9.1 (mesmo hash do VolEngine), Python (quick.py), CMake presets.

**Spec / pesquisa:** `docs/research/2026-11-16-performance-method.md` (fonte do desenho: seções 2 a 6; leia §1, §3.6 e §5). Spec: `docs/superpowers/specs/2026-11-15-portopt-design.md` §8. Material de origem: `C:/Users/x/Documents/Programs/Cpp/VolEngine/docs/performance/` e `tests/benchmarks/`.

## Global Constraints

- Valem as dos planos 1–5 (estilo do `STYLEGUIDE.md`, commits **sem** `Co-Authored-By`, interface C intocada, `-ffp-contract=off` em `portopt_core` **e em todo alvo novo** que compile o código de produção, sem `long double` em `src/`, gcc-release apenas para os testes).
- **Correção primeiro.** Nenhuma banda é afrouxada para acomodar uma otimização. A suíte inteira, as baterias de fuzz e o `style` ficam verdes em todo commit.
- **Risco 0 = bit a bit idêntico.** Toda otimização desta fase é de risco 0: mantém, para cada elemento do resultado, a sequência de operações de ponto flutuante. Cada uma leva um teste `memcmp` contra a implementação antiga, mantida no teste como referência e compilada com `-ffp-contract=off`, em n variados (inclusive n que não é múltiplo da largura do laço e condição 1e8).
- **Medir antes e depois na mesma sessão**, sempre; nunca comparar contra tabela de outro dia. Números oficiais vêm do build `gcc-bench` genérico (sem `-march=native`). Máquina e build registrados em cada relatório.
- **Fora de escopo (decisão do autor necessária, planos próprios):** kernels que reassociam (multi-acumulador, SIMD de reduções, classe B), multiversão por ISA (`target_clones` falha no MinGW GCC 16.2: o despacho seria manual), warm start do QP no `RISK_PARITY` (muda o caminho), reduzir esforço do `RISK_PARITY`, lote com threads (`po_solve_batch`, spec §5.4). Rejeitadas: BLAS, `-ffast-math`, paralelismo dentro de uma solução, mexer na família implícita de 2ⁿ linhas.
- Regra do STYLEGUIDE para otimização: perfil → A/B → leitura do assembly → número em `docs/performance/` → commit.

## Review Focus

1. **Bits idênticos de verdade.** Um `memcmp` num único n não prova nada: o teste varre n (1, 2, 3, 7, 16, 17, 64, 100), cond ∈ {1e2, 1e8} e dados aleatórios, comparando saída da versão nova com a antiga. (Tasks 2 e 3.)
2. **FMA acidental.** Num build `native-release` (`-mfma` efetivo) nenhuma instrução `vfm*` fora das funções com `std::fma` declarado: `objdump` dos objetos de `portopt_core`. (Task 3.)
3. **A banda fina do `RISK_PARITY`** (cond 1e8, n = 50, sem box, vol-alvo: margem 6) não pode se mexer: a A/B de soluções compara todas as instâncias, e as de risco 0 saem iguais bit a bit. (Tasks 2 e 3.)
4. **Medição honesta.** Os experimentos exploratórios da pesquisa (L1 8,4×, Dekker 1,9×) rodaram numa máquina não identificada e uma vez cada: o plano só os usa como hipótese, o número oficial é o do harness. (Task 1.)
5. **Build genérico × nativo.** O ganho de L6 (Dekker) existe só sem FMA de hardware; medir nos dois e registrar a diferença sem otimizar para o nativo. (Task 3.)

---

### Task 0: Infraestrutura de benchmark e perfil

**Files:**
- Modify: `CMakeLists.txt` (bloco `PO_BUILD_BENCHMARKS`: Google Benchmark v1.9.1 via FetchContent, alvos `portopt_benchmarks` e `portopt_profile`), `CMakePresets.json` (`gcc-bench` genérico sem `-march=native`; `msvc-profile` com PDB)
- Create: `tests/benchmarks/` com `bench_main.cpp`, `instances.hpp` (geradores com cache: Σ de condição escolhida, μ, retornos, cotas), `bm_solve.cpp`, `bm_kernel.cpp`, `profile.cpp`, `quick.py` (`run`, `compare`, `compare-solutions`, `selftest`)
- Modify: `docs/` (`docs/performance/measuring.md` portado e adaptado), `ARCHITECTURE.md`

**Interfaces:**
- Produces:
  - `BM_Solve/<objetivo>/<variante>/<N>` para N ∈ {10, 100, 1000} e os objetivos MIN_VARIANCE, MAX_RETURN, MAX_SHARPE, UTILITY, TRACKING_ERROR, MAX_DIVERSIFICATION, RISK_PARITY (budget, com `target_return`, com `target_volatility`, com box), variantes sem/ com box, sem/com leverage, com caixa; asserção de status em cada caso; contadores `iterations`, `qp_solves`, `allocations` (substituindo `operator new` só no alvo de benchmark);
  - `BM_Kernel/<kernel>/<n>`: `cholesky`, `solve_lower`, `qr`, `qp_iteration` (via `qp::solve`), `compensated_dot`, `covariance_times`, `gauss_newton`, `empty_active_set`;
  - `quick.py compare-solutions A.json B.json`: despejo de soluções em hex-float, imprime quantas saídas são idênticas bit a bit e a pior diferença relativa.
- O filtro da suíte rápida (`QUICK_FILTER`, ~35 s, 3 repetições) cobre um caso de cada objetivo em N = 10 e 100 e os kernels em n = 100.

- [ ] **Step 1:** Testar `quick.py selftest` (compare e compare-solutions em arquivos sintéticos: iguais, 1 ulp de diferença, ordem trocada) e um teste do contador de alocações. Ver falhar.
- [ ] **Step 2:** Implementar a infraestrutura e as instâncias; **Step 3:** rodar `cmake --preset gcc-bench` e a rápida; registrar o tempo da rápida. Sem mudar nada em `src/`.
- [ ] **Step 4: Commit** — `build: benchmark suite and profile tool`

---

### Task 1: Linha de base e confirmação das hipóteses

**Files:**
- Create: `docs/performance/baseline.md`, `docs/performance/benchmarks.md` (porte do formato do VolEngine), `docs/performance/profile-1.md`

- [ ] **Step 1:** Rodar a suíte completa, 5 repetições, `gcc-bench`: tabela por objetivo × N. Registrar máquina, compilador, build, data, "sessão de desktop, não isolada".
- [ ] **Step 2:** Perfil uProf (`msvc-profile`) dos casos lentos (RISK_PARITY com vol-alvo, MAX_RETURN com leverage, N = 1000) e leitura do assembly dos kernels quentes.
- [ ] **Step 3:** Para cada hipótese H1–H10 da pesquisa (§4.1: `empty_active_set` com passo n, `std::fma` como chamada de `libm`, `gauss_newton` e a Cholesky a cada iteração do `RISK_PARITY`, linhas unitárias de bound varridas por candidato, Householder por colunas, Cholesky left-looking, alocações em laços, `Moments` em Armijo, J₀ do caminho do MAX_RETURN): confirmada, refutada ou não medida, com o número. **Reordenar as Tasks 2 a 4 conforme o perfil**: alavanca que o perfil não mostra não entra.
- [ ] **Step 4: Commit** — `docs: performance baseline and hot-spot confirmation`

---

### Task 2: Alavancas de risco 0 do QP e da álgebra linear (L1, L3, L2)

**Files:**
- Modify: `src/qp/active_set.hpp` (`empty_active_set` por linhas, pulando zeros, mesma ordem por elemento), `src/linalg/qr.cpp` (`reflect` por linhas), `src/linalg/cholesky.cpp` (right-looking / com registradores, ordem por elemento preservada)
- Test: `tests/unit/qp/`, `tests/unit/linalg/` (testes `memcmp` com a referência antiga mantida no teste)

- [ ] **Step 1: Testes** `ActiveSet.EmptyBaseIsBitIdenticalToTheOldOrder`, `Qr.ReflectIsBitIdenticalToTheColumnOrder`, `Cholesky.RightLookingIsBitIdenticalToTheLeftLooking`: n ∈ {1, 2, 3, 7, 16, 17, 64, 100}, cond ∈ {1e2, 1e8}, `memcmp` do resultado. Ver os testes passarem contra a implementação antiga (referência) e a cópia.
- [ ] **Step 2:** Uma alavanca por vez: implementar, `memcmp` verde, suíte inteira verde, A/B na rápida (`compare`) e de soluções (`compare-solutions`: todas as saídas idênticas bit a bit). **Step 3:** Registrar o ganho medido em `docs/performance/`.
- [ ] **Step 4: Commits** — um por alavanca: `perf(qp): row-oriented empty active set base`, `perf(linalg): row-oriented Householder reflect`, `perf(linalg): right-looking Cholesky`

---

### Task 3: Alavancas de risco 0 dos problemas e do produto compensado (L4, L5, L7, L6a)

**Files:**
- Modify: `src/riskparity/objective.cpp` (`gauss_newton` por simetria e blocagem; `compensated_dot` com TwoProduct de Dekker/Veltkamp em vez de `std::fma`, mantendo `-ffp-contract=off`), `src/qp/` (fast path de linha unitária de bound: sinalizador `row_is_unit` em `gather_constraints`), `src/problems/max_return.cpp`/`parametric_path.cpp` (cache de J₀ entre solves com o mesmo fator)
- Test: testes `memcmp` por alavanca; teste de bordas do Dekker (overflow de Veltkamp, subnormais) contra `std::fma`; `objdump` do `portopt_core` no `native-release` sem `vfm*` fora de `std::fma`

- [ ] **Step 1: Testes** como na Task 2, mais `CompensatedDot.DekkerMatchesFmaOnEdgeCases` (valores grandes, subnormais, zeros com sinal, 2000 casos aleatórios) e o checklist de FMA acidental (Review Focus 2).
- [ ] **Step 2:** Uma por vez, na ordem que o perfil da Task 1 fixar; cada uma com `memcmp`, suíte, A/B de tempo e de soluções. **Step 3:** Registrar ganhos e números do Dekker em build genérico e nativo (Review Focus 5).
- [ ] **Step 4: Commits** — `perf(riskparity): gauss_newton by symmetry`, `perf(qp): unit-row fast path`, `perf(problems): cache J0 along the MAX_RETURN path`, `perf: Dekker TwoProduct in the compensated dot`

---

### Task 4: Higiene de alocação, reperfil e decisões de fase 2

**Files:**
- Modify: `src/qp/goldfarb_idnani.cpp`, `src/problems/*.cpp` (reuso de buffers no laço) **só se o contador de alocações da Task 0 mostrar**; `src/riskparity/objective.cpp` (cache de `Moments` na busca de Armijo) opcional
- Create: `docs/performance/profile-2.md`

- [ ] **Step 1:** Medir o contador de alocações por solve em N = 10 e 100; aplicar L8/L9 só onde o perfil mostrar, com `memcmp`.
- [ ] **Step 2:** Reperfilar (mesma rotina da Task 1) e registrar. **Decidir, com os números, se vale** L10 (`-flto`, só se `-ffp-contract=off` valer no link: verificar), e levar ao autor as alavancas de classe B, a multiversão por ISA e o warm start do `RISK_PARITY` como planos próprios.
- [ ] **Step 3: Commit** — `perf: allocation hygiene` (se houver mudança) e `docs: second profile and phase-2 decisions`

---

### Task 5: Camada C e Python, relatório final e portão da v1

**Files:**
- Create: `tests/benchmarks/bm_capi.cpp`, `python/benchmarks/bench_overhead.py`, `docs/performance/python-vs-arch.md` (modelo do VolEngine)
- Modify: `docs/performance/benchmarks.md` (tabela antes/depois), `README.md`, `CHANGELOG.md`

- [ ] **Step 1:** Medir o custo fixo da chamada C e do ctypes (criar handles, copiar dados, `po_solve`, ler o resultado) contra o tempo do solve em N = 10, 100, 1000.
- [ ] **Step 2:** Relatório final: tabela antes/depois por objetivo × N, o que ficou de fora e por quê, `msvc-release` verde (compilar e rodar a bateria; só aqui, não a cada task), rodada `native-release` apenas para registrar a diferença do build nativo.
- [ ] **Step 3:** Verificar: `ctest --preset gcc-release` inteiro, pytest, mypy, `check_style.py`, `compare-solutions` entre o commit anterior ao plano e o final (todas as saídas de risco 0 idênticas bit a bit).
- [ ] **Step 4: Commit** — `docs: final performance report`
