# PortOpt — Plano 7: fechamento da v0.x (API, wrapper C++, wheels, release)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Levar a PortOpt até `pip install pyportopt` nas sete plataformas do spec §11, numa release 0.x reproduzível por tag, com a API revisada e o wrapper C++ header-only como segundo consumidor.

**Architecture:** Primeiro a análise de API (documento, parada para o autor) e as correções aceitas; depois o contrato de ABI, o wrapper C++ e as regras de `install` do CMake; então o repositório remoto passa a ser o C++ e o pipeline do VolEngine é portado (`wheels.yml` com containers manylinux/musllinux, `release.yml` com import library, dump de ABI, SBOM, GitHub Release e PyPI por trusted publishing). Nada muda um bit de solução.

**Tech Stack:** C++23, CMake 3.28 (`install(EXPORT)`, `CMakePackageConfigHelpers`), GoogleTest v1.17.0, Python ≥ 3.11 com hatchling e uv, GitHub Actions, Docker (manylinux_2_28, musllinux_1_2), libabigail.

**Spec:** `docs/superpowers/specs/2026-11-17-portopt-v0-release-design.md` (este plano), sobre `docs/superpowers/specs/2026-11-15-portopt-design.md` §5–7, §11.

**Material de origem:** o pipeline do VolEngine está **só no remoto** (`https://github.com/GabrielMarquezMatte/VolEngine`, tag `v1.0.0`); o checkout em `../VolEngine` é anterior a ele. Clonar no scratchpad: `git clone --depth 1 https://github.com/GabrielMarquezMatte/VolEngine <scratchpad>/ve`. "Portar X" neste plano = copiar o arquivo do VolEngine e trocar `volengine`→`portopt`, `ve_`/`VE_`→`po_`/`PO_`, `VE_BUILD_TESTS`→`PO_BUILD_TESTS`, `VE_VERSION`→`PO_VERSION`, `VOLENGINE_LIB`→`PORTOPT_LIB`, `volengine-native-`→`portopt-native-`, mais o que a task disser. Comentários que justificam algo do VolEngine (Highway, DCC, threads) são reescritos para o motivo da PortOpt ou removidos.

## Global Constraints

- Valem as dos planos 1–6: estilo do `STYLEGUIDE.md`, commits **sem** `Co-Authored-By`, `-ffp-contract=off` em todo alvo que compila código de produção, sem `long double` em `src/`, documentos de `docs/superpowers/` em português, documentos do repositório em inglês.
- **Nenhum bit de solução muda.** Correção de API mexe em nomes, validação e erros. Se tocar o núcleo: `portopt_profile --max-n 100 --json` antes e depois e `quick.py compare-solutions A B --require-identical` (docs/performance/measuring.md) verde, ou sai do plano.
- 0.x: `PO_ABI_VERSION` fica 0. Toda quebra entra no CHANGELOG como "antes → depois".
- Bibliotecas distribuídas: genéricas (sem `-march`/`-mcpu`/`-mtune`/`-mavx`); MinGW liga `libstdc++`, `libgcc`, `winpthread` estaticamente; Linux liga `libstdc++`/`libgcc` estaticamente (já no `CMakeLists.txt`).
- Linux: **GCC 15+** dentro dos containers `quay.io/pypa/manylinux_2_28_{x86_64,aarch64}` e `quay.io/pypa/musllinux_1_2_{x86_64,aarch64}`.
- macOS: `MACOSX_DEPLOYMENT_TARGET=11.0` em arm64; x64 medido (Task 8) e registrado.
- Nome de distribuição **`pyportopt`**; import `portopt`, `libportopt.*`, `find_package(PortOpt)` não mudam.
- Ações do GitHub fixadas por major (`actions/checkout@v4`, `astral-sh/setup-uv@v6`, `msys2/setup-msys2@v2`, `ilammy/msvc-dev-cmd@v1`, `actions/upload-artifact@v4`, `actions/download-artifact@v4`, `anchore/sbom-action@v0`, `softprops/action-gh-release@v2`).
- Assets: `portopt-native-<plataforma>.{dll,so,dylib}`, `portopt-native-win-x64.{dll.a,lib}`, `portopt-abi-linux-x64.xml`, `portopt.spdx.json`; plataformas `linux-x64`, `linux-arm64`, `linux-musl-x64`, `linux-musl-arm64`, `osx-arm64`, `osx-x64`, `win-x64`.
- **Anti-objetivo:** publicar sem o autor. Force-push, projeto/trusted publisher no PyPI e push da tag são do autor; o plano para e pede em cada um.

## Review Focus

1. **Container sem GCC 15.** Se `manylinux_2_28` não tiver `gcc-toolset-15` ou o Alpine do `musllinux_1_2` vier com GCC < 15, o build tem de falhar dizendo isso, não compilar com o GCC antigo nem gerar wheel com tag errada. (Task 8: `build_native.py` confere `gcc -dumpversion`.)
2. **Python 3.11, o piso do `requires-python`.** A wheel é `py3-none`; uma sintaxe de 3.12+ quebraria quem está em 3.11. (Task 9: pytest em 3.11 no `ci.yml`.)
3. **Usuário MSVC com a DLL MinGW da release.** O `.lib` gerado do `.def` tem de ligar e rodar um programa C++ compilado pelo `cl`. (Task 10: o job `windows-lib` compila e roda `cpp/tests/package/main.cpp` contra o `.lib` e a DLL.)
4. **Uso errado do wrapper.** `Data` movido e usado depois, `std::span` de tamanho errado: erro `invalid_handle` / `invalid_argument` com mensagem, nunca leitura fora do buffer. (Task 4.)
5. **Tag e versões desalinhadas.** `v0.2.0` com CMake em 0.1.0, CHANGELOG sem a seção, ou tag `v0.1.0rc1`/`v0.1.0-rc.1`: a primeira falha antes de compilar; as duas formas de rc dão `0.1.0rc1`. (Task 6.)

---

### Task 1: Análise de API (documento) — PONTO DE PARADA

**Files:**
- Create: `docs/guide/api-review.md` (inglês)

**Interfaces:**
- Produces: a lista de decisões numerada (`R1`, `R2`, …) que a Task 2 executa e que a Task 4 aplica aos nomes do wrapper.

- [ ] **Step 1:** Escrever o documento: uma tabela por superfície — interface C (cada `PO_API` e cada `#define` de `include/portopt/portopt.h`), Python (cada nome de `portopt.__all__`, seus métodos e campos), wrapper C++ (os nomes **planejados** na Task 4 deste plano). Colunas: `Item | Finding | Decision (keep / change / remove / add) | Breaks | Migration (before → after)`. Critérios: spec §5.1 (setters, enums `int32_t` + `#define`, `po_status` + `po_last_error`, handles independentes, tipos na fronteira), anti-objetivo "API inchada", consistência de nomes entre C, Python e C++, mensagens de erro, documentação no header, usos do spec §6.
- [ ] **Step 2:** Conferir, entre outros, estes pontos de partida (achados a verificar, não decisões): `po_version` não está em `PROTOTYPES` de `_native.py` e não há `portopt.__version__`; `PO_INVALID_HANDLE` e `PO_OUT_OF_MEMORY` não têm nome em `_native.py`; `Data.volatility` (Python) × `po_portfolio_volatility` (C); `po_result_iterations` é `int32_t` e as contagens são `size_t`; Python não expõe `asset_count`/`weight_count`; os nomes dos parâmetros no header (`lambda`, `delta`, `sigma`, `max_gross`) × os campos de `Spec`.
- [ ] **Step 3:** Para cada `change`/`remove`/`add`, escrever o "antes → depois" e se toca `src/core` ou além (se tocar, marcar "exige compare-solutions").
- [ ] **Step 4: Commit** — `docs: API review of the C, Python and C++ surfaces`
- [ ] **Step 5: PARAR.** Mostrar ao autor a lista de decisões e esperar a aprovação item a item. A Task 2 só implementa os itens aprovados; o documento é atualizado com a decisão final de cada um.

---

### Task 2: Correções de API aprovadas

**Files:** por item (header, `src/capi/*`, `python/src/portopt/*`, `tests/unit/capi/*`, `python/tests/*`, `CHANGELOG.md`, `docs/guide/api-review.md`).

**Interfaces:**
- Consumes: a lista aprovada da Task 1.
- Produces: o header e o pacote Python finais da 0.1.0, que as Tasks 3 e 4 tomam como base.

Para **cada** item aprovado, um commit:

- [ ] **Step 1:** Teste que falha: em `tests/unit/capi/` para o C (nome do teste = o comportamento novo, ex. `Version.PythonExposesTheLibraryVersion`), em `python/tests/` para o Python. Rodar e ver falhar.
- [ ] **Step 2:** Implementar; atualizar o comentário do header, os `PROTOTYPES` de `_native.py` e o `__all__` se o nome é público.
- [ ] **Step 3:** `ctest --preset native-release` e `uv run --project python --python 3.14 pytest python/tests -q` e mypy (comando do `ci.yml`) verdes. Se o item tocou o núcleo: `compare-solutions --require-identical` verde.
- [ ] **Step 4:** CHANGELOG em `[Unreleased]` (`Changed`/`Removed`/`Added`) com "antes → depois"; marcar o item como feito em `api-review.md`.
- [ ] **Step 5: Commit** — `feat(api)!: <item>` quando quebra, `feat(api): <item>` quando só acrescenta.

---

### Task 3: Contrato de ABI

**Files:**
- Create: `scripts/abi_contract.py` (porte), `include/portopt/portopt.abi` (gerado), `tests/abi/test_abi_contract.py` (porte, nomes trocados)
- Modify: `CMakeLists.txt` (teste `abi_contract` no bloco `PO_BUILD_TESTS`), `.github/workflows/ci.yml`

**Interfaces:**
- Produces:
  - `python scripts/abi_contract.py write | check | def --library NAME OUT | exports DEF` (sem `base-ref`/`--base` em uso até a 1.0; mantidos no script);
  - `include/portopt/portopt.abi`: uma linha por `#define PO_*`, `typedef`, `function RET po_name(TYPES)`;
  - `def --library libportopt.dll OUT` escreve a lista de exports do Windows (Tasks 8 e 10).

- [ ] **Step 1:** Portar `tests/abi/test_abi_contract.py` (os 17 testes do VolEngine com `po_`/`PO_API`; `test_header_matches_committed_contract` lê o header e o `.abi` reais). Registrar em CMake: `add_test(NAME abi_contract COMMAND Python3::Interpreter -m unittest discover -s ${CMAKE_SOURCE_DIR}/tests/abi -t ${CMAKE_SOURCE_DIR})`. Rodar `ctest -R abi_contract` e ver falhar (script ausente).
- [ ] **Step 2:** Portar `scripts/abi_contract.py` (`PLATFORM_MACROS = {"PO_API", "PO_BUILDING"}`, `--library` padrão `libportopt.dll`); `python scripts/abi_contract.py write`. Conferir que o `.abi` tem as 36 funções (ou o número após a Task 2) e nenhum erro de leitura do header.
- [ ] **Step 3:** `ctest --preset native-release` verde. Em `ci.yml` (job `build-test`, Linux): passo `python scripts/abi_contract.py check`.
- [ ] **Step 4: Commit** — `build: ABI contract of portopt.h`

---

### Task 4: Wrapper C++23 header-only

**Files:**
- Create: `cpp/include/portopt/portopt.hpp`, `cpp/tests/CMakeLists.txt`, `cpp/tests/wrapper_test.cpp`, `cpp/scripts/check_abi_coverage.py` (porte), `cpp/README.md` (curto: o que é, como ligar)
- Modify: `CMakeLists.txt` (alvo `portopt_cpp`; `add_subdirectory(cpp/tests)` sob `PO_BUILD_TESTS`)

**Interfaces:**
- Consumes: `include/portopt/portopt.h` após a Task 2 (aplicar aos nomes abaixo as renomeações aprovadas).
- Produces (namespace `po`, nomes iguais aos do Python, em `snake_case`):
  - `enum class Status : std::int32_t { ok, invalid_argument, invalid_handle, unsupported, out_of_memory, internal_error }` (valores = `PO_*`);
  - `struct Error { Status status; std::string message; }`; `template <class T> using Expected = std::expected<T, Error>;`
  - `enum class Objective : std::int32_t { min_variance, max_return, max_sharpe, utility, importance, kelly, tracking_error, max_diversification, risk_parity }`; `enum class SolveStatus : std::int32_t { optimal, infeasible, unbounded, max_iter }` (valores = `PO_*`);
  - `struct Bounds { double lower; double upper; }`; `struct AssetBounds { std::vector<double> lower; std::vector<double> upper; }`;
  - `struct Spec { Objective objective; std::variant<std::monostate, Bounds, AssetBounds> bounds; std::optional<double> target_return, risk_aversion, importance, max_volatility, leverage, target_volatility; std::vector<double> risk_budget; }` — valor, com designated initializers; vetor vazio = padrão `1/n`;
  - `class Data` (RAII sobre `po_data*`, movível, não copiável): `static Expected<Data> create(std::size_t assets, std::span<const double> mu, std::span<const double> cov)` (span vazio = `NULL`); `Expected<void> set_risk_free(double)`; `Expected<void> set_returns(std::span<const double> returns, std::span<const double> benchmark)`; `Expected<std::size_t> asset_count() const`, `weight_count() const`; métricas `Expected<double> portfolio_return(std::span<const double> w) const`, `volatility`, `utility(w, double risk_aversion)`, `tracking_error`, `sharpe_ratio`, `diversification_ratio`; `Expected<std::vector<double>> risk_contributions(w) const`; `const po_data* get() const noexcept`;
  - `struct Result { SolveStatus status; std::vector<double> weights; std::int32_t iterations; double budget_multiplier, target_return_multiplier, volatility_multiplier, leverage_multiplier; std::vector<double> bound_multipliers; double residual; }`;
  - `Expected<Result> solve(const Spec&, const Data&)` — cria o `po_spec` num handle RAII interno, aplica os campos presentes, chama `po_solve`, lê tudo e libera;
  - `std::string_view version()`; `Data::create` confere `po_abi_version() == PO_ABI_VERSION` e devolve `Status::unsupported` com `"the PortOpt library has ABI version X, but these headers need version Y"`.
- CMake: `add_library(portopt_cpp INTERFACE)` com `cpp/include` (`$<BUILD_INTERFACE>`/`$<INSTALL_INTERFACE:include>`), `cxx_std_23`, liga `portopt`; aliases `PortOpt::portopt` e `PortOpt::cpp` (`EXPORT_NAME cpp`).

- [ ] **Step 1: Testes que falham** em `cpp/tests/wrapper_test.cpp` (GoogleTest, `/W4 /WX` ou `-Wall -Wextra -Wpedantic -Werror`, ligado a `PortOpt::cpp` — a biblioteca **compartilhada**, não `portopt_capi`):
  - `Wrapper.MinVarianceMatchesTheCApiBitForBit`: n = 3, Σ fixa, `bounds = Bounds{0, 1}`; os pesos e os quatro multiplicadores do wrapper `memcmp`-iguais aos de chamar `po_*` direto.
  - `Wrapper.ErrorCarriesStatusAndMessage`: `Spec{.objective = Objective::min_variance, .risk_aversion = 2.0}` → `status == invalid_argument`, `message == "risk_aversion does not apply to MIN_VARIANCE"` (a mensagem atual de `po_solve`; após a Task 2, a que valer).
  - `Wrapper.SpanSizeMismatchIsInvalidArgument`: `Data::create(3, mu de 2, cov de 9)` → `invalid_argument`, mensagem cita `mu` e os tamanhos; idem `volatility` com `w` de tamanho ≠ `weight_count`, e `set_returns` com `returns.size() % assets != 0` ou `benchmark.size() != returns.size() / assets`.
  - `Wrapper.MovedFromDataIsInvalidHandle`: após `Data b = std::move(a)`, `a.volatility(w)` → `invalid_handle` (o C recebe NULL); `b` segue funcionando.
  - `Wrapper.InfeasibleIsAResultNotAnError`: `target_return` acima de todo μ com `Bounds{0, 1}` → `Expected` com valor, `status == SolveStatus::infeasible`.
  - `Wrapper.AbiMismatchIsUnsupported`: o helper `detail::check_abi_version(1, 0)` → `unsupported` com a mensagem acima.
  - `static_assert(!std::is_copy_constructible_v<po::Data> && std::is_nothrow_move_constructible_v<po::Data>)`.
  - Em `cpp/tests/CMakeLists.txt`: copiar `$<TARGET_RUNTIME_DLLS:portopt_cpp_tests>` para junto do executável **antes** de `gtest_discover_tests` (Windows); teste `abi_coverage` rodando `check_abi_coverage.py`.
  - Rodar `ctest --preset native-release -R "Wrapper|abi_coverage"` e ver falhar.
- [ ] **Step 2:** Implementar `portopt.hpp` (um arquivo; o wrapper faz só tradução: conferência de tamanhos dos spans, chamadas, cópia da mensagem de `po_last_error()` no instante do erro). Portar `check_abi_coverage.py` (`PO_API`, `po_\w+`, wrapper em `cpp/include/portopt`); `UNUSED_ALLOWED` vazio, ou cada exceção com o motivo.
- [ ] **Step 3:** `ctest --preset native-release` verde (MSVC local) e `ctest --preset gcc-release` verde (MinGW).
- [ ] **Step 4: Commit** — `feat(cpp): header-only C++23 wrapper over the C API`

---

### Task 5: `install`, `find_package(PortOpt)` e `cpp.yml`

**Files:**
- Modify: `CMakeLists.txt` (regras de install), `include` dirs de `portopt` para `$<BUILD_INTERFACE>`/`$<INSTALL_INTERFACE>`
- Create: `cmake/PortOptConfig.cmake.in`, `cpp/tests/package/CMakeLists.txt`, `cpp/tests/package/main.cpp`, `cpp/scripts/ci_check.sh`, `.github/workflows/cpp.yml`

**Interfaces:**
- Consumes: `portopt`, `portopt_cpp` (Task 4).
- Produces:
  - `option(PO_INSTALL ... ${PO_IS_TOP_LEVEL})`; `install(TARGETS portopt portopt_cpp EXPORT PortOptTargets)` (RUNTIME → `bin`, LIBRARY/ARCHIVE → `lib`, inclui `libportopt.dll.a` no MinGW), `install(FILES include/portopt/portopt.h cpp/include/portopt/portopt.hpp DESTINATION include/portopt)`, `install(EXPORT PortOptTargets NAMESPACE PortOpt:: DESTINATION lib/cmake/PortOpt)`, `PortOptConfig.cmake` (inclui os targets + `check_required_components`) e `PortOptConfigVersion.cmake` com `COMPATIBILITY SameMinorVersion` (0.x quebra entre menores);
  - consumidor: `cmake -S cpp/tests/package -DPORTOPT_CONSUMER_MODE=add_subdirectory|find_package [-DPORTOPT_CONSUMER_VERSION=X] [-DCMAKE_PREFIX_PATH=prefix]`, executável `package_consumer` ligado a `PortOpt::cpp`, DLLs copiadas por `$<TARGET_RUNTIME_DLLS>`, um `add_test`;
  - `main.cpp`: MAX_SHARPE de 3 ativos com `risk_free`, confere `status == optimal`, soma dos pesos = 1 a 1e-12, imprime `version()` e os pesos; sai 1 com a mensagem do `Error` em qualquer falha.

- [ ] **Step 1:** Escrever `cpp/tests/package/` e `cpp/scripts/ci_check.sh` (porte: configura e instala a raiz com `PO_BUILD_TESTS=OFF` em `build/cpp-prefix`, roda o consumidor nos dois modos, e um terceiro configure pedindo `PORTOPT_CONSUMER_VERSION=0.2` que **tem de falhar**). Rodar e ver falhar (sem install).
- [ ] **Step 2:** Implementar as regras de install e o config. `bash cpp/scripts/ci_check.sh` verde com MinGW (`gcc-release`) e com MSVC (Developer prompt, `CPP_CMAKE_ARGS=-DCMAKE_CXX_COMPILER=cl`).
- [ ] **Step 3:** `cpp.yml` (porte, gatilhos em `cpp/**`, `include/**`, `src/**`, `CMakeLists.txt`, `cmake/**`, o próprio arquivo): `linux-gcc` (GCC 15 do PPA, como o `ci.yml`), `macos-clang`, `windows-msvc` rodam `ci_check.sh`; `linux-gcc` roda também os testes do wrapper sob `-fsanitize=address,undefined -fno-sanitize-recover=all` (configure da raiz com `-DCMAKE_BUILD_TYPE=Debug`, `cmake --build ... --target portopt_cpp_tests`, `ctest -R Wrapper`; só o alvo do wrapper, não as baterias de fuzz). Validar com `actionlint` se disponível; a execução real fica para a Task 7.
- [ ] **Step 4: Commit** — `build: install rules, find_package(PortOpt) and the C++ workflow`

---

### Task 6: Versão única e nome `pyportopt`

**Files:**
- Create: `python/scripts/release_version.py`, `python/tests/packaging/test_release_version.py`
- Modify: `python/pyproject.toml` (`name = "pyportopt"`), `python/uv.lock` (`uv lock`), `.github/workflows/ci.yml`

**Interfaces:**
- Produces:
  - `python python/scripts/release_version.py TAG` → imprime a versão PEP 440 do wheel e a `x.y.z` do CMake (porte);
  - `python python/scripts/release_version.py check [TAG]` → sai 1 com a diferença se `set(PO_VERSION X ...)` do `CMakeLists.txt` ≠ `version = "X"` do `pyproject.toml`; com `TAG`, também se `base_version(TAG) ≠ X` ou se o CHANGELOG não tem `## [X] - AAAA-MM-DD`.
  - funções testáveis: `python_version(tag) -> str`, `cmake_version(tag) -> str`, `mismatches(cmake: str, pyproject: str, changelog: str, tag: str | None) -> list[str]` (recebem o texto dos arquivos).

- [ ] **Step 1: Testes que falham:** `test_plain_tag` (`v0.1.0` → `0.1.0`, `0.1.0`), `test_prerelease_tag` (`v0.1.0-rc.1` e `v0.1.0rc1` → `0.1.0rc1`, `0.1.0`), `test_tag_without_v_is_rejected`, `test_cmake_and_pyproject_must_agree`, `test_tag_must_match_the_cmake_version` (`v0.2.0` × `0.1.0` → uma linha citando os dois), `test_tag_needs_a_dated_changelog_section` (`## [Unreleased]` só → falha; `## [0.1.0] - 2026-11-20` → passa), `test_real_repository_agrees` (`check` sem tag no repositório real).
- [ ] **Step 2:** Implementar; renomear a distribuição para `pyportopt`; `uv lock`; conferir que `uv build --wheel python` gera `pyportopt-0.1.0-py3-none-<tag>.whl` e que o import continua `portopt`.
- [ ] **Step 3:** `ci.yml`: passo `python python/scripts/release_version.py check` no job Linux. `pytest python/tests` verde.
- [ ] **Step 4: Commit** — `build: single version check and the pyportopt distribution name`

---

### Task 7: Troca do repositório remoto — AÇÃO DO AUTOR

Remoto: `origin` = `https://github.com/GabrielMarquezMatte/PortfolioOptimization`, `master` remoto em `2b16551` (R, 2022). Nada aqui roda sem o ok explícito do autor **no momento**.

- [ ] **Step 1:** `git fetch origin` e conferir `git rev-parse origin/master` = `2b16551…`. Se for outro commit, **parar** e mostrar ao autor.
- [ ] **Step 2:** Preservar o R (não destrutivo): `git push origin origin/master:refs/heads/r-legacy`, `git tag r-final origin/master`, `git push origin r-final`.
- [ ] **Step 3:** `git ls-remote origin refs/heads/r-legacy refs/tags/r-final` mostra os dois em `2b16551…`. Mostrar a saída ao autor.
- [ ] **Step 4:** Pedir ao autor, citando o comando: `git push --force-with-lease=master:2b16551 origin master`. Só com o "sim" dele, rodar.
- [ ] **Step 5:** Conferir no GitHub: Actions habilitadas, `ci.yml` e `cpp.yml` disparados pelo push; esperar os dois verdes (`gh run watch`). Vermelho → corrigir em commit novo antes da Task 8.

---

### Task 8: Wheels das sete plataformas

**Files:**
- Create: `python/scripts/build_native.py`, `build_wheel.sh`, `build_in_container.sh`, `check_wheel_tag.py`, `check_native_deps.py` (portes); `python/tests/packaging/test_packaging_scripts.py`, `tests/abi/test_native_exports.py` (portes); `.github/workflows/wheels.yml` (porte)
- Modify: `.gitignore` (`dist/`, `out/`, `python/src/portopt/_lib/`), `CMakeLists.txt` (`add_test` de `tests/abi` já cobre `test_native_exports.py`)

**Interfaces:**
- Consumes: `release_version.py` (Task 6), `abi_contract.py def|exports` (Task 3), `hatch_build.py` como está.
- Produces: artefatos `wheel-<plataforma>` (`dist/*.whl`) e `native-<plataforma>` (`out/portopt-native-<plataforma>.*`, mais `.dll.a` e `.def` no Windows); `wheels.yml` com `workflow_call` input `version`.

Diferenças do porte:
- `build_native.py`: alvo `portopt`, `-DPO_BUILD_TESTS=OFF`, `LIBRARIES = ("libportopt.so", "libportopt.dylib", "libportopt.dll", "libportopt.dll.a", "portopt.def")`, motivo do `check_flags`: "the distributed binaries are generic (no -march), like gcc-bench". Novo `check_gcc(version: str) -> None`: em Linux e Windows, roda `gcc -dumpversion` e sai com `error: GCC 15 or newer is required, found <versão>` se o major for < 15.
- `build_in_container.sh`: em manylinux, `. /opt/rh/gcc-toolset-15/enable` se existir (no lugar do 14); em musllinux, `apk add --no-cache bash g++`. A conferência da versão fica no `check_gcc` acima.
- `check_native_deps.py`: prefixo `po_`.
- `build_wheel.sh`: `-DPO_VERSION=`, `unset PORTOPT_LIB`, assert `'site-packages' in str(portopt._native.LIBRARY_DIRECTORY)`, `.def` via `abi_contract.py def --library libportopt.dll`.
- `wheels.yml`: gatilhos de path como o VolEngine mais `cmake/**`.

- [ ] **Step 1:** Portar os testes dos scripts (14 em `test_packaging_scripts.py`, 5 em `test_native_exports.py`, com `po_`) e mais `test_build_refuses_gcc_older_than_15` (`check_gcc("14.2.0")` → `SystemExit` com a mensagem; `check_gcc("15.1.0")` e `check_gcc("16")` passam). Ver falhar.
- [ ] **Step 2:** Portar os scripts. Local: `bash python/scripts/build_wheel.sh` no MSYS2 UCRT64 com `ASSET=win-x64 EXPECTED_TAG=win_amd64` → wheel `pyportopt-0.1.0-py3-none-win_amd64.whl`, pytest instalado verde, `out/` com `.dll`, `.dll.a`, `.def`.
- [ ] **Step 3:** Local com Docker: `docker run --rm -v "$PWD:/work" -w /work -e ASSET=linux-x64 -e EXPECTED_TAG=manylinux_2_28_x86_64 -e VERSION quay.io/pypa/manylinux_2_28_x86_64 sh python/scripts/build_in_container.sh` e o mesmo com `musllinux_1_2_x86_64`/`linux-musl-x64`. **Se algum container não tiver GCC 15, parar e levar ao autor** (opções: outra imagem, com outra tag de plataforma; ou pedir a exceção de GCC 14 só para o container) — não é decisão do builder.
- [ ] **Step 4: Commit** — `build: wheel pipeline for the seven platforms`
- [ ] **Step 5:** Push; `wheels.yml` verde nas sete entradas (`gh run watch`). macOS x64: conferir que `MACOSX_DEPLOYMENT_TARGET=11.0` compila; se a libc++ recusar, subir o mínimo até compilar e registrar o valor e o motivo em `CONTRIBUTING.md` e no `expected_tag`. Registrar em `CONTRIBUTING.md` o GCC de cada container. Commit `docs: platform floors of the wheels` se houve mudança.

---

### Task 9: CI completa (benchmarks e Python 3.11)

**Files:**
- Modify: `.github/workflows/ci.yml`

- [ ] **Step 1:** Job `build-test`, Linux: passo "Build benchmarks" — `cmake -S . -B build/bench -G Ninja -DCMAKE_BUILD_TYPE=Release -DPO_BUILD_BENCHMARKS=ON -DPO_BUILD_TESTS=OFF` (CC/CXX = GCC 15 do ambiente), `cmake --build build/bench`, `ctest --test-dir build/bench --output-on-failure` (os dois self-tests: `benchmark_allocation_selftest`, `benchmark_quick_selftest`). Comentário: compila e confere, nunca mede (runner compartilhado).
- [ ] **Step 2:** Passo "Python tests (3.11, the requires-python floor)" com `--python 3.11`, nas três plataformas.
- [ ] **Step 3:** Push; `ci.yml` verde. **Commit** — `ci: build the benchmarks and test Python 3.11`

---

### Task 10: `release.yml` ensaiado em `dry_run`

**Files:**
- Create: `.github/workflows/release.yml` (porte), `scripts/abi_dump.sh`, `scripts/install_abi_tools.sh` (portes)

**Interfaces:**
- Consumes: `wheels.yml` (`workflow_call`, input `version`), `release_version.py check TAG`, artefatos `native-*`, `wheel-*`.

Diferenças do porte:
- Job `version` (primeiro, Linux, sem compilar): em tag, `python python/scripts/release_version.py check "$GITHUB_REF_NAME"`; `wheels` depende dele.
- `windows-lib`: `lib /def:portopt-native-win-x64.def /name:libportopt.dll /machine:x64 /out:portopt-native-win-x64.lib`, e depois o ensaio do Review Focus 3: copiar a DLL como `libportopt.dll`, `cl /std:c++latest /EHsc /W4 /WX /I include /I cpp/include cpp/tests/package/main.cpp portopt-native-win-x64.lib`, rodar `main.exe` (sai 0). Precisa de checkout.
- `abi_dump.sh`: `-DPO_BUILD_TESTS=OFF`, alvo `portopt`, `abidw --headers-dir include/portopt ... build/abi/libportopt.so`, sanidade `grep -c 'po_solve'`.
- `publish-pypi`: `environment: ${{ contains(github.ref_name, '-') && 'testpypi' || 'pypi' }}`; para `testpypi`, `uv publish --publish-url https://test.pypi.org/legacy/ --trusted-publishing always dist/*`; para `pypi`, sem `--publish-url`. Os dois passos de publicação (GitHub Release e PyPI) com `if: github.event_name == 'push' || !inputs.dry_run`.
- SBOM: `portopt.spdx.json`; `files:` com os nomes `portopt-*`.

- [ ] **Step 1:** Portar e ajustar. Rodar `bash scripts/abi_dump.sh build/abi-out/current.xml` num container `ubuntu:24.04` local (`install_abi_tools.sh` primeiro) e conferir o XML.
- [ ] **Step 2: Commit** — `ci: release workflow`; push.
- [ ] **Step 3:** `gh workflow run release.yml -f dry_run=true`; `gh run watch` até verde. Conferir nos artefatos: 7 wheels, 7 bibliotecas, `.dll.a`, `.lib`, XML de ABI, SBOM; o passo de ensaio MSVC rodou e saiu 0; os passos de publicação pulados.

---

### Task 11: README, CHANGELOG e preparação da tag

**Files:**
- Modify: `README.md`, `CHANGELOG.md`, `ARCHITECTURE.md` (seção curta: pipeline de release e contrato de ABI), `CONTRIBUTING.md` (como fazer uma release)

- [ ] **Step 1:** README: instalação (`pip install pyportopt`, `import portopt`; C/C++ pela GitHub Release — `.so`/`.dylib`/`.dll` + `.lib` para MSVC + `portopt.h` — ou por `cmake --install` e `find_package(PortOpt 0.1 REQUIRED)` com `PortOpt::cpp`), um exemplo Python (o do spec §6) e um C++ (o `main.cpp` da Task 5). Conferir que os dois rodam como estão no README (colar e executar).
- [ ] **Step 2:** CHANGELOG: `[Unreleased]` vira `## [0.1.0] - AAAA-MM-DD` (data do dia da tag), com as entradas deste plano (wrapper C++, install, wheels, `pyportopt`, as quebras da Task 2 com "antes → depois"); `## [Unreleased]` vazio acima; links de comparação no fim.
- [ ] **Step 3:** `python python/scripts/release_version.py check v0.1.0-rc.1` e `... check v0.1.0` saem 0. **Commit** — `docs: README and CHANGELOG for 0.1.0`; push; `ci.yml`, `cpp.yml`, `wheels.yml` verdes.
- [ ] **Step 4: PARAR e entregar ao autor** a lista do spec §8, nesta ordem: (1) registrar o pending trusted publisher de `pyportopt` no TestPyPI e no PyPI (repositório `GabrielMarquezMatte/PortfolioOptimization`, workflow `release.yml`, environments `testpypi` e `pypi`) e criar os dois environments no GitHub; (2) `git tag v0.1.0-rc.1 && git push origin v0.1.0-rc.1`; (3) depois do rc verde, conferir num venv limpo: `uv pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ pyportopt==0.1.0rc1` e `python -c "import portopt; print(portopt.solve)"`; (4) `git tag v0.1.0 && git push origin v0.1.0`.
