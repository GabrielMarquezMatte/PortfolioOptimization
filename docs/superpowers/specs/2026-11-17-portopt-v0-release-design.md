# PortOpt — design do fechamento da v0.x (plano 7)

Data: 2026-11-17. Complementa `2026-11-15-portopt-design.md` (o spec da biblioteca), seções 5 a 7 e 11.

## 1. Intenção

Pôr a PortOpt nas mãos de quem usa Python com `pip install pyportopt`, em todas as plataformas do spec §11, numa
release 0.x reproduzível, e deixar o wrapper C++ pronto como segundo consumidor. A interface C ainda **não** congela:
é 0.x, `PO_ABI_VERSION` continua 0 e pode quebrar entre releases (CHANGELOG já diz isso). A análise de API existe para
corrigir o que estiver errado enquanto quebrar ainda é barato, não para congelar.

Por que agora: os planos 1 a 6 fecharam a matemática, a validação e a performance; falta o caminho do código até o
usuário. O principal consumidor é Python (decisão do autor); C++ é secundário, mas entra.

## 2. Objetivos, não-objetivos, anti-objetivos

### Objetivos

1. **Análise de API** da interface C, do pacote Python e do wrapper C++ contra o spec §5–7 e as regras de
   estabilidade §5.1: nomes, códigos de erro e mensagens, posse dos handles, simetria entre setters, o que sobra
   (anti-objetivo "API inchada") e o que falta para os usos do §6. Saída: um documento de achados com decisão por
   item; as correções aceitas entram neste plano.
2. **Wrapper C++23 header-only** `cpp/include/portopt/portopt.hpp` (spec §7): `std::expected`, `std::span`, RAII nos
   handles, sobre a interface C, com testes e um consumidor de exemplo.
3. **Wheels da matriz completa** (spec §11): Windows x64, macOS arm64 e x64, Linux glibc (manylinux_2_28) e musl
   (musllinux_1_2) em x64 e arm64 — sete wheels `py3-none-<plataforma>`, cada uma com a tag conferida e testada
   instalada num ambiente limpo.
4. **Release** por tag `v*`: wheels + bibliotecas nativas (com a import library do Windows) + dump de ABI + SBOM numa
   GitHub Release, e as wheels no PyPI por trusted publishing. Pré-release quando a tag tem `-` (`v0.1.0-rc.1`).
5. **CI completa**: o `ci.yml` atual, mais `wheels.yml` (roda em todo push que toca o build, para que uma release
   nunca compile uma plataforma pela primeira vez), `cpp.yml` (wrapper) e a compilação dos benchmarks
   (`PO_BUILD_BENCHMARKS=ON` + seus dois self-tests), apontada pela review final do plano 6.
6. **Empacotamento**: versão única nos três lugares (CMake `PO_VERSION`, `pyproject.toml`, CHANGELOG) com um script
   que confere; regras de `install` no CMake (lib, header C, wrapper C++, `find_package(PortOpt)`); README com
   instalação e exemplos Python e C++; CHANGELOG com a seção da release.

### Não-objetivos

- Congelar a ABI ou subir `PO_ABI_VERSION` para 1 (isso é a 1.0, plano futuro).
- Lote com threads (`po_solve_batch`, spec §5.4), objetivos ou restrições novas.
- Alavancas de performance da fase 2 (profile-2.md): despacho por ISA, Dekker com FMA, classe B, warm start.
- conda-forge, Homebrew, vcpkg, Conan.
- Wheels para Windows arm64, Linux 32 bits, PyPy ou versões específicas de CPython (o pacote é ctypes, um wheel por
  plataforma serve a qualquer Python ≥ 3.11).

### Anti-objetivos

- **Publicar sem o autor.** Criar o repositório remoto, o projeto no PyPI, configurar o trusted publisher e empurrar a
  tag são ações do autor (externas e irreversíveis). O plano entrega tudo pronto e ensaiado; não publica.
- **Uma plataforma compilada pela primeira vez na release.** Toda plataforma passa pelo `wheels.yml` antes.
- **Wheel que instala e não carrega.** Cada wheel é testada instalada, fora do checkout, com a suíte pytest.
- **Mudança de resultado numérico.** Nada neste plano muda um bit de solução; correções de API mexem em nomes,
  validação e erros, não na matemática. Se uma correção de API tocar o núcleo, ela é bit a bit idêntica (método do
  plano 6) ou sai do plano.
- **Correção de API sem migração escrita.** Toda quebra entra no CHANGELOG com o "antes → depois".

## 3. Restrições

- As dos planos 1–6: estilo do `STYLEGUIDE.md`, commits sem `Co-Authored-By`, `-ffp-contract=off` em todo alvo que
  compila código de produção, sem `long double` em `src/`, documentos de `docs/superpowers/` em português.
- O build MinGW das bibliotecas distribuídas liga `libstdc++`, `libgcc` e `winpthread` estaticamente (spec §11); os
  binários distribuídos são genéricos (sem `-march=native`), como o `gcc-bench`.
- `hatch_build.py` já deriva a tag da glibc e do macOS mínimo da própria biblioteca; o pipeline o usa como está.
- macOS: `MACOSX_DEPLOYMENT_TARGET=11.0` para arm64 (o mínimo possível); x64 o que o runner e a libc++ permitirem com
  `std::expected` (Xcode 15+), medido e registrado.
- Linux: GCC 15 dentro dos containers manylinux_2_28 e musllinux_1_2 (o projeto exige GCC 15+).
- Ações do GitHub com versão fixada (major), como no VolEngine.

## 4. Abordagens consideradas

| Abordagem | Troca | Decisão |
|---|---|---|
| **A. Portar o pipeline do VolEngine** (`wheels.yml` com containers manylinux/musllinux, `build_in_container.sh`, `build_wheel.sh`, `check_wheel_tag.py`, `check_native_deps.py`, `release_version.py`; `release.yml` com import library, dump de ABI, SBOM, GitHub Release e PyPI por trusted publishing; `cpp.yml`) | Herda as escolhas do VolEngine, que são as do spec §11; já publicou a v1.0.0 | **Escolhida pelo autor** |
| B. cibuildwheel | Padrão do ecossistema, mas modela extensões CPython (um wheel por versão de Python); o nosso é `py3-none` via ctypes, e diverge do VolEngine | Descartada |
| C. scikit-build-core (compilar a lib dentro do `pip wheel`) | Troca o backend hatchling e o fluxo; mudança grande sem ganho | Descartada |

## 5. Decisões já tomadas

- Release **0.x** (autor). Primeira tag proposta: `v0.1.0-rc.1` (ensaio na GitHub Release como pré-release e no
  TestPyPI), depois `v0.1.0`.
- Matriz **completa** de sete wheels (autor).
- Wrapper C++ **entra** neste plano (autor: secundário, mas interessante ter).
- Pipeline **A**, portado do VolEngine (autor).
- **Repositório:** `https://github.com/GabrielMarquezMatte/PortfolioOptimization`, público, branch padrão `master`,
  hoje com a implementação em R (último commit `2b16551`, 2022; 1 estrela, nenhum fork, nenhuma issue ou PR aberta).
  O autor quer sobrepor a implementação, sem apagar o repo: o R fica no branch `r-legacy` e na tag `r-final`, e o
  `master` local (C++) substitui o `master` remoto por force-push, com o ok explícito do autor no momento, depois de
  conferir que `r-legacy` e `r-final` existem no remoto. Público: runners arm64 nativos gratuitos.
- **Nome no PyPI:** `portopt` está ocupado por outro projeto, que instala os módulos `optimizer` e `cli` (não
  `portopt`). Só o nome de distribuição muda; o import `portopt`, `libportopt.*` e `find_package(PortOpt)` ficam.
  Nome de distribuição: **`pyportopt`** (autor; livre no PyPI em 2026-11-17): `pip install pyportopt`, `import portopt`.
- **Ambientes de publicação** (decisão do Claude): `release.yml` declara os environments `pypi` e `testpypi`; o autor
  registra o trusted publisher com esses nomes.
- **Versão de rc** (decisão do Claude, como no VolEngine): `release_version.py` deriva da tag, no build, a versão
  PEP 440 do wheel e a `x.y.z` do CMake (`v0.1.0-rc.1` → `0.1.0rc1` e `0.1.0`); nenhum commit de versão de rc. O rc vai
  ao TestPyPI; a tag final, ao PyPI.
- **Contrato de ABI** (decisão do Claude): portar também `scripts/abi_contract.py`, `abi_dump.sh` e
  `install_abi_tools.sh` do VolEngine; a lista de exports sai do header, alimenta a import library do Windows, e o CI
  confere que os exports da biblioteca batem com ela.
- **Análise de API num plano só** (decisão do Claude): a Task 1 produz o documento e para; as tasks seguintes
  implementam a lista de correções que o autor aprovar. Renomear o projeto está fora da análise.

## 6. Desenho

### 6.1 Análise de API (primeiro, porque o resto empacota o resultado)

Documento `docs/guide/api-review.md` (ou seção do ARCHITECTURE.md; escolha do builder): para cada função de
`portopt.h`, cada nome público de `portopt` (Python) e cada tipo do wrapper, uma linha com o achado e a decisão
(manter / mudar / remover / acrescentar). Critérios: spec §5.1 (setters, enums, erros, handles independentes),
anti-objetivo "API inchada", consistência de nomes entre C, Python e C++, mensagens de erro, documentação no header.
As mudanças aceitas são implementadas depois do documento e antes do wrapper C++, uma por commit, com o CHANGELOG
"antes → depois". O autor revê o documento antes de qualquer quebra (ponto de parada do plano).

### 6.2 Wrapper C++

`cpp/include/portopt/portopt.hpp` (um arquivo, ou poucos como o `volengine.hpp`; escolha do builder conforme o
tamanho): handles RAII movíveis e não copiáveis, `std::expected<T, Error>` com status e mensagem de
`po_last_error()`, `std::span` para vetores e matrizes. Testes em `cpp/tests/` (GoogleTest, contra a biblioteca
construída), incluindo um consumidor de exemplo via `find_package(PortOpt)` instalado e um via `add_subdirectory`.
`cpp.yml` roda em GCC, Clang (macOS) e MSVC, e uma vez sob ASan/UBSan.

### 6.3 Wheels

`python/scripts/` com os scripts portados do VolEngine. `wheels.yml` com sete entradas de matriz: Linux em containers
`quay.io/pypa/manylinux_2_28_{x86_64,aarch64}` e `musllinux_1_2_{x86_64,aarch64}` (arm64 em runner arm64 nativo do
GitHub), macOS arm64 e x64 em runners nativos, Windows x64 com MSYS2 UCRT64. Cada job: biblioteca genérica → wheel →
`check_wheel_tag.py` (tag esperada fixada na matriz) → `check_native_deps.py` (sem dependências além das do sistema)
→ instala num venv limpo e roda pytest fora do checkout. Artefatos: a wheel e a biblioteca nativa nomeada por
plataforma.

### 6.4 Release

`release.yml`: gatilho em tag `v*` e `workflow_dispatch` com `dry_run`. Confere a versão (`release_version.py`: tag
= CMake = pyproject = CHANGELOG), chama `wheels.yml`, gera a import library do Windows e o dump de ABI (libabigail) da
biblioteca Linux x64, publica a GitHub Release com bibliotecas, import library, dump e SBOM (pré-release se a tag tem
`-`), e então o PyPI (ou TestPyPI para `-rc`) por trusted publishing. Em `dry_run` tudo roda menos os dois passos de
publicação.

### 6.5 Empacotamento e docs

`install(TARGETS ...)`, `install(FILES include/portopt/portopt.h ...)`, wrapper C++ e `PortOptConfig.cmake` +
`PortOptConfigVersion.cmake`. README: instalação (`pip install pyportopt`; C/C++ pela GitHub Release ou pelo
`find_package`), exemplo Python e exemplo C++. CHANGELOG: seção `[0.1.0]` com data.

### 6.6 Ordem

1. Análise de API (documento) → parada para o autor revisar.
2. Correções de API aceitas.
3. Wrapper C++ e `cpp.yml`.
4. `install` no CMake e conferência de versão.
4b. Troca do repositório remoto (`r-legacy` + `r-final`, depois force-push do `master` com o ok do autor), para que
    os workflows rodem nos runners do GitHub.
5. `wheels.yml` com os scripts, verde nas sete plataformas (verificado por push no remoto).
6. Job dos benchmarks no CI.
7. `release.yml` ensaiado em `dry_run`.
8. README, CHANGELOG, preparação da tag (o autor empurra).

## 7. Deixado ao builder

- Formato e local exato do documento de análise de API.
- Divisão do wrapper C++ em arquivos.
- O macOS mínimo do x64 (medido).
- Detalhes dos scripts portados (nomes de assets: `portopt-native-<plataforma>.{dll,so,dylib}`, como o VolEngine).
- Se o ensaio da publicação usa TestPyPI (recomendado) ou só `dry_run`.

## 8. O que depende do autor (fora do plano)

- Autorizar, no momento, o force-push do `master` sobre a implementação em R.
- Registrar no PyPI e no TestPyPI o trusted publisher (pending publisher) do nome de distribuição escolhido,
  apontando para `release.yml` e os environments `pypi` / `testpypi`.
- Revisar o documento de análise de API antes das quebras.
- Empurrar a tag de release.
