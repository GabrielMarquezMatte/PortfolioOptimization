# Pesquisa: método de performance para a v1 (porte do VolEngine) e alavancas de otimização do PortOpt

> **Sobre este arquivo.** Insumo para um plano de performance, não um plano nem uma spec: nada aqui altera código,
> spec ou planos até o autor decidir. Data: 2026-11-16. Nenhum código do repositório foi editado e nada foi commitado.
> Contexto lido: `docs/performance/` e `tests/benchmarks/` do VolEngine, `CMakePresets.json`, `CMakeLists.txt`,
> `ARCHITECTURE.md` e `STYLEGUIDE.md` dele; no PortOpt, `CMakeLists.txt`, `CMakePresets.json`, `src/` inteiro das áreas
> citadas, `docs/guide/validation.md`, `docs/superpowers/specs/2026-11-15-portopt-design.md` (§4.1, §5.4, §8, §10,
> §11) e `docs/research/2026-11-16-risk-parity-and-numerics.md` (formato).
>
> **Convenção de verificação.** "Verificado" significa lido na fonte primária (artigo, documentação oficial do
> compilador ou da ferramenta, código ou README do projeto, arquivo do repositório), com a página ou seção quando a fonte
> a tem. "Não verificado" significa que a fonte não pôde ser lida (paywall, PDF sem texto, 403) e a afirmação vem de
> metadado bibliográfico, citação secundária ou memória, marcada como tal. "Derivação própria" é um argumento escrito
> aqui, conferível linha a linha, sem fonte que o enuncie para este código. **"Experimento exploratório"** é uma medição
> feita nesta sessão em arquivos soltos fora do repositório (scratchpad), em um só computador, sem isolamento, com
> `g++ 16.2.0` (MSYS2 UCRT64), `-O3 -ffp-contract=off`, alvo padrão do GCC (SSE2, sem AVX2 nem FMA), uma execução
> cada: serve para ordenar hipóteses, **não é linha de base** e deve ser refeito no harness do §3 (o CPU não foi
> identificado nesta sessão; o VolEngine registra um Ryzen 7 5700X na mesma máquina de desenvolvimento, **não
> verificado** que seja o mesmo aqui).

## Resumo

O VolEngine deixa um método completo e portável: Google Benchmark via FetchContent, uma suíte rápida de 3 repetições
com filtro fixo, A/B por mediana contra desvio padrão, executável de perfil, uProf sobre build MSVC com PDB (o binário
GCC não é simbolizável no Windows), números oficiais do build genérico. Tudo isso porta quase sem mudança; o que muda é
o conteúdo (instâncias, contadores, comparação de soluções em vez de log-verossimilhança).

Três achados de leitura do código, confirmados por experimentos exploratórios, ordenam o plano:

1. **O custo de inicializar a base do Goldfarb-Idnani é desproporcional.** `empty_active_set` (`active_set.hpp:57-73`)
   resolve n sistemas triangulares com acesso de passo n (`cholesky.cpp:72`): ≈ n³/2 multiplicações-somas com falta de
   cache. A versão orientada a linhas, que preserva a ordem de operações de cada elemento, dá os **mesmos bits** e foi
   8,4 vezes mais rápida em n = 1000 (350 ms contra 42 ms, experimento exploratório), maior que o próprio Cholesky
   (105 ms). Cada `qp::solve` paga isso, inclusive cada subproblema da sequência de QPs do RISK_PARITY e cada ponto do
   caminho do MAX_RETURN.
2. **Em build genérico, `std::fma` vira chamada de `libm`.** O `compensated_dot` (Dot2) usa `std::fma`
   (`objective.cpp:42`); sem `-mfma` o GCC emite `call fma` (verificado por compilação nesta sessão). A mesma
   TwoProduct feita por Dekker/Veltkamp deu **os mesmos bits** em 2000 casos aleatórios e foi 1,9 vez mais rápida
   (4,13 para 2,20 ns por elemento); com FMA de hardware inline, 1,11 ns.
3. **Muitos kernels têm reescritas que não alteram nenhum bit** (loop interchange e blocagem que mantêm, para cada
   elemento, a sequência de operações): Cholesky (2,0 vezes só com laço right-looking), reflexão de Householder
   (4,0 vezes), `gauss_newton` (metade do trabalho por simetria). Isso separa as otimizações em duas classes, de
   risco numérico zero e de risco real, e o plano deve esgotar a primeira antes de tocar a segunda.

O que **não** vale a pena: otimizar a família implícita de 2^n linhas (custo O(n) por candidato, não é hot spot e mexer
nela altera a ordem de geração que o caminho paramétrico compara), BLAS opcional (contra a spec §4.1/§10 e contra a
reprodutibilidade), `-ffast-math`, paralelismo dentro de uma solução.

---

## 1. Tabela resumo das alavancas

Ganho: estimativa **por leitura** mais, quando existe, o experimento exploratório (marcado EE); nenhuma está medida na
biblioteca. Risco numérico: **0** = bit a bit idêntico, verificável por `memcmp`; **B** = muda arredondamento, exige a
banda; **A** = muda o caminho do algoritmo (iterações, conjunto ativo).

| ID | Alavanca | Ganho esperado | Risco | Custo | Veredito |
|---|---|---|---|---|---|
| L1 | `empty_active_set`: base J = L⁻ᵀ por linhas, pulando zeros, mesma ordem por elemento | EE: 8,4x no kernel em n = 1000; remove o termo n³/2 de **todo** `qp::solve` | 0 | baixo | **fazer primeiro** |
| L2 | Cholesky right-looking / com registradores, ordem por elemento preservada | EE: 2,0x (105 para 53 ms, n = 1000) só com right-looking ingênuo; mais com tiling | 0 | médio | fazer |
| L3 | Reflexão de Householder por linhas (`reflect` em `qr.cpp:53`) | EE: 4,0x na fase de atualização (1500x1000) | 0 | baixo | fazer |
| L4 | `gauss_newton`: metade por simetria + blocagem (i interno, ordem por elemento) | até 2x em compute, e de memória-limitado a compute-limitado em n ≥ 300; domina RISK_PARITY | 0 (derivação própria, §4.2) | médio | fazer |
| L5 | Fast path de linha unitária (cotas) em slack e projeção do GI | tira ≈ 2n² a 3n² flops e 16 MB de leitura por candidato em n = 1000 com cotas | 0 salvo sinal de zero | médio | fazer após L1 |
| L6 | Dot2 sem chamada `libm`: TwoProduct de Dekker; depois, despacho manual `target("avx2,fma")` | EE: 4,13 → 2,20 → 1,11 ns/elem; só em `covariance_times` (O(n²) por produto) | 0 (Dekker, ver hazard 3); B se usar FMA fora de Dot2 | médio | fazer Dekker; despacho depois |
| L7 | Cache de J₀ = L⁻ᵀ entre solves com o mesmo fator (caminho do MAX_RETURN) | elimina o termo n³/2 de cada ponto do caminho (dezenas de solves) | 0 | médio | fazer após L1 |
| L8 | Higiene de alocação no laço do GI e nos problemas | pequeno em n grande; relevante em n = 10 (muitas iterações pequenas); confirmar com contador | 0 | baixo | fazer se o contador mostrar |
| L9 | Cache de `Moments` nas meias-buscas de Armijo (`value_change`, `objective.cpp:209`) | um Dot2 O(n²) por halving; pequeno | 0 | baixo | opcional |
| L10 | `-flto` no core | pequeno; chamadas entre TUs (`dot`, `solve_lower`) em n = 10 | 0 se `-ffp-contract=off` valer no link (não verificado) | baixo | medir depois |
| L11 | Kernels com reassociação (multi-acumulador, GEMM com acumulador zerado, SIMD de reduções) | 3x a 4x em Cholesky/produtos | **B**, banda fina em um caso (§5) | alto | só após L1-L7, com A/B de pesos |
| L12 | Multiversão por ISA (AVX2) só para kernels elementwise/axpy | 2x a 4x sobre SSE2 nesses kernels | 0 (elementwise); B em reduções | alto | segunda fase |
| L13 | Warm start do QP entre iterações proximais do RISK_PARITY | potencialmente o maior ganho algorítmico | **A** | alto | estudo separado, não otimização |
| L14 | Paralelismo de lote (`po_solve_batch`) | ≈ linear em núcleos, resultados idênticos por thread | 0 | médio | plano próprio (spec §5.4) |
| L15 | BLAS opcional | nenhum em n ≤ 1000 frente aos kernels próprios bem escritos (não medido) | B e não reprodutível | alto | **não** |
| L16 | Paralelismo dentro de uma solução (álgebra linear com threads, n ≥ 500) | médio em n ≥ 500 | **B** (ordem de redução) | alto | **não** na v1 |
| L17 | Armazenamento empacotado de matriz simétrica | memória ÷ 2; perde vetorização | 0 | médio | não |
| L18 | Reduzir esforço do RISK_PARITY (estágios, iterações) | grande | **A** | baixo | não é performance: é mudança de algoritmo |

---

## 2. O método do VolEngine, portado

### 2.1 O que o VolEngine tem (lido)

**Ferramentas** (`docs/performance/measuring.md`, `tests/benchmarks/`, `CMakePresets.json`, `CMakeLists.txt:269-294`):

- **Google Benchmark v1.9.1** por `FetchContent_Declare(benchmark URL .../v1.9.1.tar.gz URL_HASH SHA256=32131c08ee31eeff2c8968d7e874f3cb648034377dfc32a4c377fa8796d84981 SYSTEM)`, com
  `BENCHMARK_ENABLE_TESTING/GTEST_TESTS/INSTALL/WERROR` desligados, sob `if(VE_BUILD_BENCHMARKS)`. Executáveis:
  `volengine_benchmarks` (a suíte) e `volengine_grid_profile` (perfil de grade). Os dois recebem `-ffp-contract=off`.
- **Suíte rápida** (`quick.py`): filtro fixo `QUICK_FILTER` (regex de `BM_Loglik`, `BM_FitOne`, `BM_TccGridSample`),
  `--benchmark_repetitions=3 --benchmark_report_aggregates_only=true`, saída JSON; cerca de 35 s.
- **`quick.py compare A.json B.json`**: lê só agregados, imprime mediana A, mediana B, razão B/A e a marca `noise`
  quando |B − A| ≤ max(desvio A, desvio B). Biblioteca padrão só; tem `selftest`.
- **`compare-profile`**: compara dois JSON do executável de perfil: razão de tempo, razão de avaliações, **pior
  variação relativa da log-verossimilhança por ajuste** (e qual) e contagem de status. Regra: uma mudança de otimizador
  não pode baixar nenhuma log-verossimilhança além de 1e-9 relativo.
- **Executável de perfil** (`grid_profile.cpp`): roda a grade inteira (1080 ajustes) uma vez pelo caminho de lote;
  imprime `wall_s`, soma dos tempos, `pool_use`, tabelas por spec/distribuição/ordem, status; `--json`.
- **uProf** (`AMDuProfCLI collect --config tbp`, `report` com `--bin-path` e `--symbol-path`) sobre o preset
  `msvc-profile` (`/O2 /Ob2 /DNDEBUG /Zi`, linker `/DEBUG /OPT:REF /OPT:ICF`, benchmarks ON, testes OFF). Motivo,
  verificado no arquivo: o `report` do uProf quebra em sessão de binário DWARF e `-gcodeview` falha no assembler do
  MinGW. O build MSVC é mais lento que o GCC (20 s contra 7,7 s na mesma fatia), então as **participações são
  indicativas** e o ganho de uma mudança sempre vem do `gcc-bench`.
- **Assembly** (`ARCHITECTURE.md`, "Performance method"): `objdump -d -C -M intel --no-show-raw-insn` (ou `-S` com
  `-g`) de um objeto `gcc-bench`; `-fopt-info-vec-missed`; procurar chamadas que deveriam estar inlined, alocações em
  laço, cópias, divisões e transcendentais içáveis. Os achados vão para `docs/performance/`, nunca para comentário.
- **Preset `gcc-bench`**: `gcc-base` + `Release` + `VE_BUILD_BENCHMARKS=ON` + `VE_BUILD_TESTS=OFF`. **Sem
  `-march=native`**: o preset `native-release` existe à parte e não é fonte de números oficiais.

**Regras** (`measuring.md`, `benchmarks.md`, `STYLEGUIDE.md` "Performance Notes"):

- Medir antes e depois **na mesma sessão, um logo após o outro**; nunca comparar com tabela antiga (execuções de dias
  diferentes moveram medianas além do desvio). Diferença menor que o desvio não é resultado.
- Números oficiais vêm do build genérico porque o binário distribuído é genérico.
- Registrar a máquina em todo número; nunca derivar razão entre máquinas; benchmark contra outra biblioteca compara
  trabalho equivalente ou diz que não compara (`python-vs-arch.md`).
- Falhar alto: o benchmark faz `SkipWithError` se o ajuste não termina em `CONVERGED`/`MAX_ITER`; "o harness perigoso é
  o que mede zero trabalho e publica um número".
- Resultados bit a bit iguais para qualquer número de threads; nada de redução cuja ordem dependa do agendamento.
- Alocar uma vez fora do laço; cada otimização aceita com número ao lado; o que ela achou vai para `docs/`.
- Contadores e laço: `iterations`, `evaluations`, `fits_per_s` (`kIsRate`), status
  (`converged`, `max_iter`, `other`), e `DoNotOptimize(result.loglik)` no laço.

**O que os perfis 5a, 5b e 5c registraram e fizeram, com o que rejeitaram:**

| Plano | Perfil (uProf, MSVC) | Feito (A/B no `gcc-bench`) | Rejeitado, com medida |
|---|---|---|---|
| 5a | 27% em `ad::chain<16>`, resto em `Dual` | `Dual<12>` para k = 9..12 (0,87); memória L-BFGS 10→40 (0,50); um ajuste por bloco, o mais caro primeiro (0,98 de uso do pool); opção sem erros-padrão (0,89); gradiente adjunto escrito à mão (0,71) | `f` em `double` antes do gradiente (perda líquida estimada de 13%); interpolação no zoom (2,2% menos avaliações, tempo dentro do ruído); alocações (nenhuma no perfil); L-BFGS `ftol` 1e-15 não mudou tempo |
| 5b | `factor_entry<Dual<8>>` 65-72%, Cholesky + solve 79-82% | kernels SIMD (Highway) com Cholesky em lanes de matrizes de t diferentes: `double` 0,42x; gradiente por tangentes: ajuste DCC 0,15x a 0,25x; laço vetorial das densidades simétricas 0,40x a 0,76x | Cholesky vetorial **por elemento** (16 ms contra 8 ms do escalar: cauda mascarada e soma horizontal custam mais que os produtos em linhas de ≈ 7 entradas); "store" mascarado é lento em AMD |
| 5c | quantis Student: 80% do tempo de `static_mvt_ml_param` | Kendall por Knight O(T log T) (0,04x); warm start do perfil de forma (36 s para 3,6 s); gradiente por tangentes nas cópulas; Boost em `double` (`promote_double<false>`, 2,3x a 5x); tabela de quantis (40x); inclinação de K por trapézio | Gradiente por tangentes para erros-padrão (3 avaliações por coluna, ganho nenhum); polimento por segmento de Gauss-Legendre (não cai na raiz de `quantile`) |

Lições transferíveis: **o perfil decide a ordem** (5b: a hipótese "30-60% por kernel" estava errada, o Cholesky era 80%);
**o experimento barato antes do kernel** (spike de 5a: tape genérica 5-6x contra adjunto 1,3-1,5x); **cada rejeição é
registrada com medida**; **o ruído da máquina é admitido** (bimodalidade de `aparch11_ged/dual`).

### 2.2 O que porta como está

| Item | Estado |
|---|---|
| `quick.py` (`run`, `compare`, `selftest`) | Porta quase literal: trocar `DEFAULT_EXE`, `QUICK_FILTER` e o nome do JSON. |
| `FetchContent` do Google Benchmark (mesma URL e hash) | Literal; a opção `PO_BUILD_BENCHMARKS` já existe em `CMakeLists.txt:24` e não há bloco que a use. |
| Preset `gcc-bench` | Literal com `PO_BUILD_BENCHMARKS`/`PO_BUILD_TESTS`. |
| Preset `msvc-profile` | Literal (flags `/O2 /Ob2 /DNDEBUG /Zi`, `/DEBUG /OPT:REF /OPT:ICF`). |
| Receita uProf e a ressalva "só localiza, não mede" | Literal. |
| Regras de A/B, nunca contra tabela antiga, máquina em cada número | Literal. |
| Receita de assembly | Literal; acrescentar a checagem de `vfmadd` do §3.6. |
| `MinGW`: `-static-libgcc -static-libstdc++` nos executáveis | Literal (`CMakeLists.txt:284-286` do VolEngine; no PortOpt já está em `portopt_tests`). |

### 2.3 O que muda

| Item | Por que muda |
|---|---|
| `compare-profile` (log-verossimilhança) | Vira **`compare-solutions`**: compara dumps de soluções (§3.5) com contagem de resultados **bit a bit iguais**, pior erro relativo de peso, pior diferença de objetivo, mudança de status, razão de iterações. A regra "1e-9 relativo" do VolEngine não serve: aqui o critério é a banda KKT e, para alavancas de risco 0, igualdade exata. |
| Executável de perfil de grade | Vira `portopt_profile`: percorre as instâncias do §3 (objetivo × N × variante), imprime tempo, iterações do QP, contagem de QPs e de alocações, e o dump de soluções. Não há "lote" para estressar. |
| Contadores | Trocam `iterations`/`evaluations` por `qp_iterations` (já em `Result::iterations`), `qp_solves`, `allocations`, `descent_iterations`, `newton_polish_steps`. |
| Instâncias | Sem séries simuladas; geradores de `tests/unit/random_covariance.hpp` (`random_covariance(n, seed, condition)`, `random_vector`), que já são determinísticos (`mt19937_64`, Box-Muller próprio). |
| Núcleo de comparação | O VolEngine compara com `arch` e `rugarch`; aqui a referência é a **própria versão anterior** (A/B) e os testes KKT. |
| Python | O VolEngine tem `python-vs-arch.md`; aqui o equivalente é a sobrecarga do ctypes contra o tempo de solve (§4.2-g), não uma comparação entre bibliotecas. |
| `-ffp-contract=off` em benchmarks | O VolEngine aplica `VE_FLOAT_FLAGS` ao benchmark; aqui o benchmark que **chama** o core não precisa, mas o executável de dump e qualquer cópia de kernel de referência precisam do mesmo flag (hazard 2). |
| Estilo | `tests/style/check_style.py` varre `tests/` (`SCANNED_DIRECTORIES`): o código de benchmark segue funções ≤ 60 linhas, sem `else`, etc. |

---

## 3. A suíte de benchmarks a escrever

### 3.1 Estrutura

```
tests/benchmarks/
  quick.py                  run / compare / compare-solutions / selftest
  instances.hpp             geradores com cache (covariância n = 1000 leva segundos para gerar)
  solve_benchmarks.cpp      BM_Solve/*  (objetivos, ponta a ponta por C++ po::problems::solve*)
  kernel_benchmarks.cpp     BM_Kernel/* (micro, sobre os headers internos)
  profile.cpp               portopt_profile (grade + dump de soluções + contadores)
  alloc_counter.hpp         operator new/delete substituídos, só no executável de benchmark
```

`portopt_core` é `STATIC` com `PUBLIC` `src` e `include` (`CMakeLists.txt:61`), então o benchmark inclui
`linalg/cholesky.hpp`, `qp/goldfarb_idnani.hpp`, `riskparity/objective.hpp` diretamente, sem passar pelo C API; a
camada C é medida à parte (§3.4).

### 3.2 Casos `BM_Solve` (ponta a ponta)

Nome: `BM_Solve/<objetivo>/<variante>/n:<N>`. Instância: `random_covariance(N, seed = 1, condition = 1e2)` como padrão
benigno, `random_vector` para μ em [0,02, 0,12], `risk_free = 0,01` só nas variantes `cash`. Teto de caixa: `max(0,3, 2/N)`.
Todo caso afirma `status == Optimal` e `weights` finitos; senão `SkipWithError` (regra de falhar alto).

| Objetivo | Variantes | N = 10 | N = 100 | N = 1000 |
|---|---|:-:|:-:|:-:|
| MIN_VARIANCE | `free`, `long_only`, `box`, `long_only+target_return`, `leverage1.5`, `cash` | sim | sim | sim |
| UTILITY, IMPORTANCE, KELLY | `free`, `long_only`, `box`, `leverage1.5` | sim | sim | UTILITY só `long_only` |
| MAX_RETURN | `long_only+vol_cap`, `free+vol_cap`, `leverage1.5+vol_cap`, `cash+vol_cap` (cap = 1,5 σ_min) | sim | sim | `long_only` (estresse) |
| MAX_SHARPE, MAX_DIVERSIFICATION | `long_only`, `box`, `leverage1.5`, `free` | sim | sim | `long_only` |
| TRACKING_ERROR | `long_only`, `box` (T = 2N, T = 4N) | sim | sim | `long_only`, T = 2N |
| RISK_PARITY | ver abaixo | sim | parte | só `spinu` |

**RISK_PARITY**, o único objetivo com tempo de segundos:

| Variante | N = 10 | N = 30 | N = 50 | N = 100 | N = 1000 |
|---|:-:|:-:|:-:|:-:|:-:|
| `spinu` (sem caixa dentro, sem alvo: ponto de Spinu devolvido) | sim | sim | sim | sim | sim |
| `box_active` (caixa que corta o ponto de Spinu: descida + polimento) | sim | sim | sim | sim | estresse (`Iterations(1)`) |
| `target_return` | sim | sim | sim | sim | — |
| `target_volatility` (cond 1e2) | sim | sim | sim | sim | — |
| `target_volatility`, **cond 1e8, sem caixa** (o caso de 9,2 s em n = 50, 19,6 s em n = 100, `validation.md:1071`) | — | — | `Iterations(1)`, fora da suíte rápida | idem | — |
| `uneven_budget` (orçamento desigual) | sim | — | sim | — | — |

Cada caso reporta `state.counters`: `qp_iterations` (`Result::iterations`, que no caminho paramétrico acumula as
iterações de todos os QPs, `parametric_path.cpp:282`), `qp_solves`, `allocations`, `descent_iterations`,
`polish_steps`, `status` e, para RISK_PARITY, `residual`. `qp_solves` e as demais contagens internas exigem ou uma
definição de compilação `PO_BENCH_INSTRUMENT` num build do core só para benchmark (contador `thread_local` em
`qp::solve`) ou a recontagem por `Result::iterations` e pelo contador de alocação; o plano decide (custo: pequeno, mas
toca `src/`).

Para os casos de segundos, `->Iterations(1)` e `->Repetitions(3)` (como `BM_TccGrid`); para os de microssegundos,
`->MinTime(0.5)`; `UseRealTime()` só onde houver thread. `DoNotOptimize(result.weights.data())` e
`ClobberMemory()` no laço (ambos descritos no User Guide, §"Preventing Optimization", verificado).

### 3.3 Micro-benchmarks `BM_Kernel`

Tamanhos 10, 100, 300 e 1000 (`->Args`), com contador de `flops` (`kIsRate`) para ver a fração do pico:

| Nome | O que mede | Fonte |
|---|---|---|
| `BM_Kernel/cholesky/n` | `linalg::cholesky` em Σ SPD | `cholesky.cpp:12-50` |
| `BM_Kernel/solve_lower/n`, `solve_lower_transposed/n` | um vetor | `cholesky.cpp:52-76` |
| `BM_Kernel/basis_init/n` | `empty_active_set` (os n sistemas) | `active_set.hpp:57-73` |
| `BM_Kernel/qr_square/n`, `qr_tall/T,n` | `householder_qr` + `least_squares` | `qr.cpp:82-135` |
| `BM_Kernel/qp_iteration/n/bounds` | uma chamada de `compute_directions` + `append_active` (por `qp::solve` de um QP long-only, tempo por iteração = tempo ÷ `Solution::iterations`) | `active_set.hpp:157-190` |
| `BM_Kernel/qp_scan/n/bounds` | `most_violated_explicit` em 2n linhas unitárias densas | `goldfarb_idnani.cpp:191-206` |
| `BM_Kernel/compensated_dot/n` e `plain_dot/n` | Dot2 e a soma simples, mesmo vetor | `objective.cpp:35-50` |
| `BM_Kernel/sigma_w/n` | `covariance_times` (n Dot2) | `objective.cpp:54-63` |
| `BM_Kernel/gauss_newton/n` | J'J | `objective.cpp:187-202` |
| `BM_Kernel/linearize/n` | `linearize(..., false)` e `(..., true)` | `objective.cpp:170-185` |
| `BM_Kernel/spinu_step/n` | `newton_step` (Cholesky + dois solves) | `spinu.cpp:114-133` |

### 3.4 Camada C e Python (fora da suíte rápida)

`BM_CApi/solve/<objetivo>/n` chama `po_solve` com handles já criados, e um caso com `po_data_create` incluso, para
ver o custo do Cholesky de `Data::create` (`data.cpp:193`) isolado. O equivalente do `python-vs-arch.md` do VolEngine
(§4.2-g): script `python/benchmarks/bench_overhead.py` com tempo de uma chamada trivial pelo ctypes contra o tempo de
`solve` em n = 10, 100.

### 3.5 Filtro da suíte rápida, comparação de soluções e presets

**`QUICK_FILTER`** (alvo: ≈ 30-40 s, como o do VolEngine; 3 repetições, agregados somente):

```
^BM_Solve/(min_variance/(long_only|leverage1.5)|utility/box|max_return/long_only+vol_cap|max_sharpe/long_only|
 max_diversification/box|tracking_error/long_only|risk_parity/(spinu|box_active|target_volatility))/n:(10|100)$
|^BM_Solve/min_variance/long_only/n:1000$
|^BM_Kernel/(cholesky|basis_init|qr_square|sigma_w|gauss_newton|compensated_dot)/(100|1000)$
```

(`risk_parity/target_volatility` e `box_active` só em N = 10 e 30 na rápida, para caber no tempo.)

**`portopt_profile --json solutions.json`**: para cada instância do conjunto, grava `status`, `iterations`, `x` e
multiplicadores como **hex-float** (`%a`, exato), o conjunto ativo (`active_inequalities`, sinais de
`gross_exposure_signs`) e o objetivo calculado por `metrics::` com `compensated_dot`. **`quick.py compare-solutions A B`**
imprime: número de instâncias bit a bit iguais; pior erro relativo de peso (`max|Δw| / max|w|`) e qual instância; pior
|Δf|; instâncias cujo status ou conjunto ativo mudou; razão de iterações. Para alavanca de risco 0 o critério de
aceitação é **todas iguais bit a bit**; para B, pesos e objetivos dentro da banda da instância (`expect_user_kkt`,
`expect_risk_parity_kkt`) e status, conjunto ativo e iterações idênticos; para A, estudo próprio.

**Presets a acrescentar** a `CMakePresets.json` (hoje: `gcc-debug`, `gcc-release`, `native-release`, `msvc-release`):

- `gcc-bench`: `inherits gcc-base`, `Release`, `PO_BUILD_BENCHMARKS=ON`, `PO_BUILD_TESTS=OFF`, **sem** `-march`.
- `msvc-profile`: `inherits msvc-release`, mesmos `CMAKE_CXX_FLAGS_RELEASE`/`CMAKE_EXE_LINKER_FLAGS_RELEASE` do VolEngine.
- Build presets `gcc-bench` e `msvc-profile`; `native-release` já existe e serve para medir a diferença genérico/nativo
  (hazard 4), nunca como número oficial.

### 3.6 Guarda de correção que fica verde durante qualquer otimização

1. **`ctest` completo no preset `gcc-release` apenas** (decisão registrada na memória do autor: só `gcc-release`,
   sem `gcc-debug`): inclui as baterias fuzz (`leverage_fuzz`, `max_return_fuzz`, `ratio_fuzz`,
   `risk_parity_{descent_fuzz,fuzz,reference,solve_sweep}`) e o `style`. Custo conhecido: a parte de RISK_PARITY com
   alvo de volatilidade e sem caixa gasta ≈ 82 s só ela (`validation.md:1016`).
2. **Teste de identidade bit a bit por kernel**: para cada alavanca de risco 0, um teste unitário mantém a
   implementação antiga como referência (cópia em `tests/unit/...`, compilada com `-ffp-contract=off`, como já faz
   `compensated_dot_of`, `CMakeLists.txt:145-148`) e compara `memcmp` em matrizes aleatórias de vários n, inclusive
   condição 1e8 e n que não é múltiplo da largura do laço.
3. **A/B de soluções** (`compare-solutions`) sobre todas as instâncias do §3.2 entre a versão antes e depois: iguais bit
   a bit nas de risco 0.
4. **Sem FMA acidental**: `objdump -d` dos objetos de `portopt_core` num build `native-release` (com `-mfma`
   efetivo): nenhuma instrução `vfm*` fora das funções com `std::fma` declarado.
5. **MSVC**: `msvc-release` compila e a bateria passa (como o VolEngine fez ao fim do 5b); só depois de cada fase,
   não a cada tarefa.
6. **Bandas**: `expect_user_kkt` e `expect_risk_parity_kkt` não mudam; nenhuma banda é afrouxada para acomodar uma
   otimização (regra do plano, não do VolEngine).

---

## 4. Hipóteses de hot spot, por leitura do código

Todas **a confirmar por perfil antes de agir** (regra do VolEngine). Linhas do `HEAD` lido. Modelo de custo
(derivação própria): n variáveis, m linhas explícitas, q restrições ativas, k iterações do GI.

### 4.1 Hot spots, em ordem de suspeita

**H1. Inicialização da base do GI: ≈ n³/2 com passo n, em todo `qp::solve`.**
`goldfarb_idnani.cpp:168` chama `empty_active_set` (`active_set.hpp:57-73`), que para cada coluna c faz
`solve_lower_transposed` (`cholesky.cpp:65-76`), cuja soma interna lê `l[k * n + i]`: acesso de passo n em matriz
row-major. Cada coluna toca n² elementos úteis metade zeros (o vetor unitário tem zeros abaixo de c, e o laço soma
termos `l*0`). EE: em n = 1000, 350 ms contra 105 ms do Cholesky; versão por linhas 42 ms, bits idênticos. Aparece em:
`mean_variance.cpp:86`, `tracking_error.cpp:122`, `homogeneous.cpp:319`, **cada** `solve_with`
(`parametric_path.cpp:277`, chamado de `climb`, `bisect`, `at_minimum`, `max_return.cpp:145,161,184,235`, até 200 vezes
cada laço), `nearest_in_polyhedron` (`risk_parity.cpp:170`) e **cada iteração da descida**
(`polyhedral_descent.cpp:253`, até 200 + 2000 por descida, 50 externas com alvo de volatilidade). Nos problemas de
caminho a matriz G é fixa e J₀ é sempre a mesma (L7).

**H2. Descida proximal do RISK_PARITY: ≈ 3-5 n³ flops por iteração.** Por iteração (`polyhedral_descent.cpp:442-477`):
`model_curvature` → `gauss_newton` (`objective.cpp:187-202`: n³ multiplica-soma, lê e escreve n² doubles n vezes, em n =
1000 são 8 MB reescritos mil vezes, memória-limitado); `factor_model` (`:191-218`: cópia n² e Cholesky n³/3, até 8
tentativas com `tau` crescente); `qp::solve` (H1 mais k iterações do GI, com k ≈ número de cotas ativas, até n);
`matrix_times` (`:95-104`, n²); `armijo` (`:361-381`). `initial_proximal` (`:118-128`) paga **mais um** `gauss_newton`
por `descend`. `evaluate_kkt` (`kkt.cpp:238-248`) chama `linearize(..., true)` (`:242`), que inclui `hessian_of` →
`gauss_newton` e `second_order_term`: cada avaliação KKT custa um n³ a mais, e `newton_polish` avalia ≥ 2 por passo
(`newton_polish.cpp:245-256, 362`). Isto explica o custo medido (9,2 s em n = 50 com cond 1e8): é o **número de
iterações** (55 323, `validation.md:1071`) vezes o custo por iteração, e o primeiro fator só se reduz por otimização.

**H3. Dot2 em build genérico: chamada `libm` por termo.** `compensated_dot` (`objective.cpp:35-50`) chama
`std::fma` (`:42`); EE e compilação confirmam `call fma` sem `-mfma`. `covariance_times` (`:54-63`) faz n desses
produtos de n termos; `linearize` o chama 1 vez (`moments_of`, `:71`), mais `numerator_gradient` (`:122`) quando há
Hessiana, `value_change` 2 vezes (`:209-210`: `moments_of` recomputa `p = Σw` a cada meia-busca de Armijo, embora o
iterado não mude), `evaluate` (`polyhedral_descent.cpp:137`), `relative_gap` (`augmented_lagrangian.cpp:39`),
`volatility_of` (`risk_parity.cpp:151`), KKT (`kkt.cpp:27`). Custo por termo medido: 4,13 ns (genérico) contra 0,63 ns da
soma simples. É O(n²) por produto, então perde para H2 em n grande e pesa em n pequeno e médio, onde há muitas
iterações.

**H4. Varredura de violação com linhas de cota densas.** `add_asset_bounds` (`constraints.cpp:56-71`) cria cada cota
como linha **densa** de n doubles (`unit_row`, `:29-34`): 2n linhas, 2n² doubles (16 MB em n = 1000).
`most_violated_explicit` (`goldfarb_idnani.cpp:191-206`) calcula `slack` (um `dot` de n termos, `:122-125`) para cada
linha, a cada candidato: ≈ 2n² multiplica-somas e 16 MB lidos por candidato, k ≈ n candidatos. `compute_directions`
(`active_set.hpp:157-161`) calcula n produtos `dot(basis_column, row)` de n termos para uma linha que tem **um** termo
não nulo. `gather_constraints` (`:132-155`) ainda copia as linhas (`:137-139`). Para linha unitária o `dot` denso
dá exatamente `x_i` (soma de zeros exatos mais um termo), logo um fast path é bit a bit igual (salvo o sinal de zero).
Este é o argumento de DAQP (§IV-C, verificado: tirar proveito da estrutura de cotas reduz complexidade e memória) e do
qpmad (cotas simples tratadas à parte, verificado no README).

**H5. Custo por iteração do GI, parte linear-álgebra.** `compute_directions` (`active_set.hpp:157-173`): n produtos
n-dimensionais mais `(n - q) · n` axpys; `append_active` (`:175-190`): n - q - 1 rotações de Givens, cada uma sobre 2 colunas
de n (`rotate_basis`, `:92-100`) com um `std::hypot` (`make_givens`, `:75-83`) por rotação; `drop_active`
(`:205-225`) e `refine_toward_active` (`:139-155`, O(q·n), uma vez por candidato final).
`solve_triangle_transposed` (`:116-128`) tem passo n (`triangle[j * n + i]`). Total ≈ 5n² a 7n² flops por iteração,
vetorizável porque as colunas são contíguas (`basis_column`). `set.triangle.assign(n * n, 0.0)` e `basis.assign(n * n)`
(`:61, 68`) alocam e zeram 16 MB por solve em n = 1000 enquanto `triangle` só precisa de `q_max × q_max`.

**H6. Householder com acesso de passo.** `reflect` (`qr.cpp:53-64`) lê `x[i * stride]` com `stride = cols` para cada
coluna j (`:93-97`): n colunas × T linhas de passo n por reflexão. Usado em `tracking_error.cpp:107` (T × n, T ≥ n),
`newton_polish.cpp:218, 342` (matriz (livres + igualdades)² a cada passo de Newton) e `parametric_path.cpp:658-659`
(`multiplier_slope`, n × m, a cada `stays_beyond`/`leave_minimum`). EE: versão por linhas, mesmos bits, 4,0x. Também
`column_norm` (`qr.cpp:13-21`) lê em passo (uma vez, O(T n)). `normal_matrix` (`tracking_error.cpp:27-49`) é T n²/2,
vetorizável (laço j contíguo).

**H7. Cholesky escalar com cadeia de latência.** `factor_in_place` (`cholesky.cpp:12-38`): o laço interno
`sum -= l[i*n+k] * l[j*n+k]` é uma redução em ordem estrita; o GCC sem `-ffast-math` mantém a soma serial (a cadeia de
adições limita a ≈ latência de uma adição por termo; GCC 16 chega a vetorizar só as multiplicações e somar em ordem,
visto em `-fopt-info-vec` na amostra de teste; **não verificado** no código real). EE: 105 ms em n = 1000 (≈ 1,7e8
multiplicações-somas, 0,63 ns cada). Chamado por `Data::create` (`data.cpp:193`, uma vez por `Data`), `spinu.cpp:117` a cada passo de
Newton (≤ 200), `polyhedral_descent.cpp:209` a cada iteração (até 8 vezes), `goldfarb_idnani.cpp:562` quando não há
fator.

**H8. Alocações no laço do GI e nos problemas.** `sign_vector` e `normal` por candidato com família
(`goldfarb_idnani.cpp:258-259`), `append_generated` (`:241-249`), `discard_inactive_generated` (`:389-409`, move e
reconstrói a cada laço, `:518`), `refine_active` (`:413`), `GeneratedRows` crescendo por `insert`. Nos problemas:
cópias de `Σ` e do fator por solve (`mean_variance.cpp:19-29`, `homogeneous.cpp:75-84`; a cópia de `hessian` não é lida
quando há `hessian_factor`, `goldfarb_idnani.cpp:558-560`, mas `validate_sizes` exige `n²`, `:526`), `Path` copia o
fator (`parametric_path.cpp:604`), `hold_face`/`release_face` fazem `insert`/`erase` no meio de vetores de n² doubles
(`:228-264`). Em n pequeno é parte visível; em n grande é ruído (O(n²) contra O(n³)). Confirmar com o contador.

**H9. A família implícita de 2^n linhas não é hot spot.** `most_violated_inequality` (`:253-277`) é O(n) por candidato
(sinal de x, uma linha, um `dot`), mais `duplicates_active` O(q·n). Não otimizar; ver hazard 6.

**H10. Menores.** `model_curvature` faz `gradient[k / n] * gradient[k % n]` por elemento (`polyhedral_descent.cpp:180-183`,
n² divisões e restos inteiros); `jacobian_of` (`objective.cpp:92`) e `second_order_term` (`:147-150`) dividem por
`variance` e `variance²` em cada elemento (três divisões por elemento, içáveis só se a troca for bit a bit igual: `a / v`
**não** é `a * (1/v)`, hazard 1); `gradient_noise` (`kkt.cpp:128-139`) lê `jacobian[i*n+j]` em passo; `volatility_of`
(`homogeneous.cpp:137-152`) lê o fator em passo; `std::hypot` em `make_givens`.

### 4.2 Evidência das fontes por alavanca (a) a (g)

**(a) Goldfarb-Idnani, Powell e implementações abertas.**

- Goldfarb & Idnani, *A numerically stable dual method for solving strictly convex quadratic programs*, Math.
  Programming 27, 1-33, 1983, DOI 10.1007/BF02591962: metadado verificado (Crossref); **texto não lido** (paywall, o
  repositório já registra o mesmo em `2026-11-16-risk-parity-and-numerics.md`). A estrutura de dados do código
  (J, R triangular, rotações de Givens no `append`/`drop`) é a descrita na documentação do próprio `active_set.hpp:22-25`;
  **não verificado** no artigo que o custo por iteração é O(n²) e que o esquema de atualização é o do código; a
  hipótese de custo do §4.1 é **derivação própria** a partir do código.
- Powell, *On the quadratic programming algorithm of Goldfarb and Idnani*, Math. Programming Study 25, 46-61, 1985,
  DOI 10.1007/BFb0121074: metadado verificado (Crossref); texto **não lido** (paywall). Nenhuma afirmação sobre o
  conteúdo é feita aqui.
- **DAQP** (Arnström, Bemporad, Axehill, *A dual active-set solver for embedded quadratic programming using recursive
  LDLᵀ updates*, arXiv:2103.16236, 2021): lido na versão HTML do arXiv. Verificado: usa atualizações recursivas de
  LDLᵀ ao adicionar (Lema 3, nova coluna por um solve triangular) e remover (Lema 4, atualização de posto 1)
  restrições; o foco declarado são QPs **densos** (Observação 10); não dá fórmula de custo por iteração; "explorar a
  estrutura das cotas reduz a complexidade e a memória" (§IV-C); a complexidade de pior caso é determinável offline;
  não discute alocação. **Não é comparável** com o GI direto: usa LDLᵀ em vez da base J, o que muda o custo por
  iteração (O(q²) em vez de O(n²)? **não verificado**, nada disso está no texto lido).
- **qpmad** (Alexander Sherikov, README no GitHub): verificado: é Goldfarb-Idnani sobre Eigen, header-only; fatoração
  do Hessiano **in place**, reaproveitável; opção de guardar o fator inverso no lugar do Hessiano; **cotas simples
  tratadas à parte** das restrições gerais (preferência por ativá-las, que é mais barata); atualiza J inteiro de uma
  vez por julgar mais barato em hardware moderno; três modos de alocação (sob demanda, em compilação por template,
  `reserve()`); comparação com qpOASES/eiQuadProg remetida a outro repositório, sem números na página.
- **eiquadprog** (stack-of-tasks, README): verificado: três variantes, a `-fast` "evita alocação dinâmica" e a `-rt`
  usa tamanho fixo em compilação, "só para problemas pequenos"; sem benchmarks na página. É evidência de que a
  higiene de alocação é a técnica usada onde o GI roda em tempo real (n pequeno), não de ganho em n grande.
- **qpOASES** (Ferreau, Kirches, Potschka, Bock, Math. Prog. Computation 6, 327-363, 2014, DOI
  10.1007/s12532-014-0071-1): metadado verificado; texto **não lido**. Citado só como a referência de método ativo
  paramétrico para o caminho do MAX_RETURN (o PortOpt resolve QPs independentes por `t`, não pelo caminho homotópico).

**(b) Kernels densos sem BLAS.**

- Goto & van de Geijn, *Anatomy of high-performance matrix multiplication*, ACM TOMS 34(3), 2008 (ano do Crossref; o
  PDF lido é o preprint, "Vol. V, No. N, 20YY"): lido o resumo e §4.2: o bloco empacotado de A cabe na **L2**, não na L1,
  e o TLB limita o tamanho do bloco. É a justificativa de blocagem e empacotamento que vale a partir de n ≈ 300 aqui.
  Para os kernels deste código, a lição relevante e **derivação própria**: com acumuladores inicializados no valor
  atual do elemento (`a -= l*l`, em ordem de k), blocagem e tiling por registradores **não mudam** a sequência de
  operações de cada elemento; um GEMM com acumulador zerado (soma do bloco, depois subtração) muda os bits. As
  medições do EE (Cholesky right-looking, base por linhas, reflexão por linhas) confirmam igualdade por `memcmp`.
- BLIS (Van Zee & van de Geijn, ACM TOMS 41(3), 2015, DOI 10.1145/2764454): metadado verificado; texto não lido.
  Citado como referência de micro-kernel de registradores; **sem** afirmação sobre conteúdo.
- Eigen e OpenBLAS: não consultados nesta sessão; nenhuma afirmação sobre projeto interno deles.
- **Autovetorização em build genérico.** GCC 16.2.0, `-O3 -ffp-contract=off`: a redução `s += a[i]*b[i]` não é
  reassociada (sem `-ffast-math`); o `-fopt-info-vec` da amostra mostrou vetorização das multiplicações com soma em
  ordem. O alvo padrão é SSE2 (EE: `-mavx2` e `-mfma` `[disabled]` em `-Q --help=target`).
- **Multiversão por função.** A documentação do GCC (verificada) diz que `target_clones` gera um resolvedor em tempo
  de execução e "depende da extensão `ifunc` do padrão ELF"; **EE verificado**: neste alvo (MinGW-w64, GCC 16.2.0)
  `__attribute__((target_clones("default","avx2")))` falha com *"the call requires 'ifunc', which is not supported by
  this target"*. O despacho manual (`__attribute__((target("avx2,fma")))` com `__builtin_cpu_supports`) compila e
  vetoriza com 32 bytes (EE). Em MSVC não existe o atributo por função; seria uma unidade de tradução compilada com
  `/arch:AVX2` e despacho por `__cpuid` (**não verificado**). Distribuir variantes `-march` quebra o princípio do
  VolEngine ("o binário distribuído é genérico"). O VolEngine resolveu com Google Highway 1.2.0 (despacho em tempo de
  execução); aqui seria dependência nova, contra a spec §10 ("Álgebra linear: própria; sem Eigen").
- **BLAS opcional.** Contra: spec §4.1 (solvers próprios sem dependência, "controle total, wheels triviais") e §10;
  BLAS multithread/ISA-dispatch não garante os mesmos bits entre máquinas, o que quebra a verificação por KKT em
  banda fina e a identidade A/B. Em n ≤ 1000 o ganho esperado sobre kernels próprios bem escritos não foi medido.
  Veredito: **não**.

**(c) Determinismo de ponto flutuante.**

- GCC, `-ffp-contract` (documentação verificada): `off` desliga a contração; `fast` forma FMA "se o alvo tem suporte
  nativo"; **o padrão é `off` para C em modo padrão (`-std=c11`) e `fast` caso contrário**; contração é dentro de uma
  expressão, não entre instruções. Logo, em C++ (`-std=c++23`) o padrão é `fast`, como o `CMakeLists.txt:63-67` do
  PortOpt já registra: só vira FMA quando o alvo tem FMA, e o flag `off` do core é o que impede isso em
  `-march=native`.
- MSVC (documentação verificada): `/fp:precise` é o padrão e preserva a ordem das expressões; **a partir do Visual
  Studio 2022 não gera contrações por padrão em `/fp:precise`; antes, podia**. `/fp:fast` e `/fp:contract` habilitam. O
  comentário do `CMakeLists.txt:63` ("MSVC's default /fp:precise does not") vale só para VS ≥ 2022; o PortOpt usa o
  gerador "Visual Studio 18 2026", então vale; registrar como pré-requisito.
- `std::fma` é **definido** pelo IEEE 754 como multiplicação-soma com um só arredondamento, e produz os mesmos bits em
  qualquer máquina (em hardware ou em `libm`): uso explícito de `std::fma` num kernel é determinístico entre ISAs,
  diferente da contração do compilador, que depende do alvo. Isso permite, **por decisão explícita e fora do escopo de
  risco 0**, kernels com `std::fma` (mais exatos e mais rápidos com FMA de hardware); em build genérico, sem despacho,
  seriam mais lentos (chamada `libm`, H3). Em arm64 (GCC/Clang, plataforma da spec §11) o FMA é de linha de base e
  `std::fma` inlina; o problema de H3 é só x86-64 (**não verificado** em arm64).
- Ogita, Rump & Oishi, *Accurate Sum and Dot Product*, SIAM J. Sci. Comput. 26(6), 1955-1988, 2005, DOI
  10.1137/030601818: metadado verificado; o texto não foi lido (403). O `compensated_dot` (comentário do código,
  `objective.cpp:31-34`) já cita o Dot2; a equivalência entre o TwoProduct por FMA e o de **Dekker** (Dekker, *A
  floating-point technique for extending the available precision*, Numer. Math. 18, 224-242, 1971, DOI
  10.1007/BF01397083, metadado verificado, texto não lido) é conhecimento de memória (**não verificado** nas
  fontes), **verificado por EE**: 0 diferenças de bits em 2000 produtos de ponto de 1000 termos com magnitudes
  variando por 2⁴⁰. As condições de exatidão de Dekker (sem overflow no `split`, sem subnormais no erro) precisam ser
  tratadas por teste de bordas.
- **O que precisa ficar verde**: os testes `expect_user_kkt`, `expect_risk_parity_kkt` e as baterias fuzz, mais a
  identidade bit a bit do §3.6. A banda de RISK_PARITY tem margem 1e3 em geral, mas **6** no caso
  cond 1e8, n = 50, sem caixa (`validation.md:804-805`): é o caso em que uma mudança de arredondamento (risco B) pode
  quebrar a banda; ver hazard 1.

**(d) Prática de Google Benchmark e perfilamento.**

- User Guide (verificado): `DoNotOptimize` força o resultado em memória/registrador e age como barreira de leitura e
  escrita em compiladores GNU, mas **não impede que a expressão seja otimizada** (guardar em variável local primeiro);
  `ClobberMemory` força as escritas pendentes. `--benchmark_repetitions`, `ReportAggregatesOnly`,
  `DisplayAggregatesOnly` (só a saída de console); contadores com `kIsRate`, `kAvgThreads`; `ComputeStatistics`;
  `Complexity`/`SetComplexityN`/`BigO` (ajuste de curva por família: útil para confirmar n³); `RegisterMemoryManager`
  (contagem de alocações, **só no JSON**); o número de iterações é dinâmico (≥ 1, ≤ 1e9, até 5x o tempo mínimo);
  `--benchmark_min_time` aceita `Ns` (segundos) ou `Nx` (iterações).
- *Reducing Variance* (verificado): escalonamento de frequência por `cpupower frequency-set --governor performance`
  e `boost` desligável; fixar a CPU por `taskset -c 0`; prioridade por `nice`/`chrt`. **A página não cobre Windows**
  nem `--benchmark_enable_random_interleaving` (o VolEngine a usa em `benchmarks.md` Plano 6b; a existência da flag
  não foi confirmada na página lida). No Windows o VolEngine roda sem isolamento e aceita o desvio; recomenda-se
  intercalar A/B e repetir.
- `compare.py` (`docs/tools.md`, verificado): modos `benchmarks`, `filters`, `benchmarksfiltered`; requer `scipy`; teste
  U de Mann-Whitney com p-valor contra α; "**no mínimo 9 repetições** para o teste U ser significativo"; mostra
  Time/CPU old/new, variação relativa `(new - old)/|old|` e `OVERALL_GEOMEAN`. O `quick.py` do VolEngine usa 3
  repetições e uma regra de ruído por desvio, **mais fraca** que o U test, adequada ao laço rápido; para decidir uma
  alavanca de ganho pequeno (< 5%), usar `compare.py` com ≥ 9 repetições, ou repetir a rodada intercalada.
- **AMD uProf** (resultado de busca e a receita do VolEngine): `AMDuProfCLI collect --config tbp` (amostragem
  temporal, suportada no Windows) e `report` com caminhos de símbolos; `timechart` é subcomando à parte (potência,
  frequência). A página do guia não foi lida em profundidade (**não verificado**); vale a receita local do VolEngine.
  `perf` no Linux e VTune no Windows Intel: não consultados.

**(e) Memória, alocação e layout.**

- Reuso de workspace e alocação única: regra do STYLEGUIDE do VolEngine/PortOpt (verificada nos arquivos). Evidência
  externa: eiquadprog-fast e qpmad (acima). `std::pmr` e *small-buffer*: não consultados; para este código o ganho
  vem de tirar alocações do laço (H8), não de trocar o alocador.
- **Packed contra completo.** O núcleo usa matriz completa row-major; o empacotado economiza memória (÷ 2) mas perde
  a contiguidade de linha que a vetorização usa, e todos os códigos leem Σ por linha. Veredito: manter completo; tirar
  proveito da **simetria no cálculo** (L4: `J'J` calculado pela metade e espelhado, com bits idênticos por derivação
  própria: `(2a)·b = 2(ab) = (2b)·a` exatos em ponto flutuante binário sem overflow, e a soma em i é a mesma).
- **Triangular.** `set.triangle` é `n × n` mas só `q × q` é usada (H5); `basis` pode guardar só a parte triangular
  superior de L⁻ᵀ nas colunas ainda não rotacionadas (J₀ é triangular superior; a estrutura se perde nas rotações).

**(f) Paralelismo.**

- **Lote** (`po_solve_batch`, spec §5.4: "laço de `po_solve` no pool, sem caminho duplicado; os dados são copiados no
  handle, nada no contrato precisa mudar"). Auditoria por leitura: não há estado global mutável no `src/` além de
  `thread_local std::string last_error_message` em `capi/error.cpp:11`; `Data` é só leitura durante o solve (spec
  §5.2: "`po_solve` pode ser chamado concorrentemente com handles diferentes ou com os mesmos handles somente
  leitura"). Como cada solve aloca o próprio espaço de trabalho e usa só entradas imutáveis, o resultado é **idêntico
  para qualquer número de threads** (regra do VolEngine, "bit-identical for any number of threads"); medir por
  `BM_Batch/threads:N` análogo ao do VolEngine (speedup e `pool_use`). O pool: o VolEngine usa `bs_thread_pool` v5.0.0
  por `FetchContent`; dependência nova que a spec §10 precisaria listar.
- **Paralelismo dentro de uma solução**: reduções (produtos matriz-vetor, `J'J`, Cholesky por blocos) mudam a ordem de
  soma conforme o número de threads, quebrando a reprodutibilidade, a menos que a árvore de redução seja fixa
  (independente de threads) ou se use soma reprodutível (Demmel & Nguyen, *Parallel Reproducible Summation*, IEEE
  Trans. Computers 64(7), 2060-2070, 2015: afirmação **não verificada**, vem de resultado de busca; o artigo não foi
  lido). Não vale na v1.

**(g) Python/ctypes contra o tempo de solve (modelo: `python-vs-arch.md` do VolEngine).**

- O VolEngine mede trabalho equivalente contra `arch` e documenta "otimizadores diferentes". Aqui não há comparação
  de biblioteca; o análogo é separar custo fixo de chamada e custo de solve. Dado exploratório: uma chamada ctypes
  trivial com `argtypes` (`msvcrt.abs`, Python 3.14.7) custou ≈ 0,38 µs; `po_solve` passa pelo `create_handle`
  (`python/src/portopt/_result.py:70`), as conversões de `numpy` (`np.ascontiguousarray`, `_spec.py:38-39`) e a leitura
  do resultado vetor a vetor (`read_vector`, `_result.py:44-46`). Hipótese (**não medida**): em n = 10 o solve leva
  dezenas de µs e as camadas Python somam µs a dezenas de µs, mesma ordem; em n ≥ 100 o solve domina. Medir no §3.4.
- `Data::create` (Cholesky, `data.cpp:193`) é pago **uma vez por `Data`**, e todo solve reaproveita o fator
  (`mean_variance.cpp:25`): o padrão de uso "uma `Data`, vários `Spec`" amortiza o n³/3; é o que o frontier em Python
  (`spec §6`) faz.

---

## 5. Hazards

1. **Reassociar muda bits, e a banda é fina em um caso.** Kernels de classe B (multi-acumulador, GEMM com acumulador
   zerado, SIMD de redução) alteram o arredondamento por 1e-16 relativo por operação; amplificado por cond(Σ) até 1e8, o
   erro de peso chega a 1e-5 a 1e-4 em RISK_PARITY com vol-alvo (`validation.md`, seção do conditioning limit), e a
   margem para a banda é **6** em cond 1e8, n = 50, sem caixa (`:804-805`). Corolário: `a / v` não pode virar
   `a * (1 / v)` (muda bits), `std::hypot` não pode virar `sqrt(a*a + b*b)` (muda bits e overflow), e a ordem de
   soma de **cada elemento** é parte do contrato das alavancas de risco 0.
2. **`-ffp-contract=off` é requisito do código, não detalhe de build.** Todo alvo novo que compile cópia de kernel ou
   TU com Dot2 (executável de dump, biblioteca de variante ISA, testes) precisa do flag; sem ele, em `-march` com FMA,
   o Dot2 deixa de ser livre de erro (o erro de TwoSum/TwoProduct deixa de ser exato) e os resultados mudam de forma
   silenciosa. O oráculo de teste já duplica o flag (`CMakeLists.txt:145-148`).
3. **Dekker no lugar de `std::fma` só é exato sem overflow e sem subnormais.** O `split` multiplica por 2²⁷+1 (EE:
   `134217729.0`): falha para |a| > 2⁹⁹⁶ e perde exatidão do erro em subnormais. Tratar por teste de bordas ou por
   desvio ao caminho `std::fma` nesses intervalos; os pesos e Σ do domínio ficam longe, mas a garantia deve ser
   escrita, não suposta.
4. **Genérico contra nativo.** O `std::fma` vira chamada de `libm` só em build genérico (H3); o `native-release`
   esconde o custo. Números oficiais vêm do `gcc-bench` genérico, como no VolEngine; os wheels são genéricos. Qualquer
   ganho "de ISA" (L6 segunda etapa, L12) precisa de despacho em tempo de execução, e em MinGW **não há `ifunc`**
   (EE), logo o despacho é manual.
5. **MSVC.** Código MSVC é mais lento que o GCC (VolEngine: 20 s contra 7,7 s) e o uProf só funciona nele: perfil
   indicativo. `/fp:precise` sem contração só em VS ≥ 2022 (verificado). O refinamento com `long double` mais largo que
   `double` nos testes de RISK_PARITY é pulado em MSVC (`validation.md:808-809`): a guarda é mais fraca lá.
6. **A família implícita de 2^n linhas.** Custo O(n) por candidato (H9); mudar a regra de empate, o `sign_vector`
   (`x_i = 0` conta como +1, `:210-221`) ou a ordem de geração/descarte (`discard_inactive_generated`) altera quais
   linhas geradas ficam ativas, e `same_segment` do caminho compara `gross_exposure_signs` por igualdade exata
   (`parametric_path.cpp:627-631`): a otimização quebraria os segmentos do MAX_RETURN e as medidas de 1992 solves de
   `validation.md:299`. Não otimizar por perfil, só se um perfil o mostrar caro.
7. **Segurança de threads de caches e workspaces.** `Data` é compartilhável entre threads em leitura (spec §5.2): um
   cache de J₀ dentro de `Data` exigiria `std::call_once` ou cálculo no `create`, e +8 MB em n = 1000 por `Data`. Preferir
   cache local à função que resolve (o `Path` do MAX_RETURN) ou parâmetro explícito de `qp::solve`: sem estado
   compartilhado, sem hazard.
8. **Warm start e redução de esforço (L13, L18) não são otimizações.** Mudam o caminho de iteração, a contagem
   `iterations` (que é parte do `Result`) e o limite de iteração; podem trocar qual mínimo local o RISK_PARITY encontra
   (a validação registra 6 de 91 casos em mínimos locais). O estudo de 5a do VolEngine (memória L-BFGS 40) aceitou isso
   com os números na mão do autor; aqui exigiria a mesma decisão explícita e a varredura de 576 formas
   (`RiskParitySolveSweep`).
9. **Benchmark que mede o caminho errado.** Um caso que acaba em `MaxIterations` ou `Infeasible` é rápido e inútil:
   `SkipWithError` se `status != Optimal` (regra de falhar alto). Instâncias com cond 1e8 sem caixa têm RISK_PARITY com
   alvo de 9-20 s: ficam fora da suíte rápida, com `Iterations(1)`.
10. **Gerar a instância não é o benchmark.** `random_covariance` é O(n³) com `reflect` em passo (código de teste,
    `random_covariance.hpp:47-122`): n = 1000 leva segundos; gerar uma vez por processo, fora do laço cronometrado, e
    manter as sementes fixas para o A/B comparar a mesma instância.
11. **Alocação contada por `operator new` substituído** funciona no executável de benchmark (core estático, libstdc++
    estática) mas não combina com sanitizers; o `RegisterMemoryManager` do Google Benchmark só escreve no JSON.
12. **O gargalo se desloca.** Depois de L1 o termo n³/2 some e a iteração do GI (H5) e o `gauss_newton` (H2) passam a
    dominar: reperfilar a cada tarefa; não encadear alavancas por uma ordem fixa no papel.
13. **Fontes pagas não lidas** (GI 1983, Powell 1985, Ogita-Rump-Oishi 2005, Dekker 1971): nenhuma decisão do plano
    deve depender de afirmação sobre o conteúdo deles sem leitura.

---

## 6. Divisão de tarefas proposta

Cada tarefa é commitável sozinha e passa pela guarda do §3.6. Tarefas de otimização: perfil → A/B → leitura do
assembly → número no `docs/performance/` → commit (regra do STYLEGUIDE).

| # | Tarefa | Saída |
|---|---|---|
| 0 | **Infra**: bloco `PO_BUILD_BENCHMARKS` no `CMakeLists.txt` (Google Benchmark v1.9.1, mesmo hash), alvos `portopt_benchmarks` e `portopt_profile`, presets `gcc-bench` e `msvc-profile`, `quick.py` (`run`, `compare`, `compare-solutions`, `selftest`) | `tests/benchmarks/`, presets |
| 1 | **Instâncias e `BM_Solve`**: geradores com cache, casos do §3.2, afirmação de status, contadores, filtro da rápida | suíte rodando; tempo da rápida medido |
| 2 | **`BM_Kernel`** (§3.3) e contador de alocações (`operator new` substituído); decidir `PO_BENCH_INSTRUMENT` para `qp_solves` | micro-suíte |
| 3 | **`portopt_profile` e o dump de soluções** (hex-float), mais teste da própria comparação (`selftest`) | ferramenta de A/B de soluções |
| 4 | **Relatório de linha de base** `docs/performance/baseline.md` (+ `measuring.md` e `benchmarks.md` portados): máquina, build `gcc-bench` genérico, 5 repetições, tabela por objetivo × N; perfil uProf (`msvc-profile`) dos casos lentos; **confirma ou refuta H1-H10** | `docs/performance/` |
| 5 | **L1** base do GI por linhas (+ teste `memcmp` contra a versão antiga) | A/B, commit |
| 6 | **L3** reflexão de Householder por linhas | A/B, commit |
| 7 | **L2** Cholesky right-looking / registradores, ordem por elemento (+ teste `memcmp`) | A/B, commit |
| 8 | **L4** `gauss_newton` por simetria e blocagem | A/B, commit |
| 9 | **L5** fast path de linha unitária (decisão: sinalizador `row_is_unit` em `gather_constraints`, `compute_directions` com índice) | A/B, commit |
| 10 | **L7** cache de J₀ no caminho do MAX_RETURN (e solves com o mesmo fator) | A/B, commit |
| 11 | **L6a** Dot2 por Dekker em build genérico (+ teste de bordas e `memcmp` contra `std::fma`) | A/B, commit |
| 12 | **L8 e L9** higiene de alocação no GI e nos problemas, cache de `Moments` na busca de Armijo (só se H8 e H3 mostrarem no perfil) | A/B, commit |
| 13 | **Reperfil** e relatório intermediário; decidir L6b (despacho manual `target("avx2,fma")`), L10 (`-flto`), L12 | `profile-N.md` |
| 14 | **Camada C e Python**: `BM_CApi`, `bench_overhead.py`, linha no relatório | medida do custo fixo |
| 15 | **Gate de v1**: relatório final `docs/performance/benchmarks.md`, tabela antes/depois, MSVC `msvc-release` verde, rodada `native-release` só para registrar a diferença | v1 |

Fora deste plano (planos próprios, não otimizações): lote com threads (`po_solve_batch`, spec §5.4), warm start do
QP no RISK_PARITY (L13), mudança de esforço do RISK_PARITY (L18), kernels de classe B (L11) e multiversão por ISA (L12)
se o perfil da tarefa 13 os justificar. Cada um precisa de decisão explícita do autor sobre risco numérico.

---

## 7. Bibliografia e estado de verificação

| Fonte | URL | Estado |
|---|---|---|
| VolEngine `docs/performance/*.md`, `tests/benchmarks/*`, `CMakePresets.json`, `CMakeLists.txt`, `ARCHITECTURE.md`, `STYLEGUIDE.md` | locais (`C:/Users/x/Documents/Programs/Cpp/VolEngine`) | **Verificado** (lidos) |
| PortOpt `src/`, `CMakeLists.txt`, spec, `validation.md` | locais | **Verificado** (lidos; `validation.md` não é UTF-8, trechos lidos com U+FFFD, números e linhas conferidos) |
| Google Benchmark, *User Guide* | https://raw.githubusercontent.com/google/benchmark/main/docs/user_guide.md | **Verificado** (pelo extrator; seções citadas) |
| Google Benchmark, *Reducing Variance* | https://raw.githubusercontent.com/google/benchmark/main/docs/reducing_variance.md | **Verificado**; não cobre Windows nem interleaving |
| Google Benchmark, *Tools* (`compare.py`) | https://raw.githubusercontent.com/google/benchmark/main/docs/tools.md | **Verificado**; flag `--alpha` não documentada na página |
| GCC, *Optimize Options* (`-ffp-contract`) | https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html | **Verificado** (parte lida; o final da página, com `-ffast-math`, não foi lido) |
| GCC, *Function Multiversioning* | https://gcc.gnu.org/onlinedocs/gcc/Function-Multiversioning.html | **Verificado**; lista de SO não dada na página; ausência de `ifunc` no MinGW **verificada por compilação** (EE) |
| Microsoft, `/fp` | https://learn.microsoft.com/en-us/cpp/build/reference/fp-specify-floating-point-behavior | **Verificado** |
| Arnström, Bemporad, Axehill, DAQP, arXiv:2103.16236 (2021) | https://arxiv.org/abs/2103.16236 (HTML: https://ar5iv.labs.arxiv.org/html/2103.16236) | **Verificado** (HTML); sem fórmula de custo por iteração |
| qpmad (Sherikov) | https://github.com/asherikov/qpmad | **Verificado** (README) |
| eiquadprog (stack-of-tasks) | https://github.com/stack-of-tasks/eiquadprog | **Verificado** (README) |
| Goto & van de Geijn, *Anatomy of high-performance matrix multiplication*, ACM TOMS 34(3), 2008 | https://www.cs.utexas.edu/~flame/pubs/GotoTOMS_revision.pdf; DOI 10.1145/1356052.1356053 | **Verificado** (resumo e §4.2 do preprint; ano e volume por Crossref) |
| Goldfarb & Idnani, Math. Programming 27, 1-33, 1983 | https://doi.org/10.1007/BF02591962 | Metadado **verificado** (Crossref); texto **não verificado** (paywall) |
| Powell, Math. Programming Study 25, 46-61, 1985 | https://doi.org/10.1007/BFb0121074 | Metadado **verificado**; texto **não verificado** |
| Ogita, Rump, Oishi, SIAM J. Sci. Comput. 26(6), 1955-1988, 2005 | https://doi.org/10.1137/030601818 | Metadado **verificado**; texto **não verificado** (403) |
| Dekker, Numer. Math. 18, 224-242, 1971 | https://doi.org/10.1007/BF01397083 | Metadado **verificado**; texto **não verificado** |
| Ferreau et al., qpOASES, Math. Prog. Comp. 6, 327-363, 2014 | https://doi.org/10.1007/s12532-014-0071-1 | Metadado **verificado**; texto **não verificado** |
| Van Zee & van de Geijn, BLIS, ACM TOMS 41(3), 2015 | https://doi.org/10.1145/2764454 | Metadado **verificado**; texto **não verificado** |
| Demmel & Nguyen, *Parallel Reproducible Summation*, IEEE Trans. Computers 64(7), 2060-2070, 2015 | DOI 10.1109/TC.2014.2345391 (de resultado de busca) | **Não verificado** (nem o DOI nem o conteúdo) |
| AMD uProf, guia do usuário | https://docs.amd.com/r/en-US/57368-uProf-user-guide/ | **Não verificado** (só resultado de busca; vale a receita do VolEngine) |
| Experimentos exploratórios desta sessão (Dot2 fma contra Dekker contra simples; Cholesky, base e reflexão por linhas; `target_clones`; ctypes) | scratchpad da sessão, fora do repositório | **Derivação própria / EE**: uma execução, sem isolamento; refazer no harness |

### Números dos experimentos exploratórios (referência única, não linha de base)

| Experimento | Resultado |
|---|---|
| Dot2, n = 1000, L1, genérico (`call fma`) | 4,13 ns/elem; versão Dekker 2,20 ns/elem; soma simples 0,63 ns/elem |
| Dot2, `-mavx2 -mfma` (FMA inline) | 1,11 ns/elem; Dekker 2,19; simples 0,62 |
| Dot2 por FMA contra Dekker, bits | 0 diferenças em 2000 casos (também com `-mfma`) |
| Cholesky atual contra right-looking (n = 100, 300, 1000) | 0,07/0,06 ms; 2,32/1,56 ms; 105/53 ms; triângulo inferior bit a bit igual |
| Base J: atual contra por linhas (n = 100, 300, 1000) | 0,40/0,04 ms; 9,55/0,99 ms; 350/42 ms; bit a bit igual |
| Reflexão (20 aplicações) 120x100, 600x300, 1500x1000 | passo 0,18/5,08/42,6 ms; linhas 0,09/1,19/10,8 ms; bit a bit igual |
| `target_clones` em MinGW GCC 16.2.0 | erro: *requires 'ifunc', which is not supported by this target* |
| ctypes, chamada trivial com `argtypes`, Python 3.14.7 | ≈ 382 ns |
