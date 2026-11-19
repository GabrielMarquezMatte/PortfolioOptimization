# PortOpt — design dos bits portáveis e do cruzamento do MAX_RETURN (plano 8)

Data: 2026-11-18. Complementa `2026-11-15-portopt-design.md` (o spec da biblioteca) e
`2026-11-17-portopt-v0-release-design.md` (plano 7). Pesquisa de apoio:
`docs/research/2026-11-18-max-return-path-and-portable-generators.md` (citada como "pesquisa §N").

## 1. Intenção

A v0.1.0-rc.1 saiu (GitHub Release e PyPI), mas o `ci.yml` está vermelho no Linux e no macOS: seis testes que passam
no Windows falham lá (CHANGELOG, "Known limitations"). Duas causas, e o autor quer as duas resolvidas:

1. **Os testes não sorteiam o mesmo problema em toda plataforma.** Os geradores usam `std::log`/`exp`/`cos`/`pow` e
   `std::uniform_*_distribution`, cujos bits mudam entre glibc, UCRT, MinGW e Apple libm, e entre libstdc++, libc++ e
   MSVC STL (pesquisa §4.1, §4.4). A própria biblioteca chama `std::hypot` (`src/qp/active_set.hpp:132`, Givens do QP)
   e `std::log` (`src/riskparity/spinu.cpp:72`), então também ela pode diferir no último bit entre plataformas.
2. **O `MAX_RETURN` é frágil num caso real.** O shape 3963545680699380072 (n 24, cond(Σ) = 1e8, médias quase iguais,
   cash, leverage) vira `MAX_ITER` sob uma perturbação de 3,6e-13 relativa de Σ, em qualquer plataforma; no Windows ele
   só passa porque está pinado com uma banda relaxada (`kKnownResidual`). Causa (pesquisa §0, §1): o cruzamento
   `vol = s` é interpolado pela **diferença** de duas soluções do Goldfarb–Idnani, cujo arredondamento (o GI parte de
   `t Σ⁻¹ v`, ~7e8 no caso) gera 1,5e-7 de ruído em vol(Q(t)); e `same_segment` compara **vetores de sinal**, que
   alternam entre duas representações válidas da mesma face, então o bracket nunca "fica num conjunto ativo".

**Por que agora:** a v0.1.0 final espera este plano (decisão do autor). Sai com o CI verde nas três plataformas, os
mesmos bits de solução em todas elas e o `MAX_RETURN` sem caso pinado.

## 2. Objetivos, não-objetivos, anti-objetivos

### Objetivos

1. **Biblioteca sem libm transcendental.** `src/` só usa `+ − × ÷`, `sqrt` e `fma` (operações que o IEEE 754 exige
   corretamente arredondadas). `hypot` vira uma versão própria (escala por potência de 2, `sqrt`, desescala); o `log` do
   Spinu usa o `log` portátil. **Condição do autor:** sem perda de performance (§6.3).
2. **Geradores de teste portáveis.** A mesma semente dá os mesmos bits de μ, Σ e de todo número sorteado em Windows
   (MinGW e MSVC), Linux (GCC) e macOS (Apple Clang, arm64), conferido por teste de valores-ouro.
3. **Mesmos bits de solução em toda plataforma**, conferido por um teste de soluções-ouro (um caso por objetivo).
4. **Cruzamento do `MAX_RETURN` por tangente exata com certificado KKT** (pesquisa §5a), removendo `kAdjacentEnd`,
   `kCrossingSanity` e o pin `kKnownResidual`.
5. **CI verde** (`ci.yml`, `cpp.yml`, `wheels.yml`) nas três plataformas, sem teste desabilitado nem pinado para isso.
6. `validation.md`, CHANGELOG e o resto da documentação medidos de novo e coerentes com as novas instâncias.

### Não-objetivos

- Tags de release (rc.2 e v0.1.0 são do autor, ao fim).
- Arredondamento correto (CORE-MATH): basta determinismo com erro ≤ 1 ulp (pesquisa §4.2).
- Mexer nos outros objetivos além do que a troca de `hypot`/`log` muda de arredondamento.
- Homotopia estilo qpOASES (ratio test entre faces), no caminho ou no cruzamento: a bissecção sobre pontos do GI
  continua decidindo o bracket; a tangente entra só no cruzamento, numa face.
- Windows arm64 e outras plataformas fora da matriz do plano 7.

### Anti-objetivos

- **Esconder falha de plataforma:** `ctest -E`, `GTEST_SKIP` por plataforma, `#ifdef` de SO num teste ou banda
  afrouxada para o CI ficar verde. Nenhuma banda de `expect_user_kkt`/`expect_risk_parity_kkt` afrouxa.
- **Novo pin por nome.** Um caso que falha com as instâncias novas é defeito a explicar e corrigir, não a fixar.
- **Regressão de performance** pela troca da libm (condição do autor).
- **Duas cópias da matemática portável** (uma em `src/`, outra em `tests/`).

## 3. Restrições

- As dos planos 1–7: `STYLEGUIDE.md`, commits sem `Co-Authored-By`, `-ffp-contract=off` em todo alvo que compila
  código de produção **e nos testes** (já está), sem `long double` nem `-ffast-math`, `docs/superpowers/` em português,
  documentos do repositório em inglês.
- Testes só no preset `gcc-release` localmente (preferência do autor); MSVC (`native-release`) também, porque o
  objetivo 3 é justamente MinGW = MSVC = Linux = macOS.
- Performance: método de `docs/performance/measuring.md` (A/B na mesma sessão, `quick.py run` duas vezes antes).
- Licença: o porte do fdlibm mantém o aviso da Sun no topo do arquivo (exigência da licença, compatível com MIT).

## 4. Abordagens consideradas

| Parte | Abordagem | Troca | Decisão |
|---|---|---|---|
| Geradores | **1A. Matemática própria** (`exp`/`log` portados do fdlibm 5.3 com `std::bit_cast`; normais pelo método polar de Marsaglia, sem sin/cos; inteiros e reais uniformes a partir de `engine()` cru) | Mesmo formato de instância de hoje; ~150 linhas + guarda de bits | **Escolhida pelo autor** |
| Geradores | 1B. Geradores só com aritmética (normais de Irwin–Hall etc.) | Muda a distribuição; espectro log-uniforme ainda pede `pow` | Descartada |
| Geradores | 1C. Instâncias gravadas em arquivo | Arquivos grandes; o fuzz perde amplitude | Descartada |
| MAX_RETURN | **2A. Tangente exata a partir de um ponto** (null-space da face, ancorada num solve de face em t_a, referencial do shift) + certificado KKT num ponto | Muda bits só nos casos que chegam ao cruzamento; mantém o caminho | **Escolhida pelo autor** |
| MAX_RETURN | 2B. Solver de reserva (SOCP/Newton) quando o caminho quebra | Segundo método no núcleo | Descartada |
| MAX_RETURN | 2C. Reparametrizar o caminho | Muda todo bit do MAX_RETURN; sem garantia no cond 1e8 | Descartada |
| `src/` libm | **Tirar** `hypot` e `log` de `src/` | Mesmos bits em toda plataforma; muda bits do QP | **Escolhida pelo autor, se não custar performance** |
| `src/` libm | Deixar e registrar | Fuzz pode divergir no último bit entre CIs | Descartada |

## 5. Decisões já tomadas

- As duas partes entram, e a v0.1.0 final espera o plano (autor).
- 1A + 2A (autor). Tirar a libm de `src/`, condicionado a performance (autor).
- **Uma cópia só da matemática portável, em `src/`** (decisão do Claude): um header interno (ex.
  `src/core/portable_math.hpp`, namespace `po::portable`) com `exp`, `log`, `hypot`; os testes o incluem (já incluem
  headers de `src/`). Funções que só os testes usam (`exp10`, normal polar, uniformes) ficam em
  `tests/unit/random_covariance.hpp` ou num header de testes ao lado, por cima das de `src/`.
- **Guarda de bits da biblioteca** (decisão do Claude): o teste de soluções-ouro (§6.4) com os bits em hex-float
  gravados no repositório. Uma mudança legítima de bits os regrava no mesmo commit, com o `compare-solutions` no corpo
  do commit.
- **Fixture da Σ do Linux** (decisão do Claude): a Σ do shape 3963545680699380072 como o gerador atual a produz no
  Linux (extraída por Docker no `manylinux_2_28` com GCC 15, ou de um artefato de CI) gravada em hex-float num arquivo
  de teste, com o μ e as frações 0,3–0,9. É o teste de regressão do mecanismo, independente dos geradores, e continua
  válido depois que eles mudarem. Idem a Σ do Windows (as duas versões do mesmo problema).
- **Ordem** (decisão do Claude): cruzamento do MAX_RETURN → libm fora de `src/` → geradores portáveis → CI. O
  cruzamento é validado contra os números atuais do caso pinado antes de as instâncias mudarem; `validation.md` é
  remedido uma vez só, com a biblioteca final.

## 6. Desenho

### 6.1 Cruzamento do MAX_RETURN (pesquisa §1, §2, §5a)

1. **Segmento pela face.** `same_segment` compara linhas explícitas ativas e, para a família da leverage, se ela está
   ativa, o suporte com sinais e o conjunto de zeros D; os vetores de sinal deixam de entrar.
2. **Solve de face** (rotina nova em `parametric_path.cpp`, ao lado de `multiplier_slope`): para a face e o shift de um
   ponto, N nas coordenadas livres (zeros eliminados exatamente), QR de Householder de N (`linalg::householder_qr`),
   base Z do espaço nulo, QR de LᵀZ com o fator de `Data`; devolve `dw = Z H_r⁻¹ Zᵀ v` e `w(t)` num t dado, e os
   multiplicadores por mínimos quadrados em `N u = Σw − t v`, no referencial do shift. Não resolver o KKT inteiro por
   eliminação: cond ≈ 1e13, erro 2e-8 a 6e-8 em dw (pesquisa §0).
3. **Cruzamento** (fluxo decidido pelo Claude). A bissecção sobre pontos do GI continua como hoje. Quando as duas
   pontas do bracket estão na **mesma face** (item 1), `w_a` = solve de face em t_a; raiz de `‖Lᵀ(w_a + h dw)‖ = s`
   pela fórmula sem cancelamento de `crossing_fraction`; t* = t_a + h; certificar (item 4). Certificou: resposta.
   Não certificou: a bissecção segue. Pontas de faces diferentes: a bissecção segue, como hoje.
4. **Certificado** num ponto (o QP é convexo, KKT basta; pesquisa §2): primal com a tolerância do GI; dual das linhas
   explícitas; dual da família na forma invariante à representação (`|g_i| ≤ 1` nos zeros, ζ ≥ 0), tudo no referencial
   do shift, com tolerância dual `kNoiseFactor·n·ε·(t*|v|max + |Σ|max|w|max)`.
5. **Em doubles adjacentes** (o que hoje é `at_adjacent`): tentar o candidato da face de t_a, depois o da de t_b;
   nenhum certifica: `MAX_ITER`. **Sem homotopia com ratio test** (decisão do Claude): nenhum caso conhecido chega lá
   (com a face como identidade, o caso medido fica numa face só; pesquisa §2) e seria o maior pedaço de código sem
   teste que o exercite. Se a triagem da §6.5 revelar um caso que precise dela, é um defeito real do solver (§6.5 (c)).
6. **Resultado** como hoje: multiplicadores do usuário `u/t*`, `λ_vol = −s/t*`.
7. **Remoções:** `kAdjacentEnd`, `kCrossingSanity`, `at_adjacent` no que ficar morto, `kKnownResidual` e
   `pinned_against_reference`. Os comentários de `validation.md` que justificavam essas constantes saem ou viram
   história.

**Aceitação:** o caso pinado `OPTIMAL` com `vol/s − 1` na faixa de ε, retorno igual ao da referência e
`expect_user_kkt` na banda **cheia**; a fixture da Σ do Linux `OPTIMAL` nas frações 0,3–0,9 (hoje `MAX_ITER`); a
bateria de fuzz do MAX_RETURN inteira na banda cheia, sem pin. `compare-solutions` contra o antes: só casos de
`max_return` mudam, e cada mudança está dentro da banda.

### 6.2 Matemática portável (pesquisa §4, §5b)

- `exp`, `log`: porte do fdlibm 5.3 (`e_exp.c`, `e_log.c`), `std::bit_cast` no lugar do type punning, constantes em hex
  literal, domínio restrito ao que o código usa (finito; x > 0 no log) com `assert` fora dele. Erro medido ≤ 0,80 ulp.
- `hypot(a, b)`: escala por potência de 2 escolhida pelo expoente do maior, `sqrt(a² + b²)`, desescala; zeros e
  subnormais tratados. Mesmo contrato que `make_givens` usa hoje (sem overflow/underflow intermediário).
- Testes: `exp10(x) = exp(x · ln10)` (≤ ~35 ulp, suficiente para spread log-uniforme); `standard_normal` pelo método
  polar de Marsaglia (só `log` e `sqrt`); `uniform_index(n)` por rejeição sobre `engine()`; nenhum
  `std::uniform_*_distribution` em `tests/`. A elipse de `risk_parity_reference.hpp` (369–370), que usa sin/cos,
  precisa de sin/cos portáteis para |x| ≤ alguns π (porte de `__kernel_sin/__kernel_cos` com Cody–Waite) ou de uma
  parametrização racional da elipse; escolha do builder.
- MSVC: `#pragma fp_contract(off)` no header sob `_MSC_VER`. Clang: o `-ffp-contract=off` do CMake (a pragma de clang
  não foi verificada).
- **Guarda:** teste de valores-ouro (~30 entradas: `exp`, `log`, `hypot`, `exp10`, os primeiros normais de uma semente,
  uma entrada de `random_covariance(4, seed, 1e8)`), comparados por `std::bit_cast`. Incluir alguns `std::fma` de casos
  difíceis: o `fma` de software do MinGW sem `-mfma` é chamada de biblioteca e precisa estar correto.

### 6.3 Performance (condição do autor)

A/B com `quick.py` antes e depois da troca de `hypot`/`log` em `src/`, na mesma sessão, `gcc-bench`. Aceita se nenhum
caso piora além do `noise` (repetir intercalado ou `compare.py` com ≥ 9 repetições para diferenças abaixo de 5%). Se
piorar: **parar e levar ao autor** com os números (opções a apresentar: otimizar a função própria; manter a libm onde
custa e registrar a exceção ao objetivo 3). A troca do cruzamento do MAX_RETURN também é medida (só `max_return/*`),
com a mesma regra e a mesma parada (decisão do Claude: o solve de face roda só no bracket final, numa face, e deve
custar pouco perto da bissecção que poupa).

**Plataformas** (decisão do Claude): A/B completo em `gcc-bench` no Windows (a máquina de medição do projeto) e o
mesmo A/B do `quick.py` no Linux, em Docker na mesma máquina, contra a glibc, cujo `hypot` é mais rápido que o da
MinGW (o quick suite já passa pelos Givens do QP e pelo Spinu; nenhum benchmark novo). O macOS não é medido (não há máquina; runner compartilhado
não dá tempo estável); fica registrado em `docs/performance/`.

### 6.4 Guarda de bits da biblioteca

`tests/unit/solve/determinism_test.cpp` (nome do builder): um caso por objetivo (n pequeno, Σ e μ literais, não
gerados), pesos e os quatro multiplicadores comparados bit a bit com hex-floats gravados. Roda em todo CI; é o que
prova o objetivo 3.

**Se este teste ou o de valores-ouro (§6.2) divergir numa plataforma por outro motivo que a libm** (ex. `std::fma` de
software da MinGW, contração do Apple Clang), decisão do Claude: corrigir na causa, em `src/` ou no CMake — uma flag
ou pragma por compilador, ou mais uma função em `po::portable` (ex. um `fma` próprio pelo algoritmo de Dekker/Veltkamp,
se o da plataforma estiver errado) — com a causa no commit e o A/B da §6.3 se tocar caminho quente. Proibido no teste:
`#ifdef` de SO, skip, tolerância no lugar de bits. Causa não encontrada: parar e levar ao autor.

### 6.5 Troca das instâncias

1. Trocar as chamadas listadas na pesquisa §5b-4 (`random_covariance.hpp`, `random_problems.hpp`,
   `max_return_fuzz_test.cpp:120`, `qr_test.cpp:67`, `tracking_error_test.cpp:354`, `risk_parity_fuzz_test.cpp:202`,
   `risk_parity_reference.hpp:552`, `spinu_test.cpp:58`, `max_sharpe_test.cpp` 221/254/286, `objective_test.cpp:303`).
2. Rodar tudo em MinGW e MSVC; os testes que dependem das instâncias atuais (`validation.md` fala em sete) e qualquer
   falha nova são **triados um por um**: (a) contagem ou valor medido gravado no teste → remedir e atualizar, com o
   motivo no commit; (b) referência de resolução limitada (ex. scan da elipse) → melhorar a referência; (c) defeito real
   do solver → corrigir se for local e pequeno, senão **parar e levar ao autor**. Nunca afrouxar banda (anti-objetivo).
3. Remedir as contagens, piores casos e tempos de `validation.md` com as instâncias novas.

### 6.6 CI e documentação

- `ci.yml`, `cpp.yml`, `wheels.yml` verdes nas três plataformas (push; `gh run watch`).
- CHANGELOG: a seção `[0.1.0]` (ainda não lançada) ganha as mudanças (`Changed`: bits de solução do QP e do
  `MAX_RETURN`, com o `compare-solutions` resumido; `Fixed`: o caso do `MAX_RETURN`) e perde o item de "Known
  limitations" sobre plataformas. A data da seção fica para o dia da tag (autor).
- `validation.md` (fuzz battery, known residual, other platforms), `ARCHITECTURE.md` se citar `hypot`/libm, e o
  comentário de `kMinimumVolatilityTolerance`/afins que citem as constantes removidas.

### 6.7 Ordem

1. Fixtures das duas Σ (Linux e Windows) do caso pinado; teste que falha com a Σ do Linux.
2. Cruzamento do MAX_RETURN (§6.1) → caso pinado e fixture verdes, pin removido; A/B de `max_return`.
3. Matemática portável em `src/` + testes de valores-ouro; troca de `hypot`/`log` em `src/`; A/B (§6.3);
   `compare-solutions` (banda); teste de soluções-ouro (§6.4).
4. Geradores portáveis (§6.5), triagem, `validation.md` remedido.
5. Push, CI verde nas três plataformas, CHANGELOG e docs.
6. Entrega ao autor: rc.2 (`v0.1.0-rc.2`), conferência e tag final.

## 7. Deixado ao builder

- Local e nome exatos dos headers portáveis e do teste de soluções-ouro.
- sin/cos portáteis ou parametrização racional na referência da elipse.
- Como extrair a Σ do Linux (Docker local ou artefato de CI).
- A estrutura interna do solve de face e do ratio test, dentro do §6.1.
- Divisão em commits (um por mudança de bits, para o `compare-solutions` de cada uma ficar no seu commit).

## 8. O que depende do autor

- Decidir se uma regressão de performance medida (§6.3) é aceitável, caso apareça.
- Decidir o que fazer com uma divergência de bits entre plataformas cuja causa não foi encontrada (§6.4).
- Decidir o que fazer com um defeito real do solver revelado pelas instâncias novas, se não for local e pequeno.
- Empurrar `v0.1.0-rc.2` e, depois do rc verde, `v0.1.0`.
