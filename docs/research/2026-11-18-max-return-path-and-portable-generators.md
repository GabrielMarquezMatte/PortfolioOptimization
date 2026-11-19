# Pesquisa: o cruzamento do MAX_RETURN por tangente e geradores de teste portáveis

> **Sobre este arquivo.** Nota de apoio à decisão no formato de `docs/research/2026-11-16-multiplier-accuracy.md`
> (português, termos técnicos e títulos em inglês). Não altera código, plano nem spec. Data: 2026-11-18.
> Contexto lido: `src/problems/max_return.cpp` (`crossing_fraction`/`on_segment` via `bisect`, `at_adjacent`,
> `kAdjacentEnd` na linha 27), `src/problems/parametric_path.{hpp,cpp}` (`Face`, `Shift`, `same_segment` na linha 627,
> `crossing_fraction` nas linhas 729–740, `on_segment` nas linhas 742–754), `src/qp/goldfarb_idnani.cpp`
> (`initial_state`, linhas 209–221), `src/problems/constraints.cpp` (linhas 235–241: a alavancagem como família
> `0'x - ||x||_1 >= -L`), `docs/guide/validation.md` (linhas 1434–1460, "The known residual" e "Other platforms"),
> `tests/unit/random_covariance.hpp`, `tests/unit/fuzz/random_problems.hpp`, `tests/unit/fuzz/max_return_fuzz_test.cpp`
> (linha 330, `kKnownResidual`) e `CMakeLists.txt` (linhas 63–66, 193–196: `-ffp-contract=off` fora do MSVC).
>
> **Convenção de verificação.** "Verificado" = lido na fonte primária (artigo, preprint do autor, documentação ou
> código de referência), com seção/página/arquivo. "Não verificado" = a fonte não foi lida (paywall, bloqueio
> anti-bot, não encontrada) e a afirmação vem de citação secundária ou de memória. "Derivação própria" = argumento
> escrito aqui, conferível linha a linha. "Medido aqui" = experimento numérico desta nota (seção 0), reproduzível pelos
> scripts indicados.

## Resumo

| # | Questão | Veredito | Recomendação | Impacto |
|---|---|---|---|---|
| 1 | Tangente exata (KKT do conjunto ativo) vs diferença de duas soluções do QP | **A diferença é inutilizável no caso residual; a tangente resolve.** Medido: entre soluções GI a 100–700 ulp de t, `(x_b − x_a)/Δt` erra a direção por 10⁶–10⁷ (relativo); a tangente em null-space erra 1e-15. A causa do ruído não é cond(Σ) do problema reduzido: é o ponto de partida do Goldfarb-Idnani, `t Σ⁻¹ v` (6,9e8 no caso), cujo ε·módulo = 1,5e-7 é exatamente o ruído de vol(Q(t)) que `validation.md` mede. O erro do GI fica **ao longo** do caminho (posição em t errada por 1,8e-7), não fora da face (3e-16). | Calcular dw por null-space na face (QR de N, QR de L'·Z, como qpOASES eq. (22)), com o `v` do shift (exato por Sterbenz). Ancorar num **solve de face em t_a** (não no x do GI nem no intercepto w0 quando houver cancelamento). Não resolver o KKT cheio [Σ N; N' 0] por eliminação: cond(K) ≈ 1e13, dw erra 2e-8..6e-8. | `crossing_fraction`, `on_segment` (`parametric_path.cpp` 729–754); nova rotina "face solve". |
| 2 | Certificar o cruzamento t* com conjunto ativo fixo | **Basta o KKT pontual em t\*** (QP convexo: KKT é suficiente). Com a família da alavancagem, o certificado tem de ser **independente da representação** por vetores de sinal: ζ ≥ 0 e \|g_i\| ≤ 1 nos zeros (Derivação própria). Medido: margem 1,7e6 contra ruído 8e-7 no referencial do shift. | Identificar o segmento pela **face** (linhas explícitas ativas, suporte, sinais, conjunto de zeros), não pelos vetores de sinal (`same_segment`). Tolerâncias: primal = `feasibility_tolerance` do GI; dual = `kNoiseFactor·n·ε·(t\|v\|max + \|Σ\|max\|w\|max)`. Se nenhuma das duas faces certifica: passos de homotopia (ratio test) a partir da face inferior; depois MAX_ITER. | Remove `kAdjacentEnd` (1e-5) e o pin `kKnownResidual`. |
| 3 | Prior art (CLA, Best, qpOASES) | **Todos calculam a direção do KKT do conjunto ativo**, nunca por diferença de duas soluções independentes: CLA eq. (5) e dw/dλ explícito (Niedermayer nota 9); qpOASES eq. (9)/(22) ancorado no ponto corrente. O "convex combination" do CLA (Bailey-LdP §4.2) interpola **pontos de virada exatos do mesmo F**, não soluções ruidosas. Ties: qpOASES §3.8 (tolerâncias no ratio test) e §3.10 (ramping); CLA: perturbar μ (Markowitz 1959, via Niedermayer nota 7). Nenhuma fonte trata t ~ 1e16; qpOASES normaliza τ ∈ [0, 1] e o CLA parte de λ = ∞ (vértice) sem resolver QP em t grande. | Adotar o esquema qpOASES/Best no segmento final: direção do KKT, ratio test com `epsNum/epsDen`, ancoragem no ponto corrente. | Desenho do algoritmo (seção 5a). |
| 4 | exp/log/sin/cos portáveis para os geradores | **Viável com ~150 linhas** portando fdlibm `e_exp.c`/`e_log.c` (licença Sun permissiva). Medido: protótipo dá o mesmo hash de 600 000 resultados em g++ 16.2 (-O0/-O2/-O3 -march=native, `-ffp-contract=off`) e MSVC 14.51 (`/O2 /fp:precise`, com e sem `/arch:AVX2`); `-ffp-contract=fast`, `-ffast-math` e `/fp:contract` mudam os bits. Erro máximo 0,80 ulp. O runtime do MinGW arredonda mal 3,3% dos `cos` do gerador (log 0,15%, pow 0,17%). | Header `tests/unit/portable_math.hpp` com `exp`/`log` fdlibm via `std::bit_cast` e literais hex; normais por **Marsaglia polar** (dispensa cos/sin); 10^x = exp(x·ln10) (determinístico, ≤ 35 ulp); teste de valores-ouro em hex. | Todos os testes que usam os geradores (instâncias mudam; já previsto em `validation.md` 1449–1460). |
| 4b | `<random>`: distribuições | **Implementation-defined** ([rand.dist.general]/3). libstdc++, libc++ e MSVC STL usam algoritmos diferentes em `generate_canonical` e `uniform_int_distribution` (verificado no código). Só os engines são especificados (`mt19937_64`: [rand.predef]/4). MinGW (libstdc++) e MSVC já sorteiam instâncias diferentes nos testes que usam distribuições. | Usar só `engine()` cru + mapeamento próprio (como `uniform_open`); inteiros por rejeição própria. | `max_sharpe_test.cpp` (221, 254, 286), `objective_test.cpp` (303). |
| — | Achado fora do escopo | `src/` **usa** libm: `std::hypot` em `make_givens` (`src/qp/active_set.hpp:132`) e `std::log` em `src/riskparity/spinu.cpp:72`. Medido: o `hypot` do MinGW arredonda mal 16,8% das entradas (até 1,09 ulp). O solver não é bit-portável por construção. | Trocar `std::hypot` por `sqrt` com escala por potência de 2 (determinístico, sqrt é correto); avaliar o `log` do Spinu. | Fora deste plano; registrar. |

---

## 0. Experimentos feitos para esta nota (medido aqui)

Scripts no scratchpad da sessão
`C:\Users\x\AppData\Local\Temp\claude\5e03637a68c95476\b7c2ab9c-a650-43cd-b1ef-e78dc77ddf16\scratchpad\` (pode ser
apagado): `probe_path.cpp`, `tangent_exp.py`, `libm_probe.cpp`, `hypot_probe.cpp`, `scripts/cr_check.py`,
`portable_math.hpp`, `portable_driver.cpp`, `scripts/acc_check.py`. Toolchain: g++ 16.2.0 (MSYS2 UCRT64),
MSVC 14.51.36231 (VS 18), `build/gcc-release/libportopt_core.a` no commit `d84ae99`, Python com `numpy` + `mpmath`.

**E1 — o caso residual ao longo do caminho** (`probe_path.cpp`, ligado a `libportopt_core.a`, com
`random_covariance.hpp` do repositório, logo o mesmo Σ dos testes no Windows). Shape seed 3963545680699380072, n = 24,
cond 1e8, escala anual, médias quase iguais (spread 5,02e-7), sem box, cash (rf = 0,02), alavancagem
2,1537482526778033, fração 0,7: cap = 858,42881845239617 (σ_min = 0: tudo em cash). `po::solve` reproduz o pin:
`OPTIMAL`, vol/s − 1 = −1,0059e-7, λ_vol = −1,2888e-13 (`validation.md` 1434–1448: −1,006e-7 e −1,3e-13). Depois a
sonda refaz o `climb` de `max_return.cpp` com `solve_at`:

- O número de linhas geradas ativas sobe até 23 em t ≥ 1e14; bracket [6,4158e15; 1,2832e16].
- Em **toda** a bissecção a face é a mesma: suporte {20, 22}, ambos positivos, 22 zeros, face `Generated` (Kind 3).
  Mesmo assim `same_segment` alterna entre verdadeiro e falso; nas duplas adjacentes 6660585191236981/…982,
  vol/s − 1 = −1,99e-8 / +1,16e-7 e `same_segment` falso.
- Varredura de 800 doubles consecutivos (ulp(t) = 1): vol/s − 1 ∈ [−1,78e-7; +2,03e-7], 390 decréscimos (não
  monótono), **1 suporte e 2 conjuntos de vetores de sinal** alternando.

**E2 — referência exata e variantes em double** (`tangent_exp.py`, 60 dígitos). Na face {w_D = 0, w_20 + w_22 = L}
o caminho é afim e a raiz de vol = s é fechada:

| Quantidade | Valor |
|---|---|
| t\* exato | 6660584005560090,48 (o bracket do GI está em 6660585191236981: erro de 1,8e-7 relativo em t) |
| w\*_S | (0,88915447110145999, 1,26459378157634335); \|w0\| = (1,36; 0,80), t\*\|dw\| = 0,47 (sem cancelamento) |
| ζ (multiplicador da alavancagem) | 2,0536e14; no referencial do shift, ζ_v = ζ − t\*·mid(μ~) = 1,56e8 |
| min (1 − \|g_i\|) nos 22 zeros | 8,43e-9 relativo a ζ (todos os g_i ≈ 0,999999: dual quase degenerado) |
| menor distância entre g_i ordenados | 3,44e-9 |
| início do GI, \|t\* Σ⁻¹ v\|max | 6,9e8; ε × isso = 1,5e-7 |
| hessiana reduzida z'Σz | 1,138e6 |

| Variante (double) | erro em t\* | \|w − w\*\|max | vol/s − 1 |
|---|---|---|---|
| resultado da biblioteca (extremo adjacente) | — | 6,95e-8 | −1,01e-7 |
| x do GI no extremo inferior (sem nada) | — | 1,37e-8 | — |
| tangente, intercepto w0 + t·dw (KKT 3×3 por eliminação, `v` do shift) | 5,3e-16 | 1,6e-13 | −1,4e-17 |
| tangente ancorada no x do GI em t_a | 2,1e-7 | 1,5e-16 | +1,3e-16 |
| tangente ancorada num solve de face em t_a (3×3) | 1,5e-14 | 6,3e-16 | +1,5e-16 |
| pipeline null-space em n = 24 (QR de N com os 23 vetores de sinal, QR de L'Z), ancorado | 1,2e-14 | 4,3e-15 | −2,0e-16 |

| Direção dw (double) | erro relativo |
|---|---|
| KKT 3×3 com μ~ sem shift | 6,1e-7 (cancelamento das médias quase iguais) |
| KKT 3×3 com o `v` do shift | 3,3e-13 |
| KKT cheio 47×47 por eliminação (vetores de sinal / forma geométrica) | 2,3e-8 / 6,0e-8; cond(K) = 7,4e12 / 1,3e13 |
| null-space, n = 24, QR de N (vetores de sinal) | 7,3e-16 (componente fora da face 2,5e-16 relativa) |
| null-space nas coordenadas livres (QR de σ, QR de L'E_F Z_F) | 1,4e-15 |
| diferença de soluções GI a Δt = 700 / 200 / 100 | 8,5e5 / 1,1e7 / 8,2e6 |

Os dois conjuntos de vetores de sinal da varredura têm 23 vetores, concordam com o suporte, e **ambos certificam**
(menor ζ_j = 6,4e4 > 0 nos dois); diferem em 2 vetores cada, todos na coordenada 18. Certificado no referencial do
shift calculado em double no ponto da tangente: ζ_v igual nas duas coordenadas do suporte (diferença 0),
min_D(ζ_v − ρ_i) = 1,73e6 contra o modelo de ruído n·ε·(t\|v\| + \|Σ\|\|w\|) = 8,4e-7; o outro lado da desigualdade
vale 4,1e14.

**E3 — o libm do runtime nos geradores** (`libm_probe.cpp` + `cr_check.py`, 20 000 sorteios do fluxo
`uniform_open` do gerador, comparados com o arredondamento correto em 200 bits; g++ 16.2 UCRT64):

| Função | mal arredondada | erro máximo |
|---|---|---|
| `std::log(u)`, u ∈ (0, 1) | 0,145% | 0,529 ulp |
| `std::cos(2πu)` | 3,255% | 0,790 ulp |
| `std::pow(10, x)`, x ∈ [−8, −5] | 0,170% | 0,510 ulp |
| `std::hypot` (entradas aleatórias, expoentes ±10) | 16,8% | 1,085 ulp |

Constant folding (`constexpr` com `std::exp`/`std::log`, extensão do GCC) vs tempo de execução nos 23 autovalores de
`random_covariance` a cond 1e8: 0 diferenças neste toolchain (não prova ausência em outros).

**E4 — protótipo portável** (`portable_math.hpp`: `exp` e `log` do fdlibm, `std::bit_cast`, literais hex conferidos
contra as palavras hex do fdlibm por `acc_check.py`). Hash FNV de 600 000 resultados (log em (0,1), exp em [−40, 40],
exp(x·ln10) em x ∈ [−10, 0]):

| Compilador e flags | hash |
|---|---|
| g++ -O0 / -O2 / -O3 -march=native, todos com `-ffp-contract=off` | `c420b3634749b00a` |
| MSVC `/O2 /fp:precise` e `/O2 /fp:precise /arch:AVX2` | `c420b3634749b00a` |
| g++ -O3 -march=native `-ffp-contract=fast` (o default do GCC em C++) | `5e7fb66ae8c45e8b` |
| g++ -O3 -march=native `-ffast-math` | `1fd67f0f93c19532` |
| MSVC `/O2 /fp:precise /fp:contract /arch:AVX2` | `e8771a106d557421` |

Precisão (20 000 amostras contra 200 bits): log máx. 0,783 ulp (7,1% não correto), exp máx. 0,797 ulp (9,9%),
exp(x·ln10) contra 10^x máx. 34,2 ulp (o erro de arredondar x·ln10, amplificado por \|x·ln10\| ≤ 23).

---

## 1. Tangente exata vs diferença de duas soluções

### O que as fontes dizem

- **qpOASES** (Ferreau et al. 2014, §2.2, passo 2 e eq. (9), p. 334; Verificado): a direção da homotopia sai do
  sistema KKT do working set, `K_W (Δx, −Δy_W) = (−(g(1) − g(τ)), a_W(1) − a_W(τ))`, e o ponto avança como
  `z(τ⁺) = z(τ) + Δτ Δz` (passo 5): ancorado no ponto corrente, com o lado direito relativo. §2.3.1, eq. (22), p. 337:
  o sistema é resolvido por **null-space** — TQ de A_F, Cholesky RᵀR = ZᵀH_F Z — "a null-space method that can be
  expected to be numerically stable". No código (`src/QProblemB.cpp`, linhas aproximadas, Verificado):
  `determineDataShift` (~1600) faz `delta_g[i] = g_new[i] - g[i]`, `determineStepDirection` (~2700) resolve com o
  fator R e refinamento iterativo opcional (`numRefinementSteps`), `performStep` (~2830) faz `x += tau*delta_x`.
- **CLA, Bailey & López de Prado** (2013, preprint do autor; Verificado): eq. (5), p. 11, pesos livres
  `ω_F = −Σ_F⁻¹Σ_FB ω_B + γ Σ_F⁻¹1_k + λ Σ_F⁻¹μ_F` — afim em λ (o nosso t) com o F do segmento; `computeW`
  (Snippet 9) usa `np.linalg.inv(covarF)`. Entre pontos de virada, "the efficient frontier can be simply derived as a
  convex combination between any two neighbor turning points" (p. 6, nota 7; §4.2 p. 13): interpolação **de pontos de
  virada exatos do mesmo F**, não de soluções independentes ruidosas.
- **Niedermayer & Niedermayer** (DP Bern 07-01, Verificado): o λ de cada candidato vem de fórmula fechada no F corrente
  (eq. (8)–(10), p. 10) e a nota 9 dá a inclinação explícita `dw_i/dλ = −C_i/(1'_F Σ_F⁻¹ 1_F)`. "C_i is zero for all i
  if accidentally μ_F is proportional to 1_F" (p. 10): a degenerescência das médias iguais é a mesma que o shift do
  repositório trata.

Nenhuma das fontes obtém a direção por diferença de duas soluções obtidas independentemente.

### Derivação própria: a tangente na face e o seu erro

Convenção do repositório (`parametric_path.hpp` linhas 50–58): em um segmento, `Σw − t v = N u` e `N'w = b` para as
normais ativas N (colunas), com v o vetor do shift. Daí w(t) = w0 + t·dw, u(t) = u0 + t·du com

```
[Σ  −N] [dw]   [v]        [Σ  −N] [w0]   [0]
[N'  0] [du] = [0]   e    [N'  0] [u0] = [b]
```

Pelo null-space (Z base ortonormal de null(N'), Z'ΣZ = H_r, R'R = H_r via QR de L'Z, com Σ = LL'):
`dw = Z H_r⁻¹ Z'v`, `du = N⁺(Σ dw − v)` (sistema compatível), e w(t) = w_p + Z·H_r⁻¹ Z'(t v − Σ w_p) para uma
solução particular N'w_p = b.

- **Erro da tangente**: δdw/\|dw\| ≈ ε·κ(H_r) mais o erro de Z multiplicado por Σ. Medido 7e-16..1,4e-15 no caso
  residual, com H_r escalar (aresta) e cond(Σ) = 1e8: a tangente não herda cond(Σ), só o condicionamento da face.
- **Erro do ponto**: w(t\*) = w_a + (t\* − t_a)·dw dá δw ≤ δw_a + \|t\* − t_a\|·δdw. O intercepto w0 + t\*·dw é o
  caso t_a = 0: bom quando \|w0\| e t\*\|dw\| são da ordem de \|w(t\*)\| (o caso medido: 1,36 e 0,47), ruim quando o
  segmento é curto em t frente a t (\|w0\| ≫ \|w\|, cancelamento). Ancorar num **solve de face em t_a** (a mesma
  fatoração, lado direito t_a·Z'v − Z'Σw_p) evita os dois problemas e coloca w_D = 0 exatamente.
- **Por que o GI erra**: `initial_state` (`goldfarb_idnani.cpp` 209–221) parte de x = −G⁻¹g = t Σ⁻¹ v e caminha até a
  face; o ponto final carrega ε·\|tΣ⁻¹v\| de erro absoluto — 1,5e-7 no caso, independentemente de quão bem
  condicionada é a face. E2 mostra que esse erro fica quase todo **ao longo do caminho**: o x do GI está na face
  (3e-16) mas é a solução de um t errado em 1,8e-7. Para MAX_RETURN isso basta para estragar vol(Q(t)) = s, e é a
  "rounding of 1e-15 cond(Σ)" de `validation.md` 1442–1444.
- **Diferença de duas soluções**: o erro relativo da direção é ≈ 2ε\|tΣ⁻¹v\| / (\|t_b − t_a\|·\|dw\|). Só presta
  quando o bracket já é largo em w, exatamente o oposto do que a bissecção produz. Medido: 10⁶–10⁷.
- **O shift**: v_i = μ~_i − mid(μ~) é exato (Sterbenz: valores dentro de um fator 2), e a tangente calculada com v
  é 2·10⁶ vezes mais precisa que com μ~ (3,3e-13 vs 6,1e-7 no 3×3). A tangente deve ser calculada no referencial do
  shift, como a hessiana já é.
- **Valores grandes de t**: multiplicar por t é um arredondamento relativo; t = 1e16 não é problema em si. O problema é
  qualquer passo intermediário de módulo ≫ \|w\| (o início do GI). Equivalentemente: reparametrizar com λ = 1/t, como o
  CLA, que parte do vértice (λ = ∞) e nunca forma tΣ⁻¹v.

### As linhas geradas como normais comuns (Derivação própria)

A linha gerada é `(a − s_j)'x ≥ b` com a = 0 (`constraints.cpp` 238–239): a coluna de N é −s_j. Se k vetores de sinal
estão ativos e a família vale (`||w||_1 ≤ L`), então s_j'w = L = ||w||_1 força s_j,i w_i = \|w_i\| para todo i: os
vetores concordam com sign(w) no suporte S, e w_i = 0 em toda coordenada onde dois deles discordam (o conjunto D).
Quando as diferenças s_j − s_1 geram R^D (k − 1 = \|D\|, como o GI garante ao manter normais independentes),

```
range(N) = { (c·σ_S, z_D, 0) : c ∈ R, z ∈ R^D }
```

independe de **quais** k vetores foram escolhidos. Logo dw, w(t) e o segmento dependem só da face (S, σ, D, linhas
explícitas ativas). Os multiplicadores ζ_j **dependem** da escolha; só ζ = Σ_j ζ_j (o multiplicador do usuário) e
g = Σ_j ζ_j s_j / ζ (um subgradiente de ||·||_1) são invariantes. Pela forma geométrica — colunas e_i (i ∈ D) e −σ
em S — a mesma face tem a mesma tangente (medido: 1,4e-15 contra 7,3e-16 pelos vetores de sinal).

### Veredito

A tangente pelo KKT da face, em null-space e no referencial do shift, ancorada num solve de face, é exata a ~1e-15 no
caso que hoje é pinado. A diferença de duas soluções do GI não tem conserto por tolerância: o ruído (1,5e-7) é do
método de partida, não do problema.

---

## 2. Certificar o cruzamento

### Condições (Derivação própria)

Q(t\*) é um QP convexo; as condições KKT num ponto são **suficientes** para otimalidade. Com a face F fixada e (w, u)
dados pela tangente em t\*, estacionaridade e complementaridade valem por construção (linhas ativas com igualdade,
inativas com multiplicador 0). Resta conferir **num ponto só** (nenhum argumento de intervalo é preciso):

1. **Primal.** Linhas explícitas inativas: `C_i w − c_i ≥ −feasibility_tolerance` (a mesma regra do GI,
   `goldfarb_idnani.hpp` 13–17). Família: `family_holds` (`parametric_path.cpp` 433–442). Suporte: σ_i w_i ≥ −tol em
   S (um peso que trocou de sinal muda a face). w_D = 0 exato no solve por coordenadas livres.
2. **Dual, explícitas.** u_i(t\*) ≥ −tol_d nas desigualdades ativas.
3. **Dual, família (invariante à representação).** Com ρ = t\*·v − Σw (menos os termos das linhas explícitas
   ativas), a = `amount` do shift e s_f o vetor de sinal da face do shift: no suporte ρ_i = σ_i·ζ_v (o mesmo ζ_v em
   todo S, e ζ = ζ_v − t\*·a ≥ 0); nos zeros, |g_i| ≤ 1 equivale a

   ```
   −ζ_v + 2 t* a  ≤  s_f,i · ρ_i  ≤  ζ_v ,   i ∈ D
   ```

   Todas as grandezas são "pequenas" (ordem t\|v\| + \|Σ\|\|w\|), exceto 2t\*a, que com médias acima de rf (a < 0)
   torna o lado esquerdo folgado. No caso medido, a margem do lado direito é 1,7e6 contra ruído 8e-7: o certificado
   é robusto mesmo com min(1 − \|g_i\|) = 8,4e-9 relativo a ζ, porque é calculado no referencial do shift.
   Em vez de exigir ζ_j ≥ 0 de um conjunto de vetores de sinal (que muda com o arredondamento — E1), basta isso:
   por Carathéodory, existe um conjunto de \|D\| + 1 vértices do cubo com pesos não negativos exatamente quando
   \|g_i\| ≤ 1, e em geral há vários (os dois do E1 certificam).
4. **Volatilidade.** vol(w(t\*)) − s é só arredondamento de avaliação (medido 1e-16..1e-15): a complementaridade do
   cap volta à faixa ε, sem `kAdjacentEnd`.

**Tolerâncias.** Primal: a do GI (`kRelativeTolerance·max(1, |rhs|, ||r||·|x|max)`), para que o certificado aceite o
que o solver aceitaria. Dual: o modelo de `weight_noise`/`slope_noise` (`parametric_path.cpp` 342–357) estendido,
`kNoiseFactor·n·ε·(t\*·|v|max + |Σ|max·|w|max)` — no caso, 8e-7·16 contra margem 1,7e6. Ambas no referencial do shift.

### Quando nenhuma das duas faces certifica

Ordem proposta (as duas primeiras seguem Best 1996 via qpOASES §2.2, passos 3–7):

1. Candidato da face de t_a (o extremo inferior); se certifica, aceitar.
2. Candidato da face de t_b; se certifica, aceitar.
3. Caso contrário há breakpoint(s) entre as faces. A partir de t_a, ratio test (qpOASES eq. (10)–(12), com as
   tolerâncias da eq. (30): numeradores negativos por arredondamento cortados a 0, denominadores menores que εden
   ignorados) dá o próximo breakpoint t_k; se vol(w(t_k)) ≥ s a raiz está em [t_a, t_k] na face de t_a (o passo 1
   teria certificado salvo empate); senão trocar a face (adicionar a linha bloqueante ou retirar o multiplicador que
   zerou) e repetir. Em empate (dois bloqueios no mesmo t_k), tratar um de cada vez: o seguinte vem com passo 0. A
   continuidade do caminho (H = Σ definida positiva) garante que w no breakpoint pertence às duas faces.
4. Limite de passos (por exemplo 2n); esgotado, MAX_ITER honesto.

Com a identidade de segmento pela face, o caso medido nunca chega ao passo 3: as duas extremidades estão na mesma face.

---

## 3. Prior art

| Método | Direção | Ties / degenerescência | Condicionamento | t grande | Status |
|---|---|---|---|---|---|
| Markowitz (1956) CLA | do sistema do F corrente | Markowitz (1959) propõe perturbar μ ou o método de 1956 | — | parte de λ = ∞ | Não verificado (citado em Niedermayer, p. 8 nota 7) |
| Niedermayer & Niedermayer (2007) | fórmula fechada; dw/dλ explícito (p. 10 nota 9); atualização da inversa por bordas (Lemma 1, p. 13) | "there will hardly be two or more turning points on one (σ, µ) location" (p. 8 nota 7); imprecisão numérica λ(i) = λ_current − ε resolvida pelo sinal de dw/dλ (p. 11 nota 10); C_i ≡ 0 com médias iguais (p. 10) | inversas explícitas | parte do ativo de maior μ (p. 8, §3.1) | Verificado |
| Bailey & López de Prado (2013) | eq. (5) com `np.linalg.inv` | escolha pelo maior λ com comparações estritas (Snippets 3, 7, 8); "the same portfolio is the result of two different F sets" (p. 18) | nenhuma discussão | parte do vértice de maior μ (Snippet 2) | Verificado (preprint) |
| Best (1996), "PQP algorithm" | do KKT do working set | troca de restrição por teste de dependência (passos 6–7 de qpOASES §2.2) | — | — | Não verificado (paywall; conteúdo via qpOASES §2.2, que o cita como [7]) |
| Best & Kale (2000); Ritter (1962) | — | — | — | — | Não verificado (não lidos) |
| Ferreau, Bock, Diehl (2008) | homotopia online | — | — | — | Não verificado (não lido; a versão 2014 resume) |
| qpOASES (Ferreau et al. 2014) | eq. (9), null-space eq. (22), ancorado no ponto (passo 5) | ratio test com `epsNum = −1e3ε`, `epsDen = 1e3ε` (§3.8, eq. (30)); ties evitados por **ramping** (§3.10, eq. (33)–(35)), ativado quando τ ≤ EPS (`performRamping`); flipping bounds só para H indefinida/semidefinida (§3.7) | refinamento iterativo (§3.4); **drift correction** (§3.5, eq. (25)–(27)) | τ ∈ [0, 1] normalizado; critério de término relativo no espaço dos dados (§3.6, eq. (28)) | Verificado (artigo + `QProblemB.cpp`) |

Ponto comum: todos mantêm (w, u) como função afim **exata** do parâmetro dentro do segmento e recalculam a direção a
cada troca de conjunto ativo. A drift correction de qpOASES (§3.5) é o análogo do "solve de face em t_a": repõe o
ponto exatamente sobre as linhas ativas antes de seguir.

---

## 4. Funções matemáticas portáveis para os geradores de teste

### 4.1 Por que as instâncias mudam entre plataformas

Os geradores usam `std::log`, `std::cos` (Box–Muller em `standard_normal`), `std::exp`/`std::log` (autovalores) e
`std::pow(10, x)` (spread das médias). O C e o C++ não exigem arredondamento correto dessas funções, e as
implementações discordam (medido: 3,3% dos `cos` do gerador mal arredondados no MinGW). Derivação própria: uma matriz
24×24 de normais faz 576 `cos`; ≈ 19 diferem entre libm's, e a QR de Householder em `orthogonal_factor` espalha cada
diferença por toda Q — coerente com "561 of 576 covariance entries differ" (`validation.md` 1449–1452).

### 4.2 Bibliotecas candidatas

- **CORE-MATH** (Verificado em core-math.gitlabpages.inria.fr): licença MIT; "correct rounding is claimed for all
  IEEE-754 rounding modes"; testado "on x86_64-linux, with and without the use of fma"; "they do not necessarily
  correctly deal with NaN, infinities or zero"; funções precisam de cabeçalhos auxiliares (log usa `dint.h`).
  Arredondamento correto dá os mesmos bits em qualquer lugar, mas: uso de `__int128`/builtins do GCC e suporte a MSVC
  **Não verificado** (repositório e artigo bloqueados por anti-bot); o tamanho de cada função é bem maior que o do
  fdlibm.
- **fdlibm 5.3** (Verificado, netlib `e_exp.c`, `e_log.c`, `k_cos.c`, `e_rem_pio2.c`): licença Sun "Permission to use,
  copy, modify, and distribute this software is freely granted, provided that this notice is preserved". Só
  aritmética double e manipulação de palavras; nenhum `long double` nem fma. exp: redução x = k·ln2 + r com
  ln2HI de 32 bits ("t*ln2HI is exact here"), aproximação racional com P1..P5, "error is always less than 1 ulp".
  log: x = 2^k(1+f), s = f/(2+f), R com Lg1..Lg7 (Remez, erro 2^-58,45), "less than 1 ulp"; as constantes
  hexadecimais são "the intended ones". Armadilha: `__HI(x) *(1+(int*)&x)` em `fdlibm.h` é type punning por ponteiro
  (UB em C++) e depende de `__LITTLE_ENDIAN` detectado por macros antigas; trocar por `std::bit_cast`.
- **musl** (Verificado, `src/math/exp.c`, `log.c`): código da ARM (MIT, 2018), por tabela; `log.c` linha 89 tem
  `#if __FP_FAST_FMA` com `__builtin_fma` (linha 91) e outro caminho sem fma (linha 94): **os bits dependem de o alvo
  ter fma rápido**. `exp.c` comenta "Without fma the worst case error is 0.25/N ulp larger" (linha 125): foi escrito
  para tolerar contração, logo seus bits variam com ela. Não serve para bit-identidade entre x86-64 e arm64.
- **crlibm** (Verificado, `crlibm.h` no espelho `taschini/crlibm`): LGPL v2 ou posterior; `crlibm_init` "sets FPU
  flags when needed (mostly on Intel architectures with default double extended)". Licença incompatível com copiar
  para um projeto MIT sem cuidado; dependência pesada.

### 4.3 Armadilhas

- **Contração (FMA).** GCC: "The default is -ffp-contract=off for C in a standards compliant mode (-std=c11 or
  similar), -ffp-contract=fast otherwise" (manual, Optimize Options; Verificado) — em C++ o default é `fast`. Clang:
  `-ffp-contract=on` por padrão ("enable C and C++ standard compliant fusion in the same statement"), e
  `-ffp-model=precise` (o default) mantém a contração (UsersManual; Verificado): relevante para Apple clang em arm64.
  MSVC: "Floating-point contractions aren't generated by default under /fp:precise. This behavior is new in Visual
  Studio 2022. Previous compiler versions could generate contractions by default" (doc `/fp`, msvc-170; Verificado);
  `/fp:contract` as reativa. A doc não distingue x64 de ARM64 (arm64 com MSVC: Não verificado). Medido no E4: tudo isso
  muda os bits; o CMake do repositório já põe `-ffp-contract=off` em `portopt_tests` (linhas 193–196).
- **Precisão excedente (x87).** Irrelevante em x86-64 (SSE2) e arm64; a doc do MSVC diz que intermediários "may be
  performed at machine precision" e arredonda em atribuições, casts, argumentos e retornos (Verificado).
- **`std::ldexp`/`frexp`.** Desnecessários: o fdlibm ajusta o expoente somando k<<20 à palavra alta (exato para
  resultado normal). Com `std::bit_cast` não há chamada de biblioteca.
- **Constant folding.** GCC avalia builtins matemáticos com argumentos constantes em tempo de compilação (com MPFR,
  arredondamento correto) e pode discordar do libm em tempo de execução: Não verificado na documentação; medido 0
  diferenças em 23 casos neste toolchain. Funções próprias eliminam o problema: o folding de +, −, ×, ÷ é IEEE exato
  (E4: -O0 = -O3).
- **sqrt.** Usado no Box–Muller/polar; IEEE 754 exige arredondamento correto e `sqrtsd`/`FSQRT` o entregam (Não
  verificado aqui contra o texto da norma; o solver já depende disso).

### 4.4 As distribuições de `<random>`

- "The algorithms for producing each of the specified distributions are implementation-defined."
  ([rand.dist.general]/3; Verificado). Os engines são especificados: "The 10000th consecutive invocation of a
  default-constructed object of type mt19937_64 produces the value 9981545732273789042" ([rand.predef]/4; Verificado).
  O texto atual de [rand.util.canonical] (rejeição até S < x·r^d, resultado em [0, 1)) é o algoritmo novo; libstdc++
  o põe no inline namespace `_V2` salvo `_GLIBCXX_USE_OLD_GENERATE_CANONICAL` (`bits/random.h`; Verificado).
- `uniform_real_distribution`: os três fazem `canonical·(b − a) + a` (libstdc++ via `__detail::_Adaptor`; MSVC
  `_Nrand_impl<_Ty>(_Eng) * (_Max − _Min) + _Min`), mas `generate_canonical` difere: libc++ soma em ponto flutuante
  sem rejeição e sem tratar o resultado 1 (`__random/generate_canonical.h` linhas 41–44); MSVC acumula em inteiro com
  laço de rejeição; libstdc++ `_V2`. (Verificado.)
- `uniform_int_distribution`: libstdc++ usa Lemire "nearly divisionless" (`_S_nd`, `bits/uniform_int_dist.h`) com
  fallbacks por divisão e upscaling; libc++ usa `__independent_bits_engine` + rejeição (linhas 210–230); MSVC usa
  `_Rng_from_urng_v2` (algoritmo interno não lido). (Verificado, salvo o interior do MSVC.)

Consequência: `max_sharpe_test.cpp` (linhas 221, 254, 286) e `riskparity/objective_test.cpp` (303) já sorteiam
instâncias diferentes entre MinGW (libstdc++) e MSVC no mesmo Windows.

### 4.5 Desenho recomendado

Ver seção 5b.

---

## 5. Recomendação para o plano

### (a) O cruzamento do MAX_RETURN por tangente, com certificado

1. **Identidade do segmento pela face.** `same_segment` passa a comparar: linhas explícitas ativas, e para a família,
   o suporte com sinais e o conjunto de zeros D (coordenadas onde os vetores de sinal discordam), mais o fato de a
   família estar ativa. Os vetores de sinal deixam de entrar na comparação (E1: a mesma face com 2 representações).
2. **Face solve** (nova rotina, ao lado de `multiplier_slope`): dadas a face e o shift de um ponto, montar N nas
   coordenadas livres F = complemento de D (zeros eliminados exatamente; linhas da família como −σ_S; colunas ±1 do
   budget/cash exatas), QR de Householder de N_F (já existe `linalg::householder_qr`), Z_F, QR de L'·E_F·Z_F (fator
   de Data), e resolver para (i) dw = Z H_r⁻¹ Z'v e (ii) w(t) num t dado. Multiplicadores por `least_squares` em
   `N u = Σw − t v` no referencial do shift, como `multiplier_slope` faz para b.
3. **Cruzamento.** Com o bracket [t_a, t_b] da bissecção atual (os pontos GI continuam decidindo o bracket):
   w_a = face solve em t_a (face de t_a); h pela fórmula sem cancelamento de `crossing_fraction` com start = L'w_a e
   direction = L'dw; t\* = t_a + h; w\* = w_a + h·dw. Não é preciso esperar t_a, t_b adjacentes: pode-se tentar desde
   o primeiro bracket e bissectar só se o certificado falhar.
4. **Certificado** (seção 2): primal com a tolerância do GI, dual explícito e da família na forma invariante, no
   referencial do shift, com `kNoiseFactor·n·ε·(t|v|max + |Σ|max|w|max)`. Falhou: a face de t_b; falhou: passos de
   homotopia com ratio test (qpOASES eq. (30)); esgotado: MAX_ITER.
5. **Resultado.** `at_crossing` com t\* e (w\*, u\*): multiplicadores do usuário u/t\*, λ_vol = −s/t\*, como hoje.
6. **Remoções.** `kAdjacentEnd`, `kCrossingSanity` (o certificado o substitui) e o pin `kKnownResidual`. Critério de
   aceitação: o caso 3963545680699380072/0,7 com vol/s − 1 na faixa ε e retorno igual à referência; e o mesmo caso
   com o Σ do Linux (perturbação de 3,6e-13) `OPTIMAL` nas frações 0,3–0,9 (hoje MAX_ITER, `validation.md` 1452–1455)
   — não medido aqui, mas o mecanismo (a tangente não depende do ruído do GI) prevê que sim.

Riscos: faces com \|D\| > k − 1 (as diferenças dos vetores não geram R^D) — a face geométrica é então menor que a do
GI; tratar como "certificado falhou" e cair no passo 4. A tangente só vale se H_r for bem condicionada; quando
`R` tiver pivô pequeno (o análogo do `epsFlipping` do qpOASES), registrar e cair na bissecção.

### (b) O header de matemática portável

1. `tests/unit/portable_math.hpp`, namespace `po::testing::portable`:
   - `exp(double)` e `log(double)`: porte direto de fdlibm 5.3 (`e_exp.c`, `e_log.c`), mantendo o aviso de copyright
     da Sun no topo (exigência da licença; compatível com o MIT do repositório). Acesso às palavras por
     `std::bit_cast<std::uint64_t>`; constantes como literais hex (conferidas no E4); casos especiais reduzidos ao
     domínio dos geradores (finito; x > 0 no log), com `assert` fora dele. ~150 linhas.
   - `exp10(double x)` = `exp(x * 0x1.26bb1bbb55516p+1)`: determinístico; erro ≤ ~35 ulp para \|x\| ≤ 10 (medido),
     irrelevante para um spread log-uniforme. Se for preciso mais, dividir x·ln10 em hi + lo por Dekker (só +, ×).
   - **Sem sin/cos**: `standard_normal` passa a ser o método polar de Marsaglia (sortear (u1, u2) em (−1, 1)², rejeitar
     r² ≥ 1 ou r² = 0, multiplicar por sqrt(−2 log r² / r²)). Usa só `log` e `sqrt`; o número de chamadas ao engine
     varia, mas é determinístico. (A referência do método, Marsaglia & Bray 1964, não foi lida: Não verificado; o
     algoritmo é conferível linha a linha.) Para a elipse de `risk_parity_reference.hpp` (369–370), se for preciso,
     portar `__kernel_sin/__kernel_cos` com a redução de Cody–Waite de `e_rem_pio2.c` (pio2_1 de 33 bits + pio2_1t)
     para \|x\| ≤ alguns π.
   - Inteiros e reais uniformes só a partir de `engine()` cru: `uniform_open` (já existe) e um `uniform_index(n)`
     por rejeição própria; nada de `std::uniform_*_distribution`.
2. **Compilação.** Manter `-ffp-contract=off` em `portopt_tests` (já está); MSVC `/fp:precise` (default desde VS 2022,
   sem contração); acrescentar `#pragma fp_contract(off)` sob `_MSC_VER` no header (a pragma consta da doc `/fp`).
   Para Clang, `#pragma clang fp contract(off)` (Não verificado na doc lida). Proibir `-ffast-math` nesses alvos.
3. **Guarda.** Um teste de valores-ouro: ~30 entradas fixas com os bits esperados de `exp`, `log`, `exp10`,
   `standard_normal` (primeiros 8 de uma semente) e uma entrada de `random_covariance(4, seed, 1e8)`, comparados com
   `std::bit_cast`. Falha em qualquer plataforma em que a contração ou o folding mudem algo.
4. **Migração.** Trocar `std::log/std::exp/std::cos/std::pow` em `random_covariance.hpp`, `random_problems.hpp`,
   `max_return_fuzz_test.cpp:120`, `qr_test.cpp:67`, `tracking_error_test.cpp:354`, `risk_parity_fuzz_test.cpp:202`,
   `risk_parity_reference.hpp:552`, `spinu_test.cpp:58`; `std::log10` só em nomes de propriedades
   (`static_cast<int>(std::log10(1e8))`) é seguro na prática. Todas as instâncias mudam: refazer as contagens de
   `validation.md` e os sete testes que dependem das instâncias atuais (`validation.md` 1458–1460).

---

## 6. Fora do escopo, mas registrado

- `src/qp/active_set.hpp:132` usa `std::hypot` em `make_givens`; medido 16,8% de resultados mal arredondados no
  MinGW. O C não exige arredondamento correto de `hypot` (Não verificado no texto da norma aqui). Uma versão própria
  (escala por potência de 2 escolhida pelo expoente, `sqrt(a² + b²)`, desescala) é determinística porque sqrt é
  correto. `src/riskparity/spinu.cpp:72` usa `std::log` no valor do objetivo; se o valor decide a busca linear, as
  iterações podem diferir entre plataformas.
- A frase de `validation.md` "the solver is the same, the problem is not" vale para o status observado, não para os
  bits: com `hypot` de libm's diferentes, os iterados do GI podem diferir no último bit.

---

## Bibliografia

- Bailey, D. H.; López de Prado, M. (2013). An Open-Source Implementation of the Critical-Line Algorithm for Portfolio
  Optimization. *Algorithms* 6(1), 169–196. Preprint do autor, https://www.davidhbailey.com/dhbpapers/critical-line.pdf
  — **Verificado** (pp. 5–18 do preprint).
- Niedermayer, A.; Niedermayer, D. (2007). Applying Markowitz's Critical Line Algorithm. Discussion Paper 07-01,
  Universität Bern, https://repec.vwiit.ch/dp/dp0701.pdf — **Verificado** (pp. 7–13).
- Ferreau, H. J.; Kirches, C.; Potschka, A.; Bock, H. G.; Diehl, M. (2014). qpOASES: a parametric active-set algorithm
  for quadratic programming. *Math. Prog. Comp.* 6(4), 327–363. Versão do autor,
  https://cdn.syscop.de/publications/Ferreau2014.pdf — **Verificado** (§2.1–2.3, §3.1–3.10, tabela de opções).
- qpOASES, `src/QProblemB.cpp`, https://github.com/coin-or/qpOASES — **Verificado** (nomes de funções e trechos;
  números de linha aproximados).
- Best, M. J. (1996). An Algorithm for the Solution of the Parametric Quadratic Programming Problem. In *Applied
  Mathematics and Parallel Computing*, Physica-Verlag, doi:10.1007/978-3-642-99789-1_5 — **Não verificado** (paywall).
- Best, M. J.; Kale, J. (2000), Ritter, K. (1962), Markowitz, H. (1956, 1959), Ferreau, Bock, Diehl (2008) — **Não
  verificado** (não lidos).
- Goldfarb, D.; Idnani, A. (1983). A numerically stable dual method for solving strictly convex quadratic programs.
  *Math. Prog.* 27, 1–33 — **Não verificado** (Springer bloqueado); o ponto de partida irrestrito foi verificado no
  código do repositório (`goldfarb_idnani.cpp` 209–221).
- fdlibm 5.3, https://netlib.org/fdlibm/ (`e_exp.c`, `e_log.c`, `k_cos.c`, `e_rem_pio2.c`, `fdlibm.h`) — **Verificado**.
- musl, `src/math/exp.c`, `src/math/log.c`, https://git.musl-libc.org/cgit/musl/ — **Verificado**.
- CORE-MATH, https://core-math.gitlabpages.inria.fr/ — **Verificado** (página do projeto); repositório e artigo
  (Sibidanov, Zimmermann, Glondu, ARITH 2022) — **Não verificado** (bloqueio anti-bot).
- CRlibm, `crlibm.h`, https://github.com/taschini/crlibm — **Verificado** (licença e `crlibm_init`).
- Microsoft, `/fp (Specify floating-point behavior)`, https://learn.microsoft.com/en-us/cpp/build/reference/fp-specify-floating-point-behavior — **Verificado**.
- GCC, Optimize Options (`-ffp-contract`), https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html — **Verificado**.
- Clang User's Manual (`-ffp-contract`, `-ffp-model`), https://clang.llvm.org/docs/UsersManual.html — **Verificado**.
- C++ working draft: [rand.dist.general], [rand.predef], [rand.util.canonical], https://eel.is/c++draft/ — **Verificado**.
- libstdc++ `bits/random.h`, `bits/uniform_int_dist.h`; libc++ `__random/generate_canonical.h`,
  `__random/uniform_int_distribution.h`; MSVC STL `stl/inc/random` — **Verificado** (via GitHub, ramos principais em
  2026-11-18).
- Lemire, D. (2019). Fast Random Integer Generation in an Interval, arXiv:1805.10941 — **Não verificado** (citado no
  código do libstdc++).
- Marsaglia, G.; Bray, T. A. (1964). A convenient method for generating normal variables. *SIAM Review* 6(3) —
  **Não verificado**.
