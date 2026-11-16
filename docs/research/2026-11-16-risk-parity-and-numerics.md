# Pesquisa: RISK_PARITY (plano 4) e decisões numéricas dos planos 1–3

> **Sobre este arquivo.** Não existia convenção no repositório para notas de pesquisa (não havia `docs/research/`);
> este arquivo cria a pasta e segue o formato de `docs/superpowers/` (português, termos técnicos e títulos em inglês).
> É um documento de apoio à decisão, não um plano nem uma spec: nada aqui altera o plano 4 até que o autor decida.
> Data: 2026-11-16. Contexto lido: `docs/superpowers/plans/2026-11-16-portopt-04-risk-parity.md` (citado como
> "plano", com número de linha), `docs/superpowers/specs/2026-11-15-portopt-design.md` ("spec"),
> `docs/guide/math.md`, `docs/guide/validation.md` e o código em `src/`.
>
> **Convenção de verificação.** "Verificado" significa lido na fonte primária (artigo, preprint do autor, livro,
> código de referência), com página/seção indicada. "Não verificado" significa que a fonte não pôde ser lida
> (paywall, TLS, SSRN bloqueado) e a afirmação vem de citação secundária ou de memória, marcada como tal.
> "Derivação própria" é um argumento matemático escrito aqui, conferível linha a linha, sem fonte que o enuncie
> para este problema.

## Resumo

| # | Item | Veredito | Recomendação | Impacto |
|---|---|---|---|---|
| 1 | Lagrangiano aumentado (AL) com bounds no subproblema; L-BFGS projetado | **Confirma** a arquitetura (AL com bounds no nível inferior é o LANCELOT/ALGENCAN). **Contradiz** a atribuição da regra de ρ ("×10 se ‖h‖ não cai 4×, N&W §17.4": é a regra do ALGENCAN, não a do N&W). O "L-BFGS projetado nas livres" tem armadilhas documentadas (Bertsekas 1982). | Pôr igualdades lineares + box no nível inferior e resolver o subproblema por passos Gauss-Newton com o QP de Goldfarb-Idnani (já validado), deixando no AL só vol(w) = s; tolerância ε_k → 0; multiplicadores salvaguardados; polimento final Newton-KKT com Hessiana exata. Se o L-BFGS projetado ficar: conjunto ε-ativo de Bertsekas, Armijo ao longo do arco, descarte de pares com s'y ≤ 0. | Task 3 (`box_lbfgs` pode sair inteira) e Task 4 (`augmented_lagrangian`, política linhas 300–307). Nenhum código escrito afetado. |
| 2 | SCRIP (Feng-Palomar 2015), Bai-Scheinberg-Tütüncü 2016 | **Confirma** que sequência de QPs convexos é o método de referência para restrições lineares/box (código do `riskParityPortfolio` usa `quadprog`, isto é, Goldfarb-Idnani). **Nenhuma** das fontes trata vol-alvo como igualdade. Teorema de convergência do SCRIP **não verificado** (PDF inacessível). | Reusar `qp::solve` como solver do subproblema (item 1). Vol-alvo continua no AL. | Task 3/4. |
| 3 | Existência/unicidade por ortante | **Contradiz** a spec §3.3 ("uma solução por ortante"): há no máximo 2^(n−1) soluções normalizadas (BST 2013, Lema 2.2); em cada par de ortantes opostos há uma só, e nenhuma quando a soma dos pesos do raio é 0. **Confirma**: no ortante positivo existe e é única para Σ definida positiva. | Corrigir a frase da spec/`math.md`; documentar que `OPTIMAL` com bounds pode perder uma solução exata (f = 0) de outro ortante; opcional: varrer ortantes quando poucos ativos podem trocar de sinal. | Spec §3.3, `math.md` (Task 7), Review Focus 3. |
| 4 | Spinu 2013; Griveau-Billion-Richard-Roncalli 2013 | **Confirma** que F não é autoconcordante padrão com bᵢ < 1 e que Armijo + fração até a fronteira é sólido. **Corrige**: F/min bᵢ é autoconcordante padrão (a teoria amortecida existe, mas a cota de iterações escala com 1/min bᵢ). Spinu (via GRR) usa ‖Δy/y‖∞, não o decremento, para o passo. | Critério de parada do plano (λ² ≤ 16nε(1+\|F\|), linha 215) deve ser testado **antes** de dar o passo e o passo dado mesmo assim, ou trocado por ‖Δy/y‖∞ ≤ c·ε, senão o erro relativo do peso é ~√ε/√bᵢ. | Task 2 (`spinu`). |
| 5 | Conjunto viável de um ponto; tolerâncias 1e-12 e n·ε | **Contradiz** o plano em s = σ_min: LICQ falha e, genericamente, **não existem** multiplicadores de KKT (não é questão de "mínima norma"); a derivada ∂f*/∂s é −∞. O caso n = 2 com `target_return` **não** é degenerado (LICQ vale). **Confirma** n·ε para Σbᵢ (cota de Rump/Higham + arredondamento de bᵢ, com folga 2). 1e-12 em σ_min é defensável pela tolerância de viabilidade do QP, não por arredondamento de σ. | Em s = σ_min reportar λ_vol = −∞ (ou NaN) e os demais NaN quando a estacionaridade não fecha, como o MAX_RETURN já faz (`math.md` 208–216); separar os dois casos do Review Focus 1; amarrar `kMinimumVolatilityTolerance` ao `kRelativeTolerance` do QP. | Task 4 (testes `VolatilityEqualMinimumIsMinimumVariance`, `TwoAssetsWithTargetReturnIsThePoint`), `expect_risk_parity_kkt` (Task 0, já em andamento). |
| 6a | Leverage como família implícita no Goldfarb-Idnani | **Confirma**: é o algoritmo de Tibshirani (1996, §6) para o lasso, com terminação finita na versão "só acrescenta". O argumento de terminação do GI com remoções **não verificado** no texto de 1983. | Citar Tibshirani 1996 §6 em `math.md`; nenhuma mudança de código. | `src/qp/goldfarb_idnani.cpp` 209–270 (já escrito): sem mudança. |
| 6b | Guarda de ruído do κ (MAX_SHARPE) | **Confirma com ressalvas**: o modelo é uma cota de primeira ordem no espírito de Higham; o erro retroativo rigoroso do Cholesky usa \|Rᵀ\|\|R\|, não \|Σ\| (não verificado no livro). Britten-Jones 1999: o erro amostral dos pesos tangentes é ordens acima do arredondamento. | Documentar a diferença \|Σ\| × \|Rᵀ\|\|R\| e a margem medida (0,040 κ_noise); nota para o usuário sobre erro amostral. | `src/problems/homogeneous.cpp` 21–24, 204, 316–320: sem mudança. |
| 6c | MAX_RETURN por caminho paramétrico × CLA | **Confirma**: caminho afim por partes e interpolação no segmento são a Proposição 1 e o Corolário 1 de Niedermayer & Niedermayer. Simplificação: dobra+bissecção em vez do teste de razão exato do CLA. | Citar N&N em `math.md`; opcional: teste de razão exato para o próximo breakpoint. | `src/problems/parametric_path.cpp`, `max_return.cpp`: sem mudança. |
| 6d | Homogeneização com leverage | **Confirma** (derivação própria, abaixo): a linha L·1'y − ‖y‖₁ ≥ 0 é equivalente a ‖w‖₁ ≤ L, implica 1'y ≥ ‖y‖₁/L ≥ 1/(L‖v‖∞) > 0 e torna `UNBOUNDED` impossível. | Pôr a cota κ ≥ 1/(L‖v‖∞) em `math.md`/`validation.md` como prova do `cone_bounded`. | `homogeneous.cpp` 279 (`cone_bounded`): confirma. |

---

## 1. Lagrangiano aumentado, L-BFGS-B, L-BFGS projetado e Newton com conjunto ativo

### O que as fontes dizem

- **Convenção e atualização.** Nocedal & Wright (2006), §17.3, eq. (17.36): L_A(x, λ; μ) = f(x) − Σλᵢcᵢ(x) +
  (μ/2)Σcᵢ²(x); eq. (17.39): λᵢ ← λᵢ − μcᵢ(x_k); Framework 17.3 (p. 515). **Verificado** (pp. 514–515).
  - O plano usa exatamente essa convenção (linha 52: L = f − λ'h + (ρ/2)‖h‖², λ ← λ − ρh). Com ∇f = Σλₖ∇hₖ + Σβᵢeᵢ,
    o λ do plano é ∂f*/∂(lado direito), coerente com a regra global (linha 29). O ALGENCAN usa o sinal oposto
    (Andreani-Birgin-Martínez-Schuverdt 2007, eq. (2.2): f + (ρ/2)Σ(hᵢ + λᵢ/ρ)², eq. (3.5): λ ← λ̄ + ρh), isto é,
    λ_ALGENCAN = −λ_plano. **Confirma**; só é preciso não misturar as duas ao copiar fórmulas.
- **Formulação com bounds no subproblema.** N&W §17.4, "Bound-constrained formulation", eqs. (17.48)–(17.52) e
  **Algorithm 17.4** (Bound-Constrained Lagrangian Method, o do LANCELOT), pp. 519–521: os bounds ficam no subproblema,
  resolvido por projeção de gradiente; o teste de parada do subproblema é ‖x − P(x − ∇L_A)‖ ≤ ω_k. **Verificado.**
  - A regra de ρ do Algorithm 17.4 é: se ‖c(x_k)‖ ≤ η_k, atualiza λ e mantém μ, com η_{k+1} = η_k/μ^0.9 e
    ω_{k+1} = ω_k/μ; senão μ_{k+1} = **100**·μ_k, η = 1/μ^0.1, ω = 1/μ (p. 521). As constantes "são até certo ponto
    arbitrárias" (p. 521).
  - A regra "manter ρ se a inviabilidade caiu por um fator τ, senão ρ ← γρ" é a do **ALGENCAN**: Andreani, Birgin,
    Martínez, Schuverdt (2007), Algorithm 3.1, Step 4 (p. 6), com **τ = 0,5, γ = 10** nos experimentos (§6, p. 14),
    ρ₁ = max{1e-6, min{10, 2|f(x₀)|/(‖h(x₀)‖² + ‖g(x₀)₊‖²)}}, e multiplicadores **salvaguardados** (projetados numa caixa
    [λ_min, λ_max], eq. (3.6); 1e±20 nos experimentos). **Verificado.**
- **Garantias sem qualificação de restrição.** ABMS 2007, Teorema 4.1 (p. 6): se ρ_k fica limitado, o ponto limite é
  viável; senão é estacionário da inviabilidade ou viola CPLD. Os autores dizem que multiplicadores muito grandes
  "anunciam" pontos limite sem CPLD (p. 3) e que métodos AL **não têm convergência local rápida**: com tolerâncias
  mais estritas que 1e-4 o desempenho deteriora em relação ao IPOPT (p. 17). **Verificado.**
  - Condições sequenciais (AKKT/CAKKT): Andreani, Martínez, Svaiter (2010) provam que esse AL gera sequências cujos
    limites satisfazem uma condição sequencial nova, válida em minimizadores locais **sem** qualificação de restrição.
    Título e autores verificados (registro FAPESP); o periódico (SIAM J. Optim. 20(6), 2010) e a identificação exata
    com "CAKKT" **não verificados**.
- **L-BFGS-B.** Byrd, Lu, Nocedal, Zhu (1995): projeção de gradiente + matriz L-BFGS (abstract verificado). Os
  detalhes (ponto de Cauchy generalizado ao longo do caminho projetado, depois minimização no subespaço das variáveis
  livres nesse ponto) são **de memória, não verificados** no texto. A versão 3.0 (Morales & Nocedal 2011, "Remark on
  Algorithm 778") é a usada pelo SciPy, cujo teste de parada é max|proj gᵢ| ≤ pgtol (doc. do SciPy 1.18, verificado),
  o mesmo critério do `BoxResult::converged` do plano (linha 255).
- **Projeção ingênua de um passo quase-Newton falha.** Bertsekas (1982), §2, Fig. 1 (p. 225): com D_k definida
  positiva arbitrária, x(α) = [x − αD∇f]⁺ pode não descer para nenhum α. Proposição 1 (p. 226): há descida se D é
  **diagonal com relação a I⁺(x)** = {i : xᵢ no bound e a derivada empurra para fora}. Como I⁺ é descontínuo na
  fronteira ("zigzagging or jamming", p. 228), o algoritmo usa o conjunto **ε-ativo** I_k⁺ = {i : 0 ≤ xᵢ ≤ ε_k,
  ∂f/∂xᵢ > 0}, eq. (32), e Armijo **ao longo do arco projetado**, eqs. (7b), (33)–(35). **Verificado.** Kim, Sra,
  Dhillon (2010) estendem isso a BFGS/L-BFGS com um conjunto fixo em duas etapas (título, periódico e abstract
  verificados; definição exata do conjunto **não verificada**).

### Onde confirma ou contradiz o plano

- Confirma: bounds no subproblema (linha 10), atualização λ ← λ − ρh (linha 52), critério ‖x − P(x − g)‖∞ (linha 255).
- Contradiz a atribuição: a linha 304 diz "aumentam ρ ×10 quando ‖h‖ não cai por um fator 4 (Nocedal-Wright §17.4)".
  O N&W §17.4 usa ×100 com limiares η_k; a regra "fator τ, ×γ" é do ABMS 2007 (τ = 0,5, γ = 10). Um τ = 0,25 é uma
  variante legítima, mas sem fonte.
- Lacunas no plano em relação às fontes:
  1. **Tolerância do subproblema fixa.** `BoxOptions::tolerance` (linha 252) é uma constante; a teoria (N&W Alg. 17.4,
     ω_k → 0; ABMS, ε_k → 0) exige uma sequência decrescente.
  2. **Sem salvaguarda de multiplicadores.** O ABMS projeta λ numa caixa para manter a propriedade de penalidade
     externa; sem ela, um λ que diverge (o caso s ≈ σ_min do item 5) contamina os subproblemas.
  3. **"L-BFGS nas variáveis livres" não diz como as livres são escolhidas.** Pela Proposição 1 de Bertsekas, se a
     escolha for "xᵢ exatamente no bound", a direção pode não descer; precisa ser o conjunto ε-ativo e Armijo ao longo
     do arco. Com f não convexa, o par (s, y) do L-BFGS pode ter s'y ≤ 0 com uma busca só de Armijo (sem Wolfe), e o
     par tem de ser descartado (prática padrão; referência exata **não verificada**).
  4. **Precisão final.** O ABMS relata que o AL deteriora com tolerâncias apertadas (p. 17), e o plano exige KKT a
     1e-10 (linhas 39–45). A restauração por Gauss-Newton e os multiplicadores por QR (linha 304) cobrem a
     viabilidade e os multiplicadores, mas não a estacionaridade: se o L-BFGS parou em 1e-8, o ponto não é KKT a 1e-10.

### A Hessiana exata justifica Newton?

f = ‖r(w)‖² com r = c(w) − b. ∇²f = 2JᵀJ + 2Σᵢ rᵢ∇²cᵢ. Na solução de Spinu r = 0 e a Hessiana de Gauss-Newton 2JᵀJ é
exata; com restrições ativas r ≠ 0 e o termo de segunda ordem pesa, a Hessiana pode ser indefinida (o problema é não
convexo: BST 2013, Exemplos 3.1 e 3.2, pp. 11–12, mostram ponto estacionário de máximo e mínimo local espúrio). Além
disso c é homogênea de grau 0, então Jw = 0 e JᵀJ é **singular** (o plano já registra w'∇f = 0, linha 50). Logo:

- Newton puro com Hessiana exata precisa de modificação (indefinida e singular): não serve como método global.
- Gauss-Newton com termo proximal, 2JᵀJ + τI (τ > 0), é definida positiva e é exatamente o QP do SCRIP (item 2).
- A Hessiana exata é útil **no fim**: com o conjunto ativo identificado, um ou dois passos de Newton no sistema KKT
  (Hessiana exata do Lagrangiano, linhas ativas como igualdades) dão convergência local quadrática e levam a
  estacionaridade ao piso do arredondamento, que é o que a faixa 1e-10 pede.

### Recomendação

1. **Trocar o subproblema.** No arcabouço do ABMS ("general lower-level constraints", eq. (1.1): Ω₂ é o conjunto
   "fácil"), pôr no nível inferior **orçamento, `target_return` e o box** (poliedro), e no AL só h₃ = vol(w) − s.
   O subproblema min L(w) = ‖c − b‖² − λh₃ + (ρ/2)h₃² sobre o poliedro é um mínimos quadrados não linear:
   resíduo r̃ = (c − b, √ρ·(h₃ − λ/ρ)). Cada passo é o QP min ½d'(2J̃ᵀJ̃ + τI)d + ∇L'd s.a. linhas lineares e box,
   resolvido por `qp::solve` (Goldfarb-Idnani), seguido de Armijo em L no segmento w + γd (sempre viável, poliedro
   convexo). Vantagens: as restrições lineares valem **exatamente** em todo iterado (some a restauração de
   h₁, h₂); os multiplicadores das linhas lineares e dos bounds saem do QP no ponto fixo (d = 0 ⇒
   ∇L(w) = Nu, que é a estacionaridade original com λ₃ = λ − ρh₃); nenhum módulo novo de L-BFGS para validar.
   Sem `target_volatility`, não há AL: é o SCA puro sobre o poliedro.
2. Regra de ρ: a do ABMS (τ = 0,5, γ = 10, ρ₁ pela fórmula do §6), citada como tal; ε_k → 0 (por exemplo
   ε_{k+1} = max(ε_final, 0,1·ε_k)); λ salvaguardado em [−λ_max, λ_max].
3. Polimento: com o conjunto ativo do último QP, passos de Newton no sistema KKT com a Hessiana exata (fatoração
   do sistema aumentado; o código já tem QR e Cholesky); aceitar só se a viabilidade e o resíduo KKT caem.
4. Se o L-BFGS projetado for mantido (Task 3 como está): conjunto ε-ativo de Bertsekas, Armijo ao longo do arco
   [x − αDg]⁺, descarte de pares com s'y ≤ ε_mach‖s‖‖y‖, e o teste `BoxLbfgs.AgainstQpSolver` (linha 265) já é o
   teste certo; acrescentar um caso degenerado (bound ativo com gradiente nulo, onde I⁺ é descontínuo).

---

## 2. SCRIP (Feng-Palomar 2015) e Bai-Scheinberg-Tütüncü 2016

### O que as fontes dizem

- **Feng & Palomar (2015)**, IEEE TSP 63(19):5285–5300. O PDF do autor (palomar.home.ece.ust.hk) recusou o
  certificado TLS e não foi baixado; o conteúdo do artigo **não foi verificado**. Verificado em fontes do mesmo autor:
  - Palomar (2025), *Portfolio Optimization: Theory and Application*, cap. 11, §11.7 (versão online): formulação
    min Σgᵢ(w)² + λF(w) s.a. w ∈ 𝒲, com 𝒲 contendo restrições como turnover, market-neutral e posição máxima; o
    SCA lineariza gᵢ, resolve um QP com Q^k = 2(J^k)ᵀJ^k + τI ("a regularizer to obtain a strongly convex function as
    required by SCA") e faz w^{k+1} = w^k + γ^k(ŵ − w^k). A seção não enuncia teorema de convergência e remete ao
    artigo. Nenhuma vol-alvo aparece.
  - Código de referência `riskParityPortfolio` 0.2.2 (CRAN; GitHub `dppalomar/riskParityPortfolio`, commit
    `7d8852c`), `R/riskParityPortfolio.R`, `riskParityPortfolioSCA`, linhas 128–181: Q_k = 2A_kᵀA_k + τI com
    τ = 0,05·tr(Σ)/(2N), q_k = 2A_kᵀg(w_k) − Q_k w_k, **`quadprog::solve.QP`** (o Goldfarb-Idnani) com as linhas
    `Cmat w = cvec`, `Dmat w ≤ dvec` e o box; passo w ← w + γ(ŵ − w) com γ ← γ(1 − ζγ) (γ₀ = 0,9, ζ = 1e-7); parada
    por variação relativa de w (wtol = 5e-7) **ou** de f (ftol = 1e-8), não por KKT. A formulação
    `"rc-over-var vs b"` é exatamente a f do plano. Retorno e variância entram **só como penalidade**
    (`lmd_mu`, `lmd_var`), nunca como igualdade. **Verificado.**
- **Bai, Scheinberg, Tütüncü** (preprint Optimization Online, revisão de 21/10/2013; publicado em Quantitative
  Finance 16(3):357–376, 2016). **Verificado no preprint**:
  - Formulação (2.17)/(2.18): min Σ(xᵢ(Σx)ᵢ − θ)² com θ livre, bounds e orçamento (contribuições não normalizadas,
    diferente da f do plano).
  - Lema 3.1 (p. 11): sem restrições, todo ponto estacionário de F(x, θ) é ótimo global com F = 0. Com restrições:
    Exemplo 3.1 (p. 12) tem ponto estacionário não-RP (máximo local) só com o orçamento; Exemplo 3.2 tem mínimo local
    não global causado por bounds, com uma solução RP no mesmo ortante.
  - Algoritmos 2 e 3 (ALM, pp. 16–17): linearização alternada em duas cópias x, y, cada passo um QP estritamente
    convexo; a seção remete a [1] para a análise de convergência (**não verificada**). O texto diz "convergent to a
    local optimum" (p. 4) e trata 𝒳 como "a set defined by linear constraints" (eq. (5.1)).

### Onde confirma ou contradiz o plano

- **Confirma** que a f do plano (spec §3.3 linha 137, `equal_risk_cont` do R) é a formulação de referência e que o
  caráter local é inevitável (plano linha 11; BST Ex. 3.2).
- **Confirma** que a solução padrão para restrições lineares e box é uma sequência de QPs convexos, resolvidos pelo
  mesmo algoritmo de QP que o PortOpt já tem. O plano escolheu outra via (AL + L-BFGS projetado) sem registrar por quê.
- **Não cobre** a vol-alvo: tanto SCRIP (como implementado e como descrito no livro) quanto BST supõem conjunto viável
  poliédrico/convexo. √(w'Σw) = s é uma igualdade não convexa; linearizá-la dentro do QP é SQP, não SCA, e pede função
  de mérito (N&W cap. 18, não consultado).

### Recomendação

Uma sequência de QPs de Goldfarb-Idnani é substituto sólido do AL **para orçamento, `target_return` e box** (as linhas
ficam exatas, os multiplicadores saem do QP). Para `target_volatility`, manter o AL **apenas** sobre h₃ e usar os QPs como
solver do subproblema (item 1, recomendação 1). Parada por KKT (como o plano já faz, linha 305), não pelas variações
relativas do R, que aceitam pontos não estacionários. Passo γ: Armijo em L em vez do decaimento γ ← γ(1 − ζγ) do R,
porque Armijo dá descida monótona verificável.

---

## 3. Existência e unicidade por ortante (long-short)

### O que as fontes dizem

- **Ortante positivo.** Bruder & Roncalli (2012), §2.3.3 (p. 10): a formulação min √(y'Σy) s.a. Σbᵢ ln yᵢ ≥ c, y ≥ 0
  "demonstrates that the RB portfolio … exists and is unique as long as the covariance matrix Σ is positive-definite"
  (com bᵢ > 0, nota 13). A nota 12 (p. 9) recomenda resolver sem 1'x = 1 e normalizar depois. **Verificado.**
  BST 2013, Lema 2.1 (p. 6): solução única com xᵢ > 0 e Σxᵢ = 1. **Verificado**; o lema diz "positive semidefinite",
  mas a prova usa um minimizador do problema de barreira, que não existe se Σz = 0 para algum z > 0 (então
  F(tz) → −∞): a hipótese correta é Σ definida positiva (derivação própria). Maillard-Roncalli-Teiletche (2010)
  é citado por ambos para o caso ERC (**não lido**).
- **Outros ortantes.** BST 2013, §2.2, eq. (2.13) e **Lema 2.2** (p. 7): para cada β ∈ {−1, 1}ⁿ o problema
  min ½x'Σx − cΣln(βᵢxᵢ), βᵢxᵢ > 0 é convexo com solução única x^β; x^{−β} = −x^β; a normalização x^β/Σx^β_i existe
  sse Σx^β_i ≠ 0, e cai no ortante β ou no complementar conforme o sinal; logo **"at most 2^(n−1)" soluções
  normalizadas**. O exemplo diagonal (2.14)–(2.15) e a Tabela 4.1 (4 soluções para n = 3, pp. 8 e 14) ilustram.
  Palomar (2025) §11.4.1 descreve o mesmo truque de Spinu: w = s ∘ w̃, Σ̃ = Diag(s)ΣDiag(s). **Verificado.**

### Onde confirma ou contradiz

- Spec §3.3, linha 135: "sem bounds, o problema tem uma solução por ortante" **está errado** em dois pontos:
  1. com o orçamento 1'w = 1, cada par {s, −s} contribui no máximo uma solução (no ortante s se s'ỹ_s > 0, no −s se
     < 0), então há no máximo 2^(n−1), não 2ⁿ;
  2. em um par onde s'ỹ_s = 0 (carteira "market-neutral", BST, Remark p. 7) não há solução normalizada.
  Sem o orçamento (raios), a afirmação correta é "um raio por par de ortantes opostos".
- Confirma a parte que importa: o ortante positivo sempre tem exatamente uma solução para Σ DP e bᵢ > 0 (s = 1 dá
  s'ỹ = 1'ỹ > 0 sempre), que é o caminho rápido do plano (linha 9) e o Review Focus 3 (linha 77).
- Consequência não documentada: com bounds que admitem pesos negativos e Spinu fora do box, pode existir uma solução
  exata (f = 0, mínimo global) em outro ortante dentro do box enquanto o AL partindo de Spinu devolve um ponto KKT com
  f > 0. BST §2.3 (p. 8): "up to 2^(n−1) different instances of (2.13) may have to be solved before infeasibility can be
  established".

### Recomendação

- Corrigir a frase da spec e escrever em `math.md` (Task 7): "Sem bounds há no máximo 2^(n−1) carteiras com
  contribuições exatamente b e soma 1, uma por par de ortantes opostos cuja soma dos pesos não se anula (Bai,
  Scheinberg e Tütüncü, Lema 2.2); a biblioteca devolve a de pesos positivos, que sempre existe e é única."
- Documentar no header que `OPTIMAL` com bounds negativos é local e pode perder uma solução exata de outro ortante.
- Opcional (baixa prioridade): sejam k os ativos com lo < 0 < hi. Se k ≤ 10, resolver Spinu sobre Diag(s)ΣDiag(s)
  para os 2^k padrões compatíveis com os bounds (cada um um Newton barato) e, se algum normalizado cai no box, devolvê-lo
  (f = 0 prova otimalidade global). Isso transforma um caso "local" em global a custo previsível.

---

## 4. Spinu 2013 e Griveau-Billion-Richard-Roncalli 2013

### O que as fontes dizem

- **Spinu (2013)**, SSRN 2297383: o SSRN devolveu 403 e o texto **não foi lido**. O que segue vem de GRR (2013,
  arXiv:1311.4057v1, Apêndice A.3, pp. 8–9, **verificado**) e do código de referência (**verificado**):
  - Função (5): f(y) = ½y'Cy − Σbᵢ ln yᵢ com C a **correlação** (não a covariância), gradiente Cy − b/y e Hessiana
    C + diag(b/y²); x* = σᵢ⁻¹y*ᵢ/Σσⱼ⁻¹y*ⱼ.
  - Esquema de Nesterov para funções autoconcordantes: fase amortecida x ← x − Δx/(1 + λ_f) enquanto λ_f ≥ β, fase
    quadrática com passo cheio; λ* = (3 − √5)/2; Spinu "propõe substituir λ_f por δ_f = ‖Δy/y‖∞" e toma
    β = 0,95λ* (nota 9); partida x₀ = (1'C1)^(−1/2)·1.
  - GRR (p. 2) e A.3 (p. 8) dizem que Spinu "notices that the associated Lagrange function is self-concordant".
- **Código de referência** `src/newton_nesterov.cc` (`riskParityPortfolio`, commit `7d8852c`): λ* = 0,3628676
  (= 0,95·(3−√5)/2), fase amortecida `xk - d/(1 + dx)` com `dx = (d/xk).maxCoeff()` (máximo **com sinal**, não
  ‖·‖∞), troca de fase por `lambdak < lambda_star`; fase quadrática com passo cheio e parada por
  max|cᵢ − bᵢ| < tol (tol = 1e-8), o teste do decremento comentado.
- **CCD (GRR §2).** Atualização coordenada fechada, eq. (4) (raiz positiva de uma quadrática; positividade preservada),
  convergência pelo Teorema 5.1 de Tseng (2001) (Remark 2, p. 3). Tabelas 1–3 (pp. 5): Newton e CCD convergem em 100%
  dos casos; Jacobi e SQP não; CCD mais rápido que Newton para n > 250 em Gauss/C/Fortran, Newton mais rápido em
  Matlab/R (Remark 4). Remark 4: resolver o sistema de Newton por Cholesky, não por inversa (o plano já faz, linha 215).

### Autoconcordância: derivação própria

Para φ(y) = −b ln y: φ'' = b/y², φ''' = −2b/y³, logo |φ'''| = (2/√b)(φ'')^{3/2}. A definição (nota 7 de GRR, de
Nesterov 2004) com M_f = 2 ("padrão") vale sse b ≥ 1. Com Σbᵢ = 1 e n ≥ 2, todo bᵢ < 1: **F não é autoconcordante
padrão**, o que confirma o plano (linha 51). Mas F é autoconcordante com M_F = 2/√(min bᵢ) (a parte quadrática tem
terceira derivada nula e a soma preserva o maior M), e G = F/min bᵢ é padrão. A direção de Newton de G é a de F; o
decremento é λ_G = λ_F/√(min bᵢ). Então o passo amortecido 1/(1 + λ_G) **tem** garantia, mas a cota de iterações da fase
amortecida é proporcional a (G(y₀) − G*) = (F(y₀) − F*)/min bᵢ, que com bᵢ = 1e-6 (Review Focus 4) é inútil. Se Spinu
usou o passo com λ_F e M = 2, a garantia formal não se aplica; não foi possível conferir no texto. O passo
1/(1 + ‖Δy/y‖∞) garante positividade: |Δyᵢ|/(yᵢ(1+δ)) ≤ δ/(1+δ) < 1.

### Onde confirma ou contradiz

- Confirma: Newton com Cholesky (linha 215), Armijo com fração até a fronteira (linha 51; para F fortemente convexa com
  Hessiana Lipschitz no conjunto de nível, Newton com backtracking converge globalmente e quadraticamente no fim, teoria
  padrão de Boyd & Vandenberghe §9.5, **não verificado** no texto, o PDF recusou TLS).
- **Problema no critério de parada** (linha 215): λ² ≤ 16nε(1+|F|) testa a subotimalidade F − F* ≈ λ²/2. Se o
  iterado devolvido é aquele onde o teste passou (como no Algoritmo 9.5 de B&V, de memória), o erro no ponto é
  ‖y − y*‖_H ≈ λ ≈ √(16nε) ≈ 1e-7, e como ‖δy‖_H ≥ √bᵢ|δyᵢ|/yᵢ, o erro **relativo** do peso i pode chegar a
  1e-7/√bᵢ (1e-4 com bᵢ = 1e-6). Isso falharia `Spinu.DiagonalClosedForm` a 1e-12 (linha 204) e a faixa de
  `UnevenBudgetsConditioningSweep`. Com o passo cheio dado depois do teste, a convergência quadrática leva o erro a
  ~λ², que é o desejado.
- **Armadilha do código de referência:** `dx` com sinal; se todo dᵢ/xᵢ ≤ −1, 1 + dx ≤ 0 e o passo não preserva y > 0.
  Usar ‖Δy/y‖∞ (como GRR descreve Spinu) ou a fração até a fronteira do plano.
- Normalização por correlação: Spinu trabalha com C e reescala por σ; isso melhora a escala da Hessiana (diag de C é 1)
  e é barato. O plano trabalha com Σ.

### Recomendação

1. Task 2: dar o passo de Newton **e depois** parar quando o teste vale (ou testar após o passo), e medir o erro
   relativo por peso contra a forma fechada diagonal com bᵢ de 1e-6 a 1; alternativa: parar por ‖Δy/y‖∞ ≤ c·ε, que é
   adimensional e controla o erro relativo de cada peso diretamente.
2. Resolver em correlação (y sobre D⁻¹ΣD⁻¹, D = diag(σ)) e reescalar, como Spinu.
3. Partida: (1'C1)^(−1/2)·1 (Spinu via GRR) ou √(Σb/1'Σ1)·1 (código R), ambas com y'Σy = 1 já na escala certa.
4. CCD não é necessário para n moderado; registrar como alternativa para n grande.

---

## 5. Conjunto viável de um ponto, LICQ, restauração e tolerâncias

### Dois casos diferentes no Review Focus 1

**(a) n = 2 com `target_return` e μ₁ ≠ μ₂** (teste `TwoAssetsWithTargetReturnIsThePoint`, linha 310). As igualdades
1'w = 1 e μ'w = r têm gradientes 1 e μ, linearmente independentes: **LICQ vale**, o ponto é isolado mas regular, e os
multiplicadores existem e são únicos (∇f ∈ span{1, μ} = ℝ²). Não há degenerescência; o AL converge (N&W Teorema 17.5
pede só LICQ e segunda ordem, p. 517) e a restauração por Gauss-Newton é exata. Recomendação: detectar "linhas lineares
de posto n" (vértice) e devolver o ponto diretamente, sem AL; mais geralmente, sempre que as linhas ativas fixam w.

**(b) s = σ_min** (teste `VolatilityEqualMinimumIsMinimumVariance`, linha 313). Derivação própria:

- Seja N a matriz das normais das linhas lineares ativas no ponto de mínima variância w_mv. O KKT da mínima variância
  dá 2Σw_mv = Nu, logo ∇h₃ = Σw_mv/σ_min ∈ span(N): **LICQ falha** (e MFCQ, que para igualdades é o mesmo).
- A estacionaridade do RISK_PARITY pede ∇f(w_mv) = Nλ + λ₃∇h₃ ∈ span(N). Sem bounds e sem `target_return`,
  span(N) = span{1} e a identidade w'∇f = 0 força ∇f(w_mv) = 0, falso em geral (w_mv não é carteira de risk parity).
  Então **não existem multiplicadores de KKT**: w_mv é ponto de Fritz John com multiplicador 0 em f, não ponto KKT.
  Só quando span(N) = ℝⁿ (vértice) ou por coincidência eles existem.
- Para s = σ_min + δ, o conjunto viável é uma "esfera" em torno de w_mv de raio ~√δ (no métrico de Σ), então
  f*(s) ≈ f(w_mv) − C√(s − σ_min) e **∂f*/∂s → −∞**. Pela convenção do plano (multiplicador = derivada, linha 29),
  λ_vol = −∞, exatamente o caso "frontier vertical" que o MAX_RETURN já trata (`math.md` 208–216: λ_vol = −inf e os
  demais NaN).
- A fonte primária confirma o comportamento do AL nesse caso: ABMS 2007, Exemplo 2 (p. 15): minimizar x s.a. x² = 0,
  x³ = 0, x⁴ = 0 (conjunto viável {0}, gradientes dependentes, nenhum ponto KKT); o ALGENCAN converge para x = 0 com
  ρ crescendo de 2,5e-5 a 2,5e5 e multiplicadores grandes. Isto é: o AL **não** entra em laço infinito, mas o ponto
  sai com ρ ilimitado e multiplicadores sem sentido. O plano evita o AL nesse caso (linha 303), o que é correto.

**Contradição com o plano.** A linha 61 ("`OPTIMAL`, multiplicadores de mínima norma") e a linha 313 ("o
`expect_risk_parity_kkt` aceita multiplicadores de mínima norma") supõem que existam multiplicadores; os de mínima
norma por QR deixam um resíduo de estacionaridade igual à componente de ∇f(w_mv) fora de span(N), em geral O(1). O teste
falhará, ou passará só se o helper não checar estacionaridade nesse caso.

### Restauração por Gauss-Newton

O passo de mínima norma d = −J⁺h (J das igualdades restritas às livres) é Newton para sistemas subdeterminados, com
convergência quadrática quando J tem posto linha completo (teoria padrão; referência específica **não verificada**).
Ele falha exatamente quando J perde posto, que é o caso (b). Perto de (b) (s/σ_min − 1 pequeno mas acima do limiar),
o ângulo entre ∇h₃ e span(N) é ~√(2δ) e cond(J) ~ 1/√(2δ): com δ = 1e-12, cond(J) ~ 7e5, e os multiplicadores por QR
perdem ~6 dígitos. O ABMS (p. 3) lê multiplicadores grandes como sinal de quase-falha de CQ. Isso não é um defeito do
plano, mas a faixa de KKT relativa pode não ser atingível ali: deve ser medido na Task 5, não afrouxado.

### A tolerância 1e-12 em s × σ_min

Derivação própria do erro de σ_min calculado: em w_mv, σ é estacionária ao longo das direções viáveis, então o erro de
ŵ dentro do conjunto viável afeta σ só em segunda ordem. O termo de primeira ordem vem do resíduo das linhas ativas:
se N'ŵ − b = ρ_r, então σ(ŵ) − σ_min ≈ ∇σ'(ŵ − w_mv) = u'ρ_r/(2σ_min), e como u'b = 2σ_min² (complementaridade),
δσ/σ ≈ (u'ρ_r)/(u'b), que é da ordem da tolerância relativa de viabilidade do QP, **1e-12**
(`src/qp/goldfarb_idnani.cpp:22`, `kRelativeTolerance`, usada em `feasibility_tolerance`, linha 553). A avaliação de
√(w'Σw) acrescenta ~n·u relativo (cota de produto interno). Portanto 1e-12 é o nível certo, **mas pela razão de ser a
tolerância de viabilidade do QP**, não por arredondamento de σ. Evidência do próprio repositório: `validation.md`
linhas 203–207 registra 17 de 320 casos em que σ_min de Q(0) diferia em mais de 1e-12 do σ_min do solve de
`MIN_VARIANCE` (rotas diferentes do mesmo problema). Logo um usuário que passa s = σ_min obtido de `MIN_VARIANCE`
recebe `INFEASIBLE` se as duas rotas não forem a mesma chamada.

### |Σbᵢ − 1| ≤ n·ε

- Rump (2012), BIT, eq. (1.1) (p. 1, **verificado** no preprint "accepted for publication in BIT, May 18, 2011"): a
  estimativa clássica de Higham para soma recursiva é |float(Σpᵢ) − Σpᵢ| ≤ γ_{n−1}Σ|pᵢ|, γ_k = ku/(1 − ku),
  nu < 1; eq. (1.2): Δ ≤ (n − 1)u·Ẽ vale para todo n. Para double, u = 2⁻⁵³ (§2), e `DBL_EPSILON` = 2u. O livro de
  Higham não foi consultado (número de seção **não verificado**).
- Derivação própria: com bᵢ > 0, Σ|bᵢ| ≈ 1. Se cada bᵢ é o arredondamento de um valor exato com soma 1,
  Σ|b̂ᵢ − bᵢ| ≤ u·Σbᵢ = u; a soma em ponto flutuante acrescenta ≤ (n − 1)u(1 + O(u)). Total ≤ n·u + O(nu²). Se o
  usuário normalizou b = a/Σa em ponto flutuante, soma-se mais (n − 1)u do cálculo de Σa: ≈ (2n − 1)u. Ambos ≤ 2n·u
  = n·ε. **Confirma** a cota do plano (linha 57) com folga ≥ 1 no pior caso de normalização e 2 no caso comum. A soma
  por pares (NumPy) só melhora. Ressalva: a checagem precisa ser feita sobre os bᵢ exatamente como recebidos e em
  qualquer ordem (a cota não depende da ordem). O teste `CapiSpec.RiskBudgetDomain` (1 ± nε aceito, 1 + 4nε rejeitado,
  linha 391) é coerente.

### Recomendação

1. Separar o Review Focus 1 em (a) e (b); (a) é regular e pode sair direto; (b) é o caso degenerado.
2. Em s = σ_min: pesos = w_mv, `OPTIMAL`; multiplicadores pela regra do MAX_RETURN em σ_min — λ_vol = −∞ e os demais
   NaN quando o resíduo de estacionaridade dos multiplicadores de mínimos quadrados está acima do ruído; os
   multiplicadores do QP quando o ponto é vértice (então λ_vol é indeterminado: reportar o de mínima norma só nesse
   caso e documentar). `expect_risk_parity_kkt` deve aceitar esse padrão explicitamente, como `expect_user_kkt` faz
   para o MAX_RETURN.
3. Definir `kMinimumVolatilityTolerance` igual ao `kRelativeTolerance` do QP (uma constante, uma justificativa),
   e no teste obter σ_min pela mesma rota (`solve_mean_variance`) que o RISK_PARITY usa.
4. Task 5: medir cond(J) e a faixa de KKT para s/σ_min − 1 ∈ {1e-12, 1e-10, …, 1e-2}.

---

## 6. Planos 1–3

### 6a. Leverage como família implícita de 2ⁿ linhas no Goldfarb-Idnani

- **Tibshirani (1996)**, JRSS-B 58(1), §6 "Algorithms for finding lasso solutions", pp. 277–279 (**verificado**, cópia
  JSTOR): Σ|βⱼ| ≤ t é escrito como 2^p restrições δᵢ'β ≤ t; o algoritmo começa com E = {i₀}, δ_{i₀} = sign(β̂_LS),
  resolve o mínimos quadrados com G_Eβ ≤ t1 e, enquanto Σ|β̂ⱼ| > t, acrescenta δ = sign(β̂). "This procedure must
  always converge in a finite number of steps since one element is added to the set E at each step, and there is a
  total of 2^p elements." Uma variante que remove elementos é "more efficient but it is not clear how to establish its
  convergence". Média observada de 0,5p a 0,75p iterações. Alternativa de Gay: β = β⁺ − β⁻ (2p variáveis, 2p + 1
  restrições).
- **Brodie et al. (2009)**, PNAS 106(30):12267–12272 (via arXiv:0708.0046, **verificado**): penalidade ℓ1 resolvida por
  homotopia/LARS adaptada a Aw = a (Apêndice); com pesos não negativos e orçamento, a penalidade é constante (§3), isto
  é, ℓ1 só age sobre posições vendidas (§2). Não usa a família de sinais.
- **DeMiguel, Garlappi, Nogales, Uppal (2009)**, Management Science: título verificado; conteúdo **não lido**.
- **Goldfarb & Idnani (1983)**, Math. Programming 27:1–33: bibliografia verificada (dblp/MaRDI); o texto (paywall) **não
  foi lido**. De memória, **não verificado**: o método dual mantém otimalidade dual, adiciona uma restrição violada por
  vez, e o valor de f cresce estritamente em cada passo cheio, o que impede repetir um conjunto ativo; como há um número
  finito de conjuntos ativos, termina.

Avaliação: o PortOpt (`math.md` 113–120; `src/qp/goldfarb_idnani.cpp` 209–221, 254–270) é o algoritmo de Tibshirani
com o GI no lugar do Lawson-Hanson: a linha gerada é a de s = sign(x), a mais violada. A dúvida de Tibshirani sobre
remover elementos não se aplica da mesma forma: no GI as remoções são passos parciais governados pela viabilidade dual,
e o argumento de monotonia (se correto como lembrado) vale para qualquer família finita, inclusive 2ⁿ linhas,
desde que o oráculo devolva de fato a linha mais violada (o que `sign(x)` garante, `math.md` 116–118). A regra contra
duplicar uma linha ativa (`validation.md` 268–276) cobre o caso degenerado que um oráculo inexato criaria. A divisão de
Gay é exatamente a que a spec rejeitou por singularizar a Hessiana (spec linha 167), o que bate com a exigência de G DP
do GI. **Confirma**; recomendação: citar Tibshirani §6 em `math.md` como precedente.

### 6b. Guarda de ruído do κ no MAX_SHARPE

- O modelo (`validation.md` 356–362; `homogeneous.cpp` 204, 320) propaga um erro retroativo |ΔΣ| ≤ nε|Σ| e |Δv| ≤ εm
  em primeira ordem até κ = σ_y²z'v. É o raciocínio de Higham (erro retroativo → perturbação de primeira ordem
  componente a componente), mas: (i) a cota rigorosa do erro retroativo do Cholesky é em termos de |Rᵀ||R| (Higham,
  cap. 10, de memória, **não verificado**), que pode exceder |Σ| componente a componente; (ii) a constante n é
  heurística. A validação empírica compensa: erro máximo medido 0,040 κ_noise contra o limiar 16, margem 400
  (`validation.md` 395–401). **Confirma com ressalvas.**
- **Merton (1972)** (JFQA 7(4), **não lido**): a tangente existe e está no ramo superior sse r_f < B/A (o retorno da
  carteira de mínima variância global); isso equivale a 1'Σ⁻¹(μ − r_f1) = B − A·r_f > 0, que é o sinal de κ sem bounds
  (derivação direta, coerente com `math.md` 275 e 283).
- **Britten-Jones (1999)**, J. Finance 54(2):655–671, "The Sampling Error in Estimates of Mean-Variance Efficient
  Portfolio Weights" (abstract **verificado**; o título correto usa "Estimates", não "estimated"): pesos eficientes como
  coeficientes de MQO com testes t/F exatos sob normalidade; o erro amostral em 11 índices em 20 anos é "large". Em
  consequência, o sinal de 1'Σ⁻¹μ̃ estimado é estatisticamente incerto muito antes de ser numericamente incerto; a
  guarda trata só o segundo, o que é o papel correto de uma biblioteca determinística.
- Recomendação: registrar em `validation.md` a ressalva (i) e a frase de que `UNBOUNDED`/`OPTIMAL` perto da fronteira é
  decidido pelos dados como dados, sem significado estatístico. Nenhuma mudança de código.

### 6c. MAX_RETURN: caminho paramétrico × Critical Line Algorithm

- **Niedermayer & Niedermayer (2007)**, Discussion Paper 07-01, Univ. Bern (**verificado**, pp. 1–12): Proposição 1 (p. 5),
  entre turning points as livres resolvem o problema sem bounds no subconjunto F; **Corolário 1** (p. 6): combinação
  convexa de dois turning points vizinhos é carteira de mínima variância restrita; Proposição 2 (p. 7): μ'w é afim e
  crescente em λ entre turning points; §3.1: partida no turning point de maior retorno (exige 1'l ≤ 1 ≤ 1'u, nota 8);
  §3.2: próximo λ por teste de razão exato (casos a/b, eqs. (8)–(14)); nota 7: turning points coincidentes são
  degenerados e tratados à parte (Markowitz 1956/1959, **não lido**); nota 10: λ⁽ⁱ⁾ = λ_current − ε por imprecisão
  numérica é armadilha conhecida. O CLA deles cobre só orçamento + box.
- O PortOpt (`math.md` 140–200) usa a mesma estrutura (w(t) afim por partes com breakpoints nas trocas de conjunto
  ativo; ponto reportado como interpolação das pontas do segmento, que é o Corolário 1; inclinação dos multiplicadores
  por mínimos quadrados em vez de diferença de soluções, o análogo das eqs. (8)–(9)). Diferenças: (i) dobra + bissecção
  em vez do teste de razão exato (mais solves de QP por breakpoint, mas robusto à nota 10 sem tratamento especial);
  (ii) parte de t = 0 (mínima variância) e sobe, em vez de partir do topo, o que dispensa o box limitado; (iii) cobre
  caixa, leverage e faces deslocadas, fora do escopo do CLA. **Confirma.** Recomendação: citar N&N; opcional, quando os
  dois extremos do bracket têm o mesmo conjunto ativo, calcular o próximo breakpoint pelo teste de razão em vez de
  bissectar (economia, não correção).

### 6d. Homogeneização com leverage (derivação própria)

Seja v ∈ ℝⁿ, L > 0 e considere os conjuntos
W = {w : 1'w = 1, ‖w‖₁ ≤ L, v'w > 0} e Y = {y : v'y = 1, L·1'y − ‖y‖₁ ≥ 0}.

1. w ∈ W ⇒ y = w/(v'w) ∈ Y: v'y = 1; 1'y = 1/(v'w) > 0; ‖y‖₁ = ‖w‖₁/(v'w) ≤ L/(v'w) = L·1'y.
2. y ∈ Y ⇒ y ≠ 0 (v'y = 1), e L·1'y ≥ ‖y‖₁ > 0, logo κ = 1'y > 0; w = y/κ tem 1'w = 1, ‖w‖₁ = ‖y‖₁/κ ≤ L e
   v'w = 1/κ > 0, isto é, w ∈ W.
3. As duas aplicações são inversas e a razão é v'w/√(w'Σw) = (1/κ)/(√(y'Σy)/κ) = 1/√(y'Σy). Maximizar a razão em W
   equivale a minimizar y'Σy em Y. Os bounds homogeneizados entram na mesma interseção sem alterar o argumento.
4. **Cota de κ**: 1 = v'y ≤ ‖v‖∞‖y‖₁ ≤ ‖v‖∞·L·κ, logo **κ ≥ 1/(L‖v‖∞) > 0** em todo Y. Com leverage, a linha
   1'y ≥ 0 é redundante, κ → 0 é impossível, Y é fechado e a razão atinge o supremo: nunca `UNBOUNDED`.

Isso prova o que `validation.md` 370–374 e `cone_bounded` (`homogeneous.cpp` 279) afirmam para "a leverage is set", e
dá uma cota inferior explícita de κ em todo ponto viável, que pode substituir o argumento qualitativo (não medida
contra κ_noise aqui). Cornuéjols & Tütüncü (2007, §8.2) e Choueifaty & Coignard (2008) / Choueifaty, Froidure, Reynier (2013) são
a origem da homogeneização para Sharpe e diversificação (**não lidos**; capítulo citado de segunda mão). **Confirma.**
Recomendação: pôr os passos 1–4 em `math.md` (seção MAX_SHARPE) como a prova da linha com leverage.

---

## Bibliografia

| Referência | URL | Status |
|---|---|---|
| Nocedal, J.; Wright, S. J. *Numerical Optimization*, 2ª ed., Springer, 2006, cap. 17 | https://pages.mtu.edu/~struther/Courses/OLD/Other/LessOld/5630_11/NocedalAndWright/Chpt17.pdf (cópia de curso, não do editor) | Verificado: §17.3–17.4, eqs. (17.36)–(17.53), Framework 17.3, Alg. 17.4, Teoremas 17.5–17.6, pp. 514–524. Cap. 18 não consultado. |
| Andreani, R.; Birgin, E. G.; Martínez, J. M.; Schuverdt, M. L. "On augmented Lagrangian methods with general lower-level constraints", SIAM J. Optim. 18(4):1286–1309, 2007 | https://users.ime.usp.br/~egbirgin/publications/abms2.pdf | Verificado (preprint do autor): Alg. 3.1, Teorema 4.1, §6 (τ = 0,5, γ = 10), Exemplos 1–2, p. 17. |
| Andreani, R.; Martínez, J. M.; Svaiter, B. F. "A new sequential optimality condition for constrained optimization and algorithmic consequences", SIAM J. Optim. 20(6), 2010 | https://bv.fapesp.br/en/publicacao/34470/a-new-sequential-optimality-condition-for-constrained-optimi | Título, autores e abstract verificados; periódico e relação com CAKKT não verificados. |
| Conn, Gould, Toint, LANCELOT (1992) | — | Não consultado; descrito via N&W §17.4. |
| Byrd, R. H.; Lu, P.; Nocedal, J.; Zhu, C. "A limited memory algorithm for bound constrained optimization", SIAM J. Sci. Comput. 16(5):1190–1208, 1995 | https://users.iems.northwestern.edu/~nocedal/lbfgsb.html | Bibliografia e abstract verificados; detalhes do algoritmo não verificados. |
| Zhu, C.; Byrd, R. H.; Nocedal, J. "Algorithm 778: L-BFGS-B", ACM TOMS 23(4):550–560, 1997; Morales, J. L.; Nocedal, J. "Remark on Algorithm 778", ACM TOMS, 2011 | https://users.iems.northwestern.edu/~nocedal/lbfgsb.html; https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html | Bibliografia e critério pgtol verificados. |
| Bertsekas, D. P. "Projected Newton methods for optimization problems with simple constraints", SIAM J. Control Optim. 20(2):221–246, 1982 | https://www.mit.edu/~dimitrib/ProjectedNewton.pdf | Verificado: §1–2, Fig. 1, Proposição 1, eq. (32), pp. 221–228. |
| Kim, D.; Sra, S.; Dhillon, I. "Tackling box-constrained optimization via a new projected quasi-Newton approach", SIAM J. Sci. Comput. 32(6):3548–3563, 2010 | https://optml.lids.mit.edu/papers/pqnj_sisc.pdf | Bibliografia e abstract verificados (busca); texto não lido. |
| Feng, Y.; Palomar, D. P. "SCRIP: Successive convex optimization methods for risk parity portfolio design", IEEE TSP 63(19):5285–5300, 2015 | https://palomar.home.ece.ust.hk/papers/2015/FengPalomar-TSP2015%20-%20risk_parity_portfolio.pdf | **Não lido** (certificado TLS recusado). Bibliografia verificada (CRAN). |
| Palomar, D. P. *Portfolio Optimization: Theory and Application*, Cambridge, 2025, cap. 11 | https://portfoliooptimizationbook.com/book/11-RPP.html | Verificado: §11.4, §11.4.1, §11.7. |
| Pacote R `riskParityPortfolio` 0.2.2 (CRAN, 2021-05-31) | https://cran.r-project.org/web/packages/riskParityPortfolio/refman/riskParityPortfolio.html; https://github.com/dppalomar/riskParityPortfolio (commit 7d8852c) | Verificado: `R/riskParityPortfolio.R` 18–181, `src/newton_nesterov.cc`, `src/cyclical_coordinate_descent.cc`. |
| Bai, X.; Scheinberg, K.; Tütüncü, R. "Least-squares approach to risk parity in portfolio selection", Quantitative Finance 16(3):357–376, 2016 | https://optimization-online.org/wp-content/uploads/2013/10/4089.pdf | Verificado no preprint (rev. 21/10/2013): Lemas 2.1–2.3, 3.1, Exemplos 3.1–3.2, Algs. 1–3, pp. 2–17. Versão publicada não lida (paywall). |
| Griveau-Billion, T.; Richard, J.-C.; Roncalli, T. "A fast algorithm for computing high-dimensional risk parity portfolios", arXiv:1311.4057, 2013 | https://arxiv.org/pdf/1311.4057 | Verificado: §1–3, Apêndice A.3, pp. 1–9. |
| Spinu, F. "An algorithm for computing risk parity weights", SSRN 2297383, 2013 | https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2297383 | **Não lido** (SSRN 403). Conteúdo via GRR A.3 e código R. |
| Bruder, B.; Roncalli, T. "Managing risk exposures using the risk budgeting approach", MPRA 37749, 2012 | https://mpra.ub.uni-muenchen.de/37749/2/MPRA_paper_37749.pdf | Verificado: §2.1–2.3.4, pp. 3–10. |
| Maillard, S.; Roncalli, T.; Teiletche, J. "The properties of equally weighted risk contribution portfolios", J. Portfolio Management 36(4):60–70, 2010 | — | Não lido; citado via Bruder-Roncalli e Palomar. |
| Roncalli, T. *Introduction to Risk Parity and Budgeting*, Chapman & Hall/CRC, 2013 | — | Não lido. |
| Tibshirani, R. "Regression shrinkage and selection via the lasso", JRSS-B 58(1):267–288, 1996 | https://webdoc.agsci.colostate.edu/koontz/arec-econ535/papers/Tibshirani%20(JRSS-B%201996).pdf (cópia JSTOR de curso) | Verificado: §6, pp. 277–279. |
| Goldfarb, D.; Idnani, A. "A numerically stable dual method for solving strictly convex quadratic programs", Math. Programming 27:1–33, 1983 | https://doi.org/10.1007/BF02591962 | Bibliografia verificada; texto **não lido** (paywall). |
| Brodie, J.; Daubechies, I.; De Mol, C.; Giannone, D.; Loris, I. "Sparse and stable Markowitz portfolios", PNAS 106(30):12267–12272, 2009 | https://ar5iv.arxiv.org/html/0708.0046 | Verificado (versão arXiv): §2–3, Apêndice. |
| DeMiguel, V.; Garlappi, L.; Nogales, F. J.; Uppal, R. "A generalized approach to portfolio optimization: improving performance by constraining portfolio norms", Management Science, 2009 | — | Só o título verificado. |
| Niedermayer, A.; Niedermayer, D. "Applying Markowitz's critical line algorithm", Discussion Paper 07-01, Univ. Bern, 2007 | https://repec.vwiit.ch/dp/dp0701.pdf | Verificado: pp. 1–12. |
| Markowitz, H. "The optimization of a quadratic function subject to linear constraints", Naval Research Logistics Quarterly 3, 1956 | — | Não lido. |
| Britten-Jones, M. "The sampling error in estimates of mean-variance efficient portfolio weights", J. Finance 54(2):655–671, 1999 | https://doi.org/10.1111/0022-1082.00120 | Abstract verificado. |
| Merton, R. C. "An analytic derivation of the efficient portfolio frontier", JFQA 7(4), 1972 | — | Não lido; afirmação usada é derivação direta. |
| Rump, S. M. "Error estimation of floating-point summation and dot product", BIT, 2012 (aceito em 18/05/2011) | https://www.tuhh.de/ti3/paper/rump/Ru11.pdf | Verificado: eqs. (1.1)–(1.2), §2. Volume/páginas não verificados. |
| Higham, N. J. *Accuracy and Stability of Numerical Algorithms*, 2ª ed., SIAM, 2002 | — | Não lido; a cota γ_{n−1} é citada via Rump (1.1). Afirmações sobre o cap. 10 (Cholesky) de memória. |
| Nesterov, Y. *Introductory Lectures on Convex Optimization*, Kluwer, 2004 | — | Não lido; definição de autoconcordância via GRR, nota 7. |
| Boyd, S.; Vandenberghe, L. *Convex Optimization*, 2004, §9.5–9.6 | https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf | Não lido (certificado TLS recusado); afirmações de memória. |
| Tseng, P. "Convergence of a block coordinate descent method for nondifferentiable minimization", JOTA 109(3):475–494, 2001 | — | Não lido; citado por GRR, Remark 2. |
| Cornuéjols, G.; Tütüncü, R. *Optimization Methods in Finance*, Cambridge, 2007, §8.2 | — | Não lido; seção citada de segunda mão. |
| Choueifaty, Y.; Coignard, Y. "Toward maximum diversification", JPM 35(1), 2008; Choueifaty, Y.; Froidure, T.; Reynier, J. "Properties of the most diversified portfolio", J. Investment Strategies 2(2), 2013 | — | Não lidos. |

---

## Mudanças priorizadas no plano (antes das Tasks 3 e 4)

1. **[Alta, Task 4, Task 0 em andamento] s = σ_min não tem multiplicadores de KKT.** Trocar a linha 61 e o teste da
   linha 313: pesos = w_mv; λ_vol = −∞ e demais NaN quando a estacionaridade dos multiplicadores de mínimos quadrados
   não fecha (padrão do MAX_RETURN, `math.md` 208–216); multiplicadores do QP só quando w_mv é vértice. Ajustar
   `expect_risk_parity_kkt` para aceitar esse padrão. Separar `TwoAssetsWithTargetReturnIsThePoint` (regular, LICQ
   vale) do caso degenerado no Review Focus 1.
2. **[Alta, Task 2] Critério de parada do Spinu (linha 215).** Dar o passo de Newton depois do teste (ou parar por
   ‖Δy/y‖∞ ≤ c·ε), senão o erro relativo do peso é ~1e-7/√bᵢ; medir contra a forma fechada diagonal com bᵢ de 1e-6
   a 1. Usar ‖Δy/y‖∞ com valor absoluto se houver fase amortecida (não o máximo com sinal do código R).
3. **[Alta, Tasks 3–4] Decidir o solver do subproblema antes da Task 3.** Recomendado: igualdades lineares + box no
   nível inferior, resolvidas por passos Gauss-Newton com `qp::solve` (Q = 2J̃ᵀJ̃ + τI, Armijo em L), AL só para
   vol(w) = s; isso remove `box_lbfgs` e a restauração de h₁, h₂, e dá os multiplicadores lineares do próprio QP.
   Se a Task 3 for mantida, acrescentar à interface o conjunto ε-ativo (Bertsekas eq. (32)), Armijo ao longo do arco,
   descarte de pares com s'y ≤ 0 e um teste com bound ativo de gradiente nulo.
4. **[Média, Task 4] Regra de ρ e tolerâncias do AL.** Atribuir a regra ao ABMS 2007 Alg. 3.1 (não ao N&W §17.4),
   usar τ = 0,5, γ = 10 e o ρ₁ do §6 como partida (a Task 5 mede), ε_k → 0 em vez de tolerância fixa, multiplicadores
   salvaguardados em [−λ_max, λ_max].
5. **[Média, Task 4] Polimento final por Newton-KKT** com a Hessiana exata no conjunto ativo identificado, para atingir
   a faixa 1e-10 que o AL sozinho não garante (ABMS p. 17).
6. **[Média, Task 4] `kMinimumVolatilityTolerance`** = `kRelativeTolerance` do QP, com a justificativa do item 5 em
   `validation.md`; no teste, σ_min pela mesma rota de `solve_risk_parity`.
7. **[Média, spec e Task 7] Corrigir "uma solução por ortante"** (spec §3.3 linha 135) para "no máximo 2^(n−1),
   uma por par de ortantes opostos com soma não nula; a positiva sempre existe e é única para Σ DP"; documentar que
   `OPTIMAL` com bounds negativos pode perder uma solução exata de outro ortante.
8. **[Baixa, Task 5] Medir** cond(J) e a faixa de KKT para s/σ_min − 1 de 1e-12 a 1e-2; registrar.
9. **[Baixa, docs] Citações nos guias:** Tibshirani 1996 §6 (família de leverage), Niedermayer & Niedermayer 2007
   (caminho paramétrico), a prova 6d com κ ≥ 1/(L‖v‖∞) (homogeneização com leverage), a ressalva |Rᵀ||R| × |Σ| do
   κ_noise.
10. **[Opcional] Varredura de ortantes** quando k ≤ 10 ativos podem trocar de sinal (2^k Spinu), devolvendo uma
    solução exata dentro do box se existir.
