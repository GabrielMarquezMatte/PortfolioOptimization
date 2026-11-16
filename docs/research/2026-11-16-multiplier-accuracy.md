# Pesquisa: precisão dos multiplicadores do RISK_PARITY com volatilidade-alvo

> **Sobre este arquivo.** Nota de apoio à decisão no formato de `docs/research/2026-11-16-risk-parity-and-numerics.md`
> (português, termos técnicos e títulos em inglês). Não altera código, plano nem spec. Data: 2026-11-16.
> Contexto lido: `include/portopt/portopt.h` (entrada RISK_PARITY, linhas 114–127 e 239–295),
> `src/riskparity/{kkt,newton_polish,objective,risk_parity}.cpp`, `docs/guide/validation.md` (seções "The stationarity
> band…", "Forward error…", "The multipliers of the target-volatility path…", "Future work", "The Newton polish") e o
> probe de 80 dígitos `probe6.py` do scratchpad da sessão.
>
> **Convenção de verificação.** "Verificado" = lido na fonte primária (artigo, preprint do autor, documentação ou
> código de referência), com seção/página. "Não verificado" = a fonte não foi lida (paywall, PDF ilegível, não
> encontrada) e a afirmação vem de citação secundária ou de memória. "Derivação própria" = argumento escrito aqui,
> conferível linha a linha. "Medido aqui" = experimento numérico desta nota (seção 0), reproduzível pelo script indicado.

## Resumo

| # | Questão | Veredito | Recomendação | Impacto |
|---|---|---|---|---|
| 1 | Sensibilidade de λ a δw; o erro observado é ~ cond(H)·ε? | **Confirma** a hipótese como diagnóstico: o λ de mínimos quadrados (LS) em **qualquer** w double tem erro ≈ (A')⁺H δw, de ordem ε·\|H\|\|w\|/\|termos de ∇f\| (derivação própria; medido: o LS em fl(w\*) avaliado em 80 dígitos erra tanto quanto o solver). **Refuta** a parte "nenhum w double pode dar": o par KKT (w, λ) não é mal condicionado em λ; só o estimador LS-em-w é. | Não reestimar λ por LS no ponto final; tomar λ = λ_LS + dλ do sistema de Newton-KKT no ponto final. | `newton_polish.cpp` (`settled`, `newton_step`, `newton_polish`). |
| 2 | Refinamento com resíduo em precisão estendida: quando funciona | **Confirma** (Demmel et al. 2006; Tisseur 2001; Carson-Higham 2018): a precisão limite é ‖K⁻¹‖·ψ (erro de avaliação do resíduo) + u, e a convergência exige u·κ(K)·(erro do solver/Jacobiana) ≲ 1. Para λ, o resíduo em double **já basta** até cond(Σ) = 1e8 (ψ pequeno com o Σw compensado); o fator limitante medido passa a ser a precisão da Hessiana analítica. | Fase 1 (só double): dλ final. Fase 2 (opcional, cond ≥ 1e8): resíduo double-double em w_ext = (w, dw) e 2–3 passos, **depois** de medir o erro real da Hessiana. | Fase 1: 1 solve extra por polish. Fase 2: módulo pequeno de aritmética dd. |
| 3 | Double-double/compensado, portável a MSVC | **Confirma** viabilidade: TwoSum (6 flops), TwoProductFMA (2 flops), Dot2 (25n−7 flops; Ogita-Rump-Oishi 2005) bastam; o código já tem Dot2 (`compensated_dot`). Custo do resíduo dd: < 20% de um passo de polish para n ≤ 100. | Template das três somas que cancelam (Σw, Σ(w∘r), 1'w) sobre o par (hi, lo); nada de `long double`. | Fase 2 apenas. |
| 4 | Alternativas e o que outros solvers fazem | OSQP (polishing, §4) resolve **x e y juntos** no KKT reduzido com refinamento iterativo (3 passos) — o análogo de "manter dλ". IPOPT oferece multiplicadores LS (`recalc_y`, opcional, desligado). LAPACK xGERFSX devolve **cotas de erro** certificadas. Nenhuma fonte encontrada relata precisão de multiplicadores em problemas mal condicionados como este (não verificado que não exista). | Fase 1 + documentar a faixa medida; certificação a posteriori só se a Fase 2 for feita (a razão de correções sucessivas dá a cota, Demmel et al. eq. (6)). | Docs (`portopt.h`, `validation.md`). |
| 5 | Caminho QP (Goldfarb-Idnani) | **Sem análogo direto** (derivação própria): H = Σ (não ~cond²·³), e o GI obtém multiplicadores da mesma fatoração que dá x, não por LS em x arredondado. `validation.md` não mede o erro progressivo dos multiplicadores do QP a 1e8 (só estacionaridade e diferenças de inclinação). | Medir uma vez (barato) com a mesma referência; não esperar problema. | Nenhum código. |

---

## 0. Experimento feito para esta nota (medido aqui)

Script: `C:/Users/x/AppData/Local/Temp/claude/5e03637a68c95476/86bc5840-28b8-4a5d-80a9-19ce98cfef7a/scratchpad/exp_mult.py`
(o scratchpad pode ser apagado; o script é o `probe6.py` acrescido das variantes abaixo, ~200 linhas). Ambiente:
`python/.venv`, `PORTOPT_LIB=build/gcc-release/libportopt.dll` (commit `5882596`), PATH com `C:/msys64/ucrt64/bin`.
Casos do `probe6.py`: sem box, orçamento igual, s = 1,2·σ_min, (n, seed) ∈ {(5,7), (5,11), (8,3)}, cond(Σ) ∈ {1e4, 1e6,
1e8, 1e10}. Referência (w\*, λ\*): 12 passos de Newton em `decimal` com 80 dígitos a partir do resultado do solver.
Erros relativos como no `probe6.py` (λ_vol relativo a \|λ_vol\*\|; λ_budget relativo a max(\|λ_b\*\|, \|λ_vol\*\|s)).

Variantes (todas sem mudar w, exceto E4/E5):

- **E0** λ do solver (reproduz `validation.md` 853–856).
- **E1** λ_LS em fl(w\*) (w\* arredondado corretamente a double), **avaliado em 80 dígitos**: isola o efeito de w ser double.
- **E2** λ_LS no w do solver, avaliado em 80 dígitos: isola a aritmética do LS do solver.
- **E3** um sistema de Newton-KKT no ponto do solver: resíduo num modelo da aritmética double da biblioteca (Σw e w'Σw
  corretamente arredondados, como o Dot2; o resto em float64), matriz KKT = Jacobiana de 80 dígitos arredondada a
  double, escalada como `newton_system` (linhas de A à escala de H), QR Householder; reporta λ + dλ.
- **E3-H** idem com a parte H perturbada por η = 2,2e-16·cond(Σ) (modelo do erro da Hessiana analítica: `validation.md`
  755–757 mede 3,6e-8 do maior elemento a 1e8): "global" = η·max\|H\| por elemento; "entrywise" = η\|H_ij\|.
- **E3x** como E3 com o resíduo calculado em 80 dígitos e arredondado a double (modelo de resíduo double-double).
- **E4** λ_LS em w_ext = w + dw (dw do E3x), avaliado em 80 dígitos.
- **E5** 3 passos de refinamento guardando w_ext = (hi, lo) e λ, resíduo dd, H perturbada "global".

| cond | caso | erro w | E0 λ_vol / λ_b | E1 | E3 | E3-H global | E3-H entrywise | E3x | E4 | E5 após 1 / 3 passos (λ_vol) | cond(K) escalada |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1e4 | 5,7 | 1,2e-16 | 3,5e-11 / 2,8e-11 | 2,0e-11 / 1,6e-11 | 2,1e-15 / 8,9e-15 | igual E3 | igual E3 | 2,9e-17 / 1,5e-17 | 7e-25 | 1,7e-21 / 7e-37 | 7,0e4 |
| 1e4 | 5,11 | 1,8e-16 | 9,5e-11 / 7,8e-11 | 6,1e-10 / 5,2e-10 | 6,8e-14 / 1,3e-14 | igual | igual | 1,2e-17 / 1,6e-17 | 9e-25 | 2,4e-21 / 4e-37 | 1,5e5 |
| 1e4 | 8,3 | 1,2e-15 | 8,3e-12 / 5,8e-12 | 3,7e-12 / 2,6e-12 | 7,3e-15 / 6,1e-15 | igual | igual | 6,9e-17 / 1,9e-17 | 6e-26 | 9,0e-21 / 5e-33 | 1,1e6 |
| 1e6 | 5,7 | 0 | 3,5e-8 / 8,8e-8 | 3,5e-8 / 8,8e-8 | 1,3e-14 / 2,7e-13 | 1,3e-14 / 2,7e-13 | igual | 4,1e-17 / 6,7e-17 | 1e-20 | 5,3e-18 / 7e-24 | 6,7e6 |
| 1e6 | 5,11 | 1,3e-16 | 6,9e-8 / 1,0e-7 | 1,5e-8 / 2,2e-8 | 1,0e-12 / 1,3e-12 | 1,0e-12 / 1,3e-12 | igual | 5,1e-18 / 9,8e-17 | 4e-20 | 2,2e-16 / 5e-23 | 3,1e7 |
| 1e6 | 8,3 | 2,6e-14 | 8,4e-8 / 2,1e-7 | 1,4e-7 / 3,5e-7 | 5,3e-13 / 1,5e-12 | 3,2e-14 / 1,1e-12 | 4,3e-13 / 1,2e-12 | 1,3e-17 / 1,9e-17 | 2e-20 | 2,4e-13 / 1,2e-14 | 8,5e9 |
| 1e8 | 5,7 | 2,1e-16 | 3,1e-5 / 1,4e-4 | 4,9e-5 / 2,3e-4 | 7,1e-13 / 3,4e-12 | 9,1e-11 / 1,1e-10 | 5,6e-12 / 3,2e-11 | 9,0e-17 / 5,8e-17 | 1e-16 | 9,6e-11 / 6,5e-13 | 4,5e8 |
| 1e8 | 5,11 | 2,2e-16 | 1,9e-4 / 4,5e-4 | 2,5e-4 / 5,7e-4 | 1,6e-11 / 1,8e-10 | 1,6e-10 / 4,8e-10 | 8,1e-12 / 5,8e-13 | 3,6e-18 / 8,9e-17 | 3e-15 | 9,3e-11 / 2,2e-11 | 2,7e9 |
| 1e8 | 8,3 | 3,6e-13 | 7,4e-4 / 3,6e-3 | 8,9e-4 / 4,4e-3 | 5,3e-11 / 1,2e-10 | 4,3e-10 / 3,7e-10 | 9,4e-11 / 2,0e-10 | 1,6e-14 / 7,5e-14 | 1e-15 | 1,0e-9 / 1,6e-10 | 1,7e13 |
| 1e10 | 5,7 | 1,5e-16 | 3,1e-1 / 2,4 | 2,4 / 19 | 1,5e-10 / 1,9e-9 | 5,3e-6 / 1,0e-4 | 2,2e-5 / 6,8e-5 | 2,2e-16 / 1,1e-14 | 1e-10 | 5,4e-6 / 7,3e-7 | 3,4e10 |
| 1e10 | 5,11 | 5,0e-16 | 2,1 / 8,1 | 1,7 / 6,8 | 6,6e-10 / 1,1e-8 | 1,2e-4 / 6,0e-4 | 1,8e-7 / 2,9e-5 | 1,1e-14 / 1,1e-13 | 5e-9 | 1,2e-4 / 7,3e-7 | 2,4e11 |
| 1e10 | 8,3 | 2,1e-12 | 1,0 / 17 | 1,0 / 17 | 2,2e-9 / 1,9e-9 | 1,2e-4 / 5,6e-5 | 1,8e-5 / 4,9e-5 | 2,1e-10 / 5,1e-10 | 5e-10 | 1,7e-4 / 8,7e-7 | 1,4e16 |

(E2 coincidiu com E0 em todos os casos, a 2 algarismos: a aritmética do LS do solver não contribui.)

Leitura:

1. **E1 ≈ E0 ≈ E2.** O LS no w\* arredondado corretamente, feito em 80 dígitos, erra tanto quanto o solver. O erro é
   do estimador LS num ponto double, não da aritmética nem do ponto devolvido. Isso confirma o diagnóstico da hipótese.
2. **E3 ≪ E0.** Um único sistema de Newton-KKT no ponto final, todo em double, recupera λ a 1e-15..1e-13 (1e4),
   ~1e-12 (1e6), 1e-12..2e-10 (1e8) e 1e-10..1e-8 (1e10), com w intocado. Ganho de 4 a 9 ordens. Não é preciso
   avaliar nada em w_ext para isso.
3. **A Hessiana é o próximo limite.** Com o erro de H modelado em ε·cond(Σ), E3 piora para ~1e-10 a 1e8 e para
   1e-5..1e-4 a 1e10. O E5 (dd + w_lo guardado) contrai devagar com essa H (1e10: 1e-4 → 7e-7 em 3 passos).
   Com H exata, o resíduo dd leva a 1e-14..1e-17 até 1e8 (E3x).
4. **E4 > E3x a 1e10.** Reestimar por LS em w_ext reintroduz H·(erro de w_ext); guardar dλ é melhor que reestimar.
5. A segunda solução de Newton no mesmo w double (sem guardar dw) **não** melhora λ (repetida, deu o mesmo E3-H): para
   ir além é preciso guardar w_lo (E5), como em Demmel et al. 2006 (§1, p. 328: x em "doubled working precision").

Ressalvas: a matriz KKT de E3 é a Jacobiana de 80 dígitos arredondada (não a Hessiana analítica de `objective.cpp`), o
resíduo é um modelo (não o código de `kkt.cpp`), a referência usa b = 1/n exato enquanto a biblioteca usa o double
mais próximo (n = 8 tem 1/8 exato e mostra os mesmos níveis, então essa diferença não domina). O modelo de H é uma
hipótese calibrada por uma medida de `validation.md`; a primeira tarefa de um plano é medir com a H real (seção 6).

---

## 1. Sensibilidade dos multiplicadores num ponto KKT

### O que as fontes dizem

- **Teoria de sensibilidade (Fiacco 1976/1983; Robinson 1980; Bonnans & Shapiro 2000, cap. 4–5).** Com LICQ,
  condição suficiente de segunda ordem e complementaridade estrita, o par (x(p), λ(p)) é diferenciável no parâmetro p e
  sua derivada resolve o sistema linear com a matriz KKT K = [H −A'; A 0]: d(x, λ)/dp = −K⁻¹ ∂F/∂p. **Não verificado**
  (textos não acessíveis; enunciado de memória, padrão). Nocedal & Wright (2006), cap. 12 e 18, discutem o mesmo
  sistema e o multiplicador de mínimos quadrados λ = (AA')⁻¹A∇f: **não verificado** (o PDF do cap. 18 usado na nota
  anterior não está mais no ar; o cap. 17 foi verificado na nota anterior e não trata disto).
- **Gill & Murray (1979)**, "The computation of Lagrange-multiplier estimates for constrained minimization", Math.
  Programming 17:32–60: segundo o resumo indexado, separa o erro de uma estimativa em erro de decisão (conjunto ativo),
  de truncamento (x ≠ x\*) e computacional. **Não verificado** no texto (paywall). A distinção "estimativa de primeira
  ordem (LS), erro O(‖x − x\*‖), × estimativa de segunda ordem (com a Hessiana), erro O(‖x − x\*‖²)" é de Gill, Murray &
  Wright, *Practical Optimization* (1981): **de memória, não verificado**. Se correta, é exatamente o fenômeno aqui.

### Derivação própria: o erro do LS em w double

Seja F(w, λ) = (∇f(w) − A(w)'λ, h(w)) sobre as variáveis livres, com A = [1'; μ'; (Σw/vol)'] (linhas presentes) e
(w\*, λ\*) a solução. Para um w = w\* + δw, o LS de `estimate_multipliers` dá λ_LS(w) = (A')⁺∇f(w). Como
∇f(w) − A(w)'λ\* = H_L δw + O(‖δw‖²), com H_L = ∇²f − λ_vol∇²vol:

  **λ_LS(w) − λ\* = (A')⁺ H_L δw + O(‖δw‖²)**,  logo ‖δλ‖ ≲ ‖(A')⁺‖·‖H_L‖·‖δw‖.

Para um w double genérico ‖δw‖ ≥ ~u‖w‖ (u = 2⁻⁵³), e o erro relativo fica ≈ u·ρ·cond(A), com
ρ = (\|H_L\|\|w\|)_j / (maior termo de ∇f): é a razão que `validation.md` 731–743 já mede (u·ρ = 1,6e-7 a 1e6 e 8,4e-3 a
1e8, n = 20, sem box). O observado (n = 5, 8) é 3,5e-8..2,1e-7 a 1e6 e 3e-5..3,6e-3 a 1e8: **mesma ordem**, e a
escala ~cond(Σ)^2 vem de 2J'J (J ~ Σ/v). A cota não é inferior rigorosa (a projeção pode anular-se), mas E1 mostra que
o w\* arredondado corretamente já a atinge.

### Derivação própria: o erro do λ de Newton

O passo de Newton em (w, λ) resolve K(dw, dλ) = −F. Se F tem erro de avaliação e (e não o "erro" H δw, que o sistema
atribui a dw), a parte em λ do erro é, com Z base de null(A) e Y = range(A'):

  H dw − A'dλ = −e, A dw = 0 ⇒ dw = −Z(Z'HZ)⁻¹Z'e,  **dλ_err = (A')⁺ (I − HZ(Z'HZ)⁻¹Z') e**,

e o operador entre parênteses tem imagem em range(A') (Z' o anula). Logo ‖dλ_err‖ ≲ ‖(A')⁺‖(1 + ‖H‖‖(Z'HZ)⁻¹‖)‖e‖.
A diferença com o LS está no tamanho do vetor amplificado: **u·(termos de ∇f)** no Newton (o Σw é compensado; `validation.md`
703–705 mede o resíduo a 1,4e-4 do seu ruído) contra **u·\|H\|\|w\|** no LS — a razão ρ, 7e8 a 1e6 e 4e13 a 1e8. O
fator ‖H‖‖(Z'HZ)⁻¹‖ é o mesmo condicionamento reduzido que limita o erro progressivo dos pesos (`validation.md` 777–837),
o que explica por que λ via Newton fica na ordem do erro progressivo de w (E3 vs coluna "erro w").

### Veredito

**Confirma** que o erro é intrínseco ao **LS num w double**: nenhum w double dá λ_LS melhor que u·ρ·cond(A) (E1).
**Contradiz** a extensão "o λ é mal condicionado em relação a w": o sistema KKT conjunto separa δw de δλ, e o λ do
Newton no mesmo w double tem a precisão do problema (E3). Um w_ext só é necessário para passar do limite da Hessiana.

---

## 2. Refinamento iterativo com resíduo estendido

### O que as fontes dizem

- **Demmel, Hida, Kahan, Li, Mukherjee, Riedy (2006)**, ACM TOMS 32(2):325–351 (LAWN 165). **Verificado** pp. 325–330:
  - Três precisões (§2, p. 329): ε_w (dados e fatoração), ε_x (solução; "possibly twice the working precision … if
    necessary for componentwise convergence") e ε_r ≤ ε_w² (resíduo). Erros modelados por eqs. (3)–(5).
  - Clássico (§1, p. 326): com resíduo em ~dobro da precisão, x converge a ~ε_w "as long as the condition number of A is
    not too large (sufficiently less than 1/ε_w)", a custo O(n²) sobre O(n³).
  - Convergência: ‖dx^(i+1)‖ cai por ρ = O(ε_w)‖A‖‖A⁻¹‖ por passo (p. 329); parada por ‖dx‖/‖x‖ ≤ ε_w, razão de correções
    ≥ ρ_thresh (= 1/2, Wilkinson) ou número máximo de passos; cota de erro eq. (6), B_norm = max{(‖dx‖/‖x‖)/(1 − ρ_max), γε_w},
    γ = max(10, √n) (p. 330).
  - Resultados (p. 328): para κ < 1/(γε_w) a cota é ≤ 2γε_w e limita o erro verdadeiro; para κ entre 1/(γε_w) e ε_w⁻²,
    erros pequenos em 96% (normwise) / 94% (componentwise) dos casos, mediana de 2 e no máximo 4 passos; o mau
    condicionamento extremo é detectado por estimadores de condição. Guardar x em dobro (Algoritmo 2) é o que melhora
    os casos mais mal condicionados.
- **Tisseur (2001)**, SIAM J. Matrix Anal. Appl. 22(4):1038–1057 (Newton em ponto flutuante). **Verificado**:
  eq. (2.3) (p. 1040): erro do resíduo ‖e‖ ≤ u‖F(v)‖ + ψ(F, v, u, ū) com ū a precisão do resíduo; **Teorema 2.2**
  (p. 1041): ‖v̂ − v\*‖ ≤ G‖v − v\*‖ + g, g ≈ ‖J⁻¹‖ψ + u‖v\*‖, e G depende de ‖J⁻¹E‖ (Jacobiana inexata e solver);
  p. 1043: "the limiting accuracy is essentially independent of the solver", enquanto a **taxa** depende da precisão da
  Jacobiana; **Corolário 2.3** (p. 1043): com u·κ(J) ≤ 1/8 e u‖J⁻¹‖φ ≤ 1/8, o erro cai até ‖J⁻¹‖ψ/‖v\*‖ + u. Observação
  de Wilkinson citada na p. 1044: partir de um ponto melhor que a precisão limite e iterar "will usually spoil" o ponto.
- **Carson & Higham (2018)**, SIAM J. Sci. Comput. 40(2):A817–A847. **Verificado**: Tabela 5.1 (p. A828) — erro
  progressivo limite 4p·u_r·cond(A, x) + u, convergência sob 2u_s·min(cond(A), κ∞(A)μ_i) + u_s‖E_i‖ < 1; §7.1 (p. A829):
  u_f = u, u_r = u² dá erro limite ≤ (4φ/3 + 1)u com φ ≈ 3n·u_f·cond(A), "as long as φ … is sufficiently less than 1";
  §7.2: resíduo em precisão fixa dá só u·cond(A, x). Resumo (p. A817): com working precision single e LU em meia
  precisão, LU-IR resolve até κ∞ ≤ 1e4 e GMRES-IR (GMRES precondicionado pelos fatores LU) até κ∞ ≤ 1e8, "although in
  general there is no guarantee that GMRES will converge quickly" — o mecanismo conhecido para quando κ(K) passa de 1/u_f.
- Higham, *Accuracy and Stability of Numerical Algorithms*, 2ª ed., cap. 12 (refinamento) e cap. 25 (Newton): **não lido**;
  o conteúdo relevante foi conferido nos três artigos acima, que o citam.

### Aplicação a este problema

- O "sistema" refinado é a equação não linear F(w, λ) = 0 e a Jacobiana é K; vale o Teorema 2.2 de Tisseur. A precisão
  limite em λ é ‖[K⁻¹]_λ‖·ψ. Com resíduo double e Σw compensado, ψ ≈ u·(termos de ∇f), e isso já dá E3. Com resíduo dd
  em w_ext, ψ ≈ u² ·(termos) + u·‖F‖, e o limite é ~u (E3x, H exata).
- **Regimes** (medido aqui, n = 5, 8): cond(K) escalada vai de 7e4 (1e4) a 1,4e16 (1e10, n = 8). A garantia de
  Demmel/Tisseur (u·κ ≲ 1/8) cobre até cond(Σ) ~ 1e6–1e8; acima disso o refinamento **ainda funcionou** com H exata (E3x a
  1e10: 2e-16..5e-10), como os 94–96% de Demmel et al. Não há garantia teórica ali.
- **O limite prático é a Hessiana** (o "E" de Tisseur): com H errada em η ≈ u·cond(Σ), cada passo contrai o erro de λ por
  ~η·‖K⁻¹‖‖H‖ (E5). A 1e10, η·cond(K) não é pequeno e a convergência é lenta (fator ~0,1–0,6 por passo). Nenhum método
  baseado em double consegue λ a 1e-12 em cond(Σ) = 1e10 **com a Hessiana atual** sem muitos passos; com uma H mais
  exata (termos que cancelam em Dot2) o E3x mostra que consegue.

### Recomendação

Fase 1 (double): λ_final = λ_LS + dλ do Newton no ponto final. Esperado: ~1e-12 até 1e6 e ~1e-10 a 1e8 (E3-H); a 1e10,
1e-5..1e-4 (melhor que hoje, O(1), mas não validável). Fase 2 (opcional): 2–3 passos de refinamento com resíduo dd e
w_lo guardado (E5), parados pela regra de Demmel (razão de correções ≥ 1/2 ou ‖dλ‖ ≤ u‖λ‖), que também dá uma cota
a posteriori; antes, medir o erro real de H e, se for o limitante, calcular em Dot2 os produtos da Hessiana que cancelam.

---

## 3. Aritmética double-double e compensada

### O que as fontes dizem

- **Ogita, Rump, Oishi (2005)**, SIAM J. Sci. Comput. 26(6):1955–1988. **Verificado** (preprint do autor): TwoSum
  (Alg. 3.1, de Knuth, 6 flops, Teorema 3.4: a + b = x + y exato), Split e TwoProduct (Algs. 3.2–3.3, de Dekker e
  Veltkamp, 17 flops), **TwoProductFMA** (Alg. 3.5, 2 flops: y = FMA(a, b, −x)); Sum2 (Alg. 4.4, 7(n−1) flops,
  Prop. 4.5); **Dot2** (Alg. 5.3, 25n − 7 flops com TwoProduct; menos com FMA), **Proposição 5.5**:
  \|res − x'y\| ≤ eps\|x'y\| + γ_n²\|x'\|\|y\|, e Corolário 5.7: erro relativo ≤ eps + ½γ_n²·cond(x'y), "as if computed in
  twice the working precision"; DotK (Alg. 5.10) para K-vezes. O `compensated_dot` de `objective.cpp` é o Dot2 com
  TwoProductFMA.
- **Hida, Li, Bailey**, "Library for Double-Double and Quad-Double Arithmetic" (doc do QD, fonte `qd.tex` da
  distribuição Debian 2.3.7). **Verificado**: dd = soma não avaliada de dois doubles (≥ 106 bits); Two-Prod-FMA
  descrito nas preliminares; exige aritmética IEEE com arredondamento ao par e, em x87, "the round-to-double flag must be
  enabled" (`fpu_fix_start`); "the basic routines do not have a full correctness proof"; **não dá** cota de erro nem
  contagem de flops para dd (só tempos de quad-double). Custos de dd citados de cabeça (soma ~20 flops, produto ~10 com
  FMA): **não verificados**.
- Dekker (1971), Numer. Math. 18:224–242; Knuth, TAOCP vol. 2; Shewchuk (1997), Discrete Comput. Geom. 18:305–363:
  **não lidos**; os algoritmos usados aqui foram conferidos em Ogita-Rump-Oishi, que os reproduz.

### Como calcular ∇f(w_ext) e o resíduo KKT só com double + fma (derivação própria)

A amplificação está toda nos produtos por Σ de vetores que contêm w: δ(Σw) = Σδw ~ u\|Σ\|\|w\| ≫ u\|Σw\| perto da carteira
de mínima variância (`objective.cpp`, comentário de `covariance_times`). Com w_ext = w_hi + w_lo (\|w_lo\| ≲ u\|w_hi\|):

1. p = Σw_ext: Dot2 sobre os 2n termos (Σ_ij, w_hi,j) e (Σ_ij, w_lo,j), resultado arredondado a double — erro ≈ u\|p\|.
2. v = w_ext'p e 1'w_ext − 1, μ'w_ext − r: Dot2 sobre 2n termos (o resíduo das igualdades precisa de w_lo).
3. c_i = w_i p_i/v, r = c − b, J: em double com w_hi (erro relativo u por termo, nível do "ruído de avaliação").
4. g = 2J'r: os termos Σ(w∘r) (já em `numerator_gradient`) em Dot2 sobre o par; o resto em double.
5. Estacionaridade g − A'λ: em double (λ e A com erro relativo u por termo).

Isso é o "resíduo dd" de E3x no sentido de Tisseur (ψ ≈ u·termos + u²·\|H\|\|w\|). Não é preciso um tipo dd completo:
basta um `covariance_times` e um `compensated_dot` que aceitem o par (hi, lo). Portabilidade: só `std::fma` e somas
na ordem escrita (o repositório já compila os testes com `-ffp-contract=off`, commit `5882596`); em MSVC, `std::fma` é
exato por norma, mas lento sem FMA de hardware, e a garantia de não-contração sob `/fp:precise` **não foi verificada**.

### Custo (derivação própria, n ≤ 100)

O polish faz por passo: `evaluate_kkt` com J'J (n³ flops) e o termo de segunda ordem (n²), e um QR de (n_F + m)²
(~(4/3)(n + 3)³). Um resíduo dd custa ~3 produtos Σ·(par) a ~11 flops por termo ≈ 66n² ≈ 6,6e5 flops a n = 100, contra
~2,4e6 do passo: < 30% de um passo, e dispensa recalcular H (a matriz e sua fatoração do último passo servem, como no
refinamento de Demmel/OSQP). Fase 1 custa um sistema KKT a mais (às vezes zero: ver seção 6).

### Quantos passos

Com H exata, 1 passo dd bastou em todos os casos (E3x). Com H a u·cond(Σ): 1–2 passos até 1e6, 3 passos a 1e8 para
~1e-11..1e-10, e a 1e10 o terceiro passo ainda melhora (7e-7). Regra de parada: a de Demmel et al. (p. 330).

---

## 4. Alternativas e software comparável

- **OSQP** (Stellato, Banjac, Goulart, Bemporad, Boyd 2020, Math. Prog. Comp. 12:637–672; arXiv:1711.08013v4).
  **Verificado** §4 "Solution polishing" (pp. 10–11): com o conjunto ativo adivinhado pelos sinais de y, resolve **um**
  sistema KKT para (x, y) juntos, regularizado por δ ≈ 1e-6, e corrige a regularização com refinamento iterativo
  (eq. (31)) "for a fixed number of passes" (padrão 3, §7 "OSQP"; versão arXiv mais recente, não a publicada). É o análogo direto de "guardar dλ": o multiplicador sai do
  sistema conjunto, não de um LS no x final.
- **IPOPT** (documentação de opções, coin-or.org/Ipopt/OPTIONS.html). **Verificado**: `recalc_y` ("recalculate the
  equality and inequality multipliers as least square estimates", desligado por padrão, custa uma fatoração extra);
  `least_square_init_duals`, `constr_mult_init_max` (LS só para iniciar); refinamento iterativo "on the full
  unsymmetric system" (`min_refinement_steps` = 1, `max_refinement_steps` = 10, `residual_ratio_max` = 1e-10). Nenhuma
  opção controla ou relata a precisão dos multiplicadores finais além de `dual_inf_tol`.
- **LAPACK xGERFSX** (Demmel et al. 2006): devolve cotas de erro normwise e componentwise por solução. Um análogo aqui
  seria reportar uma cota de λ pela eq. (6) (razão de correções sucessivas). Só faz sentido com a Fase 2 (a Fase 1 dá
  uma correção só, sem razão).
- **quadprog / KNITRO**: não consultados nesta nota.
- Literatura específica sobre multiplicadores precisos em problemas mal condicionados: Gill & Murray 1979 (acima,
  não verificado). Nada mais encontrado nas buscas feitas (o que não prova que não exista).

Opções para a biblioteca:

| Opção | Precisão de λ (medida aqui, n = 5, 8) | Custo | Risco |
|---|---|---|---|
| (a) Documentar o limite (hoje) | 5e-10 / 3e-7 / 4e-2 / O(1) a 1e4/1e6/1e8/1e10 | 0 | nenhum |
| (b) Fase 1: λ_LS + dλ final, double | ~1e-14 / ~1e-12 / ~1e-10 / ~1e-4 | 1 solve KKT | o motivo do LS (s ≈ σ_min) — ver seção 6 |
| (c) Fase 2: + resíduo dd, w_lo, 2–3 passos | ~u / ~u / 1e-13..1e-10 / ~1e-6 (com a H atual) | (b) + ~1 passo | código novo pequeno; H limita a 1e10 |
| (d) (c) + H com somas compensadas | ~u até 1e10 (E3x) | (c) + O(n³) dd | não medido com a H real |
| (e) Cota certificada (eq. (6) de Demmel) | — | sobre (c) | API nova; cota só heurística acima de 1/(γu) |

Os pesos continuam double em todas: w_lo é interno; o resultado reporta w_hi.

---

## 5. Caminho QP (Goldfarb-Idnani): há análogo?

Derivação própria: no QP, ∇f = Σw (ou Σw − tμ), H = Σ, e o erro de um λ_LS em w arredondado seria
(A')⁺Σδw ~ u\|Σ\|\|w\|/\|Σw\|, no máximo ~u·cond(Σ) (não ~u·cond(Σ)²·³). Além disso, o GI não reestima λ em x: os
multiplicadores saem da fatoração do conjunto ativo junto com x (`validation.md` 285–288 descreve a correção conjunta de
x e multiplicadores), que é o equivalente de "guardar dλ". O que `validation.md` mede no QP: estacionaridade relativa
≤ 1e-7 a 1e8 (fuzz), formas fechadas de multiplicador a 1e-12, multiplicadores de `MAX_RETURN`/leverage contra
diferenças finitas (1e-5/1e-6, limitados pela diferença), e inclinações de multiplicadores por diferença de soluções
com arredondamento "1e6 eps of the multiplier scale at … 1e8" (linhas 208–213), contornado por LS sobre as normais.
**Não há** medida do erro progressivo dos multiplicadores do QP contra uma referência estendida a 1e8; a estrutura
sugere ≤ u·cond(Σ). Recomendação: uma medida com a mesma referência quando o plano rodar; nenhuma mudança esperada.

---

## 6. O que um plano "accurate multipliers" deve conter (priorizado)

1. **[Primeiro, barato] Medir E3 com a Hessiana real.** Um teste C++ desabilitado (padrão `DISABLED_Measure…`) nos casos
   de `RiskParityFuzz.TargetVolatilityMultipliersByCondition` e nos do probe (sem box, s = 1,2σ_min, orçamento igual,
   n = 5 e 8, seeds 7/11/3, cond 1e4..1e10): no ponto devolvido, montar `newton_system` com a `evaluate_kkt` atual,
   resolver, e comparar λ_LS e λ_LS + dλ com o refinamento `long double` de `risk_parity_refinement.hpp` (e, para 1e10,
   com o script de 80 dígitos, que é independente). Critério: confirma se λ_LS + dλ ≤ ~1e-10 a 1e8 (prevê E3-H);
   refuta se ficar na ordem de E0. Medir também o erro da Hessiana analítica contra a de 80 dígitos (η real) e o
   resíduo de estacionaridade de (w, λ_LS + dλ) em unidades do piso (ver item 4).
2. **Precisão-alvo por faixa** (se o item 1 confirmar E3-H): ≤ 1e-12 relativo até cond(Σ) = 1e6; ≤ 1e-9 a 1e8 (com folga
   10 sobre o medido); a 1e10, Fase 1 dá ~1e-4: manter "não validado" ou aplicar a Fase 2. As faixas dos pesos não mudam.
3. **Implementação, Fase 1** (`src/riskparity/newton_polish.cpp`):
   - Após o laço de `newton_polish`, no `result.point` (multiplicadores já `settled` por LS), montar `newton_system`,
     resolver por QR e somar dλ às igualdades mantidas, **sem mudar w**. Quando o laço terminou por um passo recusado,
     esse sistema já foi resolvido em `newton_step` (a partir do ponto aceito): basta devolver o `step` e reutilizá-lo.
     Unidades: a coluna b da matriz é −scale_b·∇h_b, então dλ_b = scale_b·z_b, com scale_b = `hessian_scale` / max\|linha\|.
   - **Multiplicadores de bound** (fraca em held weights): β_j = (∇f − A'λ)_j em w sofre o mesmo H δw. Com o dw livre do
     sistema, usar β_j = (∇f − A'λ_new)_j + (H_L dw)_j (linha j de H_L sobre as livres); verificar sinal como em `settled`.
   - Salvaguarda (o motivo do LS em `validation.md` 883–889, s = σ_min(1 + 1e-9)): aceitar dλ só se o QR tem posto
     completo, ‖dλ‖ ≤ c‖λ‖ (c pequeno, p. ex. 1e-2) e uma segunda solução a partir de (w, λ_new) devolve ‖dλ₂‖ ≤ ‖dλ‖/2;
     senão, manter λ_LS. A aceitação do passo de w continua como hoje (o dλ final não entra na decisão de passos).
4. **Portão KKT**: `is_kkt_point` e `expect_risk_parity_kkt` avaliam a estacionaridade em w com os multiplicadores do
   resultado. Com λ\*, o resíduo em w é H δw inteiro (inclui a componente em range(A') que o LS absorvia), que pode
   exceder 2× o piso onde o erro progressivo de w passa de ~u (ex.: 1e8, n = 8, erro de w 3,6e-13). Medir no item 1;
   se exceder, aplicar o portão com λ_LS (o que já é feito hoje) e reportar λ_new, e no teste checar a estacionaridade
   com os multiplicadores LS recalculados pelo próprio teste e os multiplicadores do resultado **contra a referência**.
5. **Fase 2 (só se o alvo a 1e8/1e10 exigir)**: `covariance_times`/`compensated_dot` sobre o par (hi, lo)
   (`objective.cpp`), resíduo KKT em w_ext (`kkt.cpp`), 2–3 passos reutilizando a fatoração do último passo, parada de
   Demmel (‖dλ‖ ≤ u‖λ‖ ou razão ≥ 1/2); antes, se o item 1 mostrar η ≈ u·cond(Σ), calcular em Dot2 os termos de H que
   cancelam. Sem `long double` (MSVC).
6. **Testes que fixam o resultado**:
   - `RiskParityFuzz.TargetVolatilityMultipliersByCondition`: novos tetos medidos (hoje 1e-9/1e-6/1e-1).
   - Novo: o mesmo com box `LongOnly`/`Capped` e vol-alvo (multiplicadores de bound com o termo H dw).
   - Regressão: s = σ_min(1 + δ), δ ∈ {1e-2, …, 1e-10} (tabela de `validation.md` 901–909) continua `OPTIMAL` e passa
     `expect_risk_parity_kkt`; registrar quantas vezes a salvaguarda manteve λ_LS.
   - `RiskParity.MultipliersAreDerivatives` estendido a cond 1e6 (diferenças finitas limitadas a ~1e-6, então só pega
     erros grosseiros).
   - `validation.md`: substituir a tabela de 839–861 e o "Future work"; `portopt.h` 124–127: novas faixas.
   - Caminho QP: uma medida dos multiplicadores do GI contra a referência estendida a 1e8 (item 5 desta nota).

---

## Bibliografia

| Referência | URL | Status |
|---|---|---|
| Demmel, J.; Hida, Y.; Kahan, W.; Li, X. S.; Mukherjee, S.; Riedy, E. J. "Error bounds from extra-precise iterative refinement", ACM TOMS 32(2):325–351, 2006 (LAWN 165) | https://people.eecs.berkeley.edu/~wkahan/p325-demmel.pdf | Verificado: pp. 325–330 (resumo, §1, §2 eqs. (3)–(6), critérios de parada). Seções 3–6 não lidas. |
| Tisseur, F. "Newton's method in floating point arithmetic and iterative refinement of generalized eigenvalue problems", SIAM J. Matrix Anal. Appl. 22(4):1038–1057, 2001 | https://eprints.maths.manchester.ac.uk/979/1/tisseur3.pdf | Verificado: eqs. (2.2)–(2.8), Teorema 2.2, Corolário 2.3, pp. 1040–1044. |
| Carson, E.; Higham, N. J. "Accelerating the solution of linear systems by iterative refinement in three precisions", SIAM J. Sci. Comput. 40(2):A817–A847, 2018 | https://eprints.maths.manchester.ac.uk/2629/1/cahi18.pdf | Verificado: resumo, Tabela 1.1, Tabela 5.1, §7.1–7.2 (pp. A817–A829). |
| Ogita, T.; Rump, S. M.; Oishi, S. "Accurate sum and dot product", SIAM J. Sci. Comput. 26(6):1955–1988, 2005 | https://www.tuhh.de/ti3/paper/rump/OgRuOi05.pdf | Verificado (preprint): Algs. 3.1–3.5, Teorema 3.4, Alg. 4.4, Prop. 4.5, Alg. 5.3, Prop. 5.5, Cor. 5.7. |
| Hida, Y.; Li, X. S.; Bailey, D. H. "Library for Double-Double and Quad-Double Arithmetic" (doc. do QD) | https://sources.debian.org/data/main/q/qd/2.3.7-2.1/docs/qd.tex | Verificado: representação, Two-Prod-FMA, nota x86, ausência de cotas dd. Custos por operação dd não verificados. |
| Stellato, B.; Banjac, G.; Goulart, P.; Bemporad, A.; Boyd, S. "OSQP: an operator splitting solver for quadratic programs", Math. Prog. Comp. 12:637–672, 2020 | https://arxiv.org/pdf/1711.08013 | Verificado (arXiv): §4, eqs. (27)–(31), padrão de 3 passos de refinamento. |
| IPOPT, documentação de opções | https://www.coin-or.org/Ipopt/OPTIONS.html | Verificado: `recalc_y`, `least_square_init_duals`, `constr_mult_init_max`, opções de refinamento, `dual_inf_tol`. |
| Gill, P. E.; Murray, W. "The computation of Lagrange-multiplier estimates for constrained minimization", Math. Programming 17:32–60, 1979 | https://doi.org/10.1007/BF01588224 | Bibliografia e resumo indexado (erro de decisão/truncamento/computacional); texto **não lido**. |
| Gill, P. E.; Murray, W.; Wright, M. H. *Practical Optimization*, Academic Press, 1981 | — | **Não lido**; estimativas de 1ª/2ª ordem de memória. |
| Fiacco, A. V. *Introduction to Sensitivity and Stability Analysis in Nonlinear Programming*, Academic Press, 1983; Robinson, S. M. "Strongly regular generalized equations", Math. Oper. Res. 5(1), 1980; Bonnans, J. F.; Shapiro, A. *Perturbation Analysis of Optimization Problems*, Springer, 2000 | — | **Não lidos**; teorema de sensibilidade enunciado de memória. |
| Nocedal, J.; Wright, S. J. *Numerical Optimization*, 2ª ed., 2006, caps. 12 e 18 | — | **Não lidos** nesta nota (o PDF do cap. 18 não está mais acessível). |
| Higham, N. J. *Accuracy and Stability of Numerical Algorithms*, 2ª ed., SIAM, 2002, caps. 12 e 25 | — | **Não lido**; conteúdo conferido via Tisseur e Carson-Higham. |
| Dekker, T. J. "A floating-point technique for extending the available precision", Numer. Math. 18:224–242, 1971; Knuth, D. E. TAOCP vol. 2 | — | **Não lidos**; algoritmos conferidos em Ogita-Rump-Oishi. |
| Shewchuk, J. R. "Adaptive precision floating-point arithmetic and fast robust geometric predicates", Discrete Comput. Geom. 18:305–363, 1997 | — | **Não lido**; não usado nas conclusões. |
| Moler, C. B. "Iterative refinement in floating point", J. ACM 14(2), 1967; Wilkinson, J. H. *Rounding Errors in Algebraic Processes*, 1963 | — | **Não lidos**; citados via Demmel et al. e Tisseur (p. 1044). |
