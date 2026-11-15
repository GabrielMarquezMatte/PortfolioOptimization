# PortOpt — documento de design

Data: 2026-11-15
Status: revisado (builder check); aguardando aprovação do autor

## 1. Intenção

PortOpt é uma biblioteca nativa de otimização de portfólio, chamável de qualquer linguagem por uma interface C
estável. As equações vêm do pacote R do autor (`R/PortfolioOptimization.R` no repositório antigo
`GabrielMarquezMatte/PortfolioOptimization`, cujo remote será apagado e recriado para este projeto); as derivações
são do próprio autor.

O problema que ela resolve não é velocidade: é portabilidade. Em R as funções ficam presas ao R. Aqui o núcleo é
C++23, o único contrato é a interface C, e Python (primeiro) e C++ (wrapper header-only) consomem esse contrato;
outras linguagens virão depois sobre o mesmo header.

A prioridade do autor é **uma API sucinta e fácil de usar, trabalhada agora para não mudar depois**. Três
conceitos: **spec** (o problema: objetivo e restrições), **dados** (as entradas numéricas) e **solve** (resolve e
devolve um resultado).

O projeto segue a disciplina do VolEngine (`C:\Users\x\Documents\Programs\Cpp\VolEngine`): ARCHITECTURE.md,
STYLEGUIDE.md, CONTRIBUTING.md, CHANGELOG.md, specs e plans em `docs/superpowers/`, `.clang-format`,
`.clang-tidy`, `tests/style/check_style.py`, regras de estabilidade da ABI. Onde este documento não diz nada, vale
o que o VolEngine faz.

## 2. Objetivos, não-objetivos e anti-objetivos

### Objetivos

1. Resolver os problemas da seção 3 a partir de μ e Σ (e, no tracking error, de uma matriz de retornos e um
   benchmark) fornecidos pelo usuário, devolvendo pesos, status e multiplicadores.
2. Métricas de carteira (retorno, volatilidade, contribuições de risco, Sharpe, razão de diversificação, utilidade,
   tracking error) como funções sobre pesos e dados.
3. Interface C estável a partir da 1.0, consumida por Python (ctypes + numpy) e por um wrapper C++23 header-only.
4. Resultados corretos com precisão de máquina onde a matemática permite, verificados por condições de KKT e
   formas fechadas (seção 9).

### Não-objetivos da v1

- Algoritmo genético e o `momentum_utility_portfolio` (decisão do autor sobre o GA; o objetivo sai junto, ver
  seção 13).
- `pesos_carteira` (pesos aleatórios), `merton_portfolio` (simulação de consumo) e o wrapper `portfolios`: não são
  otimizações.
- Estimar μ ou Σ, calcular retornos a partir de preços, anualizar qualquer coisa. Os dados chegam prontos.
- Restrições lineares de grupo (setor), turnover, custos de transação, cardinalidade. Restrições novas entram
  depois como setters novos, sem quebra.
- Lote com threads (entra depois, sem quebra; ver seção 6.4).
- Precisão `float`.

### Anti-objetivos

- **Quebrar a interface C depois da 1.0.** Nenhuma assinatura muda, nenhum valor de enum é removido ou renumerado,
  nenhuma semântica documentada muda.
- **API inchada.** Antes de qualquer função nova: "o usuário consegue fazer isso com o que já existe?" Se sim, não
  entra.
- **Gambiarras numéricas no contrato.** Nada de variância `eps` para o ativo livre de risco, nada de grid search
  onde existe formulação exata, nada de bounds ignorados em silêncio.
- **Falhar em silêncio.** Problema inviável ou combinação sem sentido devolve status ou erro explícito, nunca NULL
  silencioso como o `purrr::possibly` do R.

## 3. Problemas (v1)

### 3.1 Notação

- N ativos de risco; w ∈ ℝᴺ os pesos dos ativos de risco; μ ∈ ℝᴺ; Σ ∈ ℝᴺˣᴺ simétrica definida positiva;
  σᵢ = √Σᵢᵢ.
- **Caixa (ativo livre de risco)**: quando os dados têm `risk_free(rate)`, existe o ativo N+1 com retorno rf,
  variância e covariâncias zero. Seu peso é w₀ = 1 − 1'w. Retorno da carteira: rf + (μ − rf·1)'w. Variância:
  w'Σw. Σ continua só dos ativos de risco e invertível; nenhum `eps`. Sem caixa, vale o orçamento 1'w = 1.
- Os pesos devolvidos têm N entradas sem caixa e N+1 com caixa (caixa por último).
- Unidades: todas as entradas estão nas mesmas unidades do usuário (μ, Σ, rf, linhas da matriz de retornos). A
  biblioteca não converte nada.

### 3.2 Restrições (setters da spec)

| Setter | Significado | Padrão |
|---|---|---|
| `bounds(lo, hi)` | lo ≤ wᵢ ≤ hi para todo ativo, caixa incluído | sem bounds: (−∞, +∞) |
| `asset_bounds(lo[], hi[])` | por ativo; vetores de tamanho N, ou N+1 com caixa (validado no solve) | — |
| `leverage(L)` | Σ\|wᵢ\| ≤ L sobre os ativos de risco (o caixa fica fora; o endividamento é controlado pelos bounds do caixa) | sem limite |
| `target_return(r)` | retorno da carteira = r (igualdade) | — |
| `max_volatility(s)` | √(w'Σw) ≤ s | — |
| `target_volatility(s)` | √(w'Σw) = s (igualdade) | — |
| `risk_budget(b[])` | contribuições de risco-alvo, bᵢ > 0, Σbᵢ = 1 | 1/N |
| `risk_aversion(δ)` | δ > 0 da utilidade | 5 (padrão do R) |
| `importance(λ)` | λ ∈ [0, 1) | 0.5 (padrão do R) |

Não existe flag `short`: venda a descoberto é bounds com mínimo negativo, e o padrão já é sem bounds.

### 3.3 Objetivos

Cada objetivo é um valor de enum. A matriz abaixo diz quais setters cada um aceita; qualquer outra combinação é
`INVALID_ARGUMENT` no solve, com mensagem que nomeia o setter.

| Objetivo | Problema | μ | bounds | leverage | target_return | vol | caixa | parâmetro |
|---|---|---|---|---|---|---|---|---|
| `MIN_VARIANCE` | min w'Σw | se target_return | sim | sim | sim (= fronteira) | — | sim | — |
| `MAX_RETURN` | max retorno s.a. vol ≤ s | sim | sim | sim | — | `max_volatility` obrigatório | sim | — |
| `MAX_SHARPE` | max (retorno − rf)/vol | sim | sim | sim | — | — | peso 0 (ver abaixo) | — |
| `UTILITY` | max retorno − (δ/2) w'Σw | sim | sim | sim | — | — | sim | δ |
| `IMPORTANCE` | max λ·retorno − (1−λ) w'Σw | sim | sim | sim | — | — | sim | λ |
| `KELLY` | ver abaixo | sim | sim | sim | — | — | sim | — |
| `TRACKING_ERROR` | min ‖Rw − b‖² | se target_return | sim | sim | sim | — | sim (coluna constante rf) | — |
| `MAX_DIVERSIFICATION` | max (w'σ)/√(w'Σw) | não | sim | sim | — | — | peso 0 | — |
| `RISK_PARITY` | min Σ(cᵢ − bᵢ)², cᵢ = wᵢ(Σw)ᵢ / w'Σw | se target_return | sim | — | sim | `target_volatility` | — | b |

Notas por objetivo:

- **Fronteira** não é objetivo: é `MIN_VARIANCE` com `target_return`, resolvido uma vez por ponto.
- **Vol-alvo** (o `port_vol_target` do R) é `MAX_RETURN` com `max_volatility`. Se s está abaixo da volatilidade
  da carteira de mínima variância (com as mesmas restrições), o status é `INFEASIBLE`; o R devolvia a mínima
  variância em silêncio. Se bounds ou leverage tornam o conjunto viável compacto e a carteira de retorno máximo tem
  volatilidade abaixo de s, a restrição fica inativa: o resultado é essa carteira, `OPTIMAL`, com multiplicador de
  volatilidade zero.
- **MAX_SHARPE** com caixa: a tangente é por definição só de ativos de risco, então o caixa tem peso fixo 0. Se os
  bounds do caixa excluem 0, o status é `INFEASIBLE` (nada é ignorado em silêncio). O rf dos dados é o limiar; sem
  `risk_free`, o limiar é 0. Se nenhuma carteira viável tem excesso de retorno positivo, o status é `INFEASIBLE`.
  Se existe excesso positivo mas o supremo não é atingido (a homogeneização só chega a κ = 1'y → 0, por exemplo sem
  bounds e com 1'Σ⁻¹(μ − rf) ≤ 0), o status é `UNBOUNDED` e os pesos são NaN.
- **MAX_DIVERSIFICATION** com caixa: a razão não depende do caixa, que tem peso fixo 0, com a mesma regra de bounds
  do Sharpe. O caso de supremo não atingido (com venda a descoberto e 1'Σ⁻¹σ ≤ 0) também é `UNBOUNDED`.
- **KELLY** (formulação do R, `kelly_port`): max (μ − rf)'x − ½ x'Σx s.a. 1'x = 1/(1+rf), w = (1+rf)x, bounds
  divididos por (1+rf). Substituindo x = w/(1+rf), é exatamente `UTILITY` com δ = 1/(1+rf) sobre μ (δ = 1 sem
  caixa) e é implementado como essa redução; bounds e leverage valem sobre w. Com caixa (equivale ao
  `use_rf = TRUE` do R), sem restrições os ativos de risco ficam w = (1+rf)Σ⁻¹(μ − rf) e o caixa com o resto. Sem
  `risk_free`, rf = 0.
- **IMPORTANCE** é exatamente `UTILITY` com δ = 2(1−λ)/λ, e com λ = 0 é a mínima variância; implementado como essa
  redução. KELLY e IMPORTANCE continuam como valores de enum porque são o vocabulário do autor (decisão do autor,
  2026-11-15); internamente não há solver próprio para eles.
- **TRACKING_ERROR** recebe R (T×N) e b (T) nos dados. Sem anualização: o R antigo dividia rf por 252 e anualizava
  μ e Σ; aqui tudo está nas unidades do usuário. Com caixa, o caixa é eliminado por w₀ = 1 − 1'w, e o resíduo fica
  (R − rf·1·1')w + rf·1 − b sobre N variáveis (uma coluna constante rf deixaria R'R singular com rf = 0). A matriz
  que precisa ter posto coluna completo é a dos retornos em excesso, R − rf·1·1' (ou R sem caixa); se não tiver,
  `INVALID_ARGUMENT`. Σ não é usada: `cov` pode ser NULL nos dados.
- **RISK_PARITY**: sem bounds, o problema tem uma solução por ortante; a biblioteca devolve a de pesos positivos,
  que é a que o R encontrava partindo de pesos iguais. Com `target_return`, `target_volatility` ou bounds ativos, o
  problema é o do `equal_risk_cont` do R: min Σ(cᵢ − bᵢ)² com essas restrições como igualdades (vol-alvo com
  1'w = 1, sem escalar a carteira). Ele é sobredeterminado e o resíduo ótimo em geral é maior que 0: é um
  compromisso, não um risk parity exato (decisão do autor, 2026-11-15). O status é `OPTIMAL` quando o solver
  converge para um ponto de KKT, qualquer que seja o resíduo, e `MAX_ITER` quando não converge; o resíduo
  Σ(cᵢ − bᵢ)² é lido do resultado.

## 4. Abordagem de solução

### 4.1 Abordagens consideradas

1. **Solvers próprios, sem dependência** (escolhida pelo autor): QP denso dual active-set de Goldfarb-Idnani (o
   algoritmo do `quadprog`, que gerou os resultados do R), Newton para risk parity, Lagrangiano aumentado para o
   risk parity com restrições extras. Controle total, multiplicadores exatos, wheels triviais; precedente do
   VolEngine (L-BFGS e álgebra linear próprios). Custo: escrever e testar o QP, com cuidado em degenerescência.
2. **HiGHS (MIT) + NLopt**: menos código numérico, mas dependência pesada no build e nos wheels, multiplicadores e
   tolerâncias nas convenções deles, e o Newton do risk parity continua necessário.
3. **OSQP (ADMM)**: integração simples, mas precisão de 1e-4 a 1e-6 inviabiliza a validação por KKT e formas
   fechadas. Descartada.

### 4.2 Como cada objetivo é resolvido

Quase tudo se reduz a um único QP estritamente convexo:

| Caso | Método |
|---|---|
| Sem bounds, sem leverage, só igualdades | forma fechada: sistema KKT por Cholesky de Σ |
| MIN_VARIANCE, UTILITY, IMPORTANCE, KELLY, TRACKING_ERROR com bounds | QP de Goldfarb-Idnani |
| MAX_SHARPE | homogeneização (Cornuejols-Tütüncü): min y'Σy s.a. (μ − rf)'y = 1, restrições homogeneizadas com κ = 1'y eliminado (lo·1'y ≤ yᵢ ≤ hi·1'y; Σ\|yᵢ\| ≤ L·1'y), 1'y > 0; w = y/1'y. Substitui o grid de 500 pontos do R |
| MAX_DIVERSIFICATION | homogeneização (Choueifaty): igual ao Sharpe com σ no lugar de μ − rf |
| MAX_RETURN com vol ≤ s | sequência de QPs de utilidade com busca em δ até √(w'Σw) = s (o builder escolhe bissecção ou outra raiz) |
| leverage | família implícita de restrições s'w ≤ L, s ∈ {±1}ᴺ. O Goldfarb-Idnani é dual e adiciona a restrição mais violada a cada passo; a mais violada é s = sign(w), com violação ‖w‖₁ − L. Exato e mantém a Hessiana definida positiva (dividir w = w⁺ − w⁻ deixaria a Hessiana singular, que o Goldfarb-Idnani não aceita). O builder pode escolher outro tratamento exato |
| RISK_PARITY sem target_return/target_volatility, bounds inativos | Newton sobre a forma convexa de Spinu: min ½ y'Σy − Σ bᵢ log yᵢ, y > 0; w = y/1'y |
| RISK_PARITY com target_return, target_volatility ou bounds ativos | Lagrangiano aumentado (com L-BFGS-B ou equivalente para os bounds) sobre Σ(cᵢ − bᵢ)², partindo da solução de Spinu |

### 4.3 Peças do núcleo

- `linalg`: Cholesky, solves triangulares, atualizações de posto usadas pelo Goldfarb-Idnani. Próprio, sem Eigen.
- `qp`: Goldfarb-Idnani denso, com igualdades, desigualdades explícitas e famílias implícitas de restrição (o
  oráculo de "mais violada" do leverage). Devolve x, multiplicadores, conjunto ativo, iterações e status.
- `problems`: monta, para cada objetivo, o QP (ou a sequência de QPs, ou o problema de risk parity) a partir de
  spec e dados, e traduz a solução de volta para pesos e multiplicadores nas variáveis do usuário.
- `riskparity`: Newton de Spinu e Lagrangiano aumentado.
- `metrics`: as funções da seção 5.3.
- `capi`: handles e fronteira de erro, como no VolEngine.

Os nomes de pastas são sugestão; o builder decide, respeitando o limite de 8 módulos por pasta.

## 5. Interface C

Prefixo `po_` / `PO_`, header `include/portopt/portopt.h`. Os nomes exatos ficam com o builder; a forma abaixo é
o contrato.

### 5.1 Regras de estabilidade

As do VolEngine (ARCHITECTURE.md, "Stability rules of the C API"), integralmente:

- Nenhuma struct pública; configuração por setters em handles opacos. Configuração nova é função nova, nunca
  parâmetro novo nem `_ex`.
- Enums são `typedef int32_t` + `#define`; valores só são acrescentados.
- Toda função falível devolve `po_status` e põe a mensagem em `po_last_error()` (thread-local); nenhuma exceção
  cruza a fronteira (`std::bad_alloc` vira `PO_OUT_OF_MEMORY`, o resto `PO_INTERNAL_ERROR`).
- Sem estado global mutável. Handles independentes: liberar um nunca invalida outro.
- Tipos na fronteira: `int32_t`, `size_t`, `double`, ponteiros.
- `PO_ABI_VERSION` no header e `po_abi_version()` em runtime; contrato congela na 1.0; 0.x pode quebrar.

### 5.2 Esboço

```c
typedef int32_t po_status;      /* PO_OK, PO_INVALID_ARGUMENT, PO_INVALID_HANDLE, PO_UNSUPPORTED,
                                   PO_OUT_OF_MEMORY, PO_INTERNAL_ERROR */
typedef struct po_spec po_spec;
typedef struct po_data po_data;
typedef struct po_result po_result;

/* objective: PO_OBJECTIVE_MIN_VARIANCE .. PO_OBJECTIVE_RISK_PARITY */
po_status po_spec_create(int32_t objective, po_spec** out);
po_status po_spec_set_bounds(po_spec* s, double lower, double upper);
po_status po_spec_set_asset_bounds(po_spec* s, const double* lower, const double* upper, size_t n);
po_status po_spec_set_leverage(po_spec* s, double max_gross);
po_status po_spec_set_target_return(po_spec* s, double r);
po_status po_spec_set_max_volatility(po_spec* s, double sigma);
po_status po_spec_set_target_volatility(po_spec* s, double sigma);
po_status po_spec_set_risk_budget(po_spec* s, const double* b, size_t n);
po_status po_spec_set_risk_aversion(po_spec* s, double delta);
po_status po_spec_set_importance(po_spec* s, double lambda);
void po_spec_free(po_spec* s);

/* Copies mu and cov (n x n, row-major, contiguous). Either may be NULL; the solve checks that the
   objective has what it needs (cov is only unused by TRACKING_ERROR). cov, when given, is checked here:
   finite, symmetric, positive definite. */
po_status po_data_create(size_t n, const double* mu, const double* cov, po_data** out);
po_status po_data_set_risk_free(po_data* d, double rate);          /* adds the cash asset */
po_status po_data_set_returns(po_data* d, const double* returns,   /* t x n, row-major; copied */
                              const double* benchmark, size_t t);
void po_data_free(po_data* d);

po_status po_solve(const po_spec* s, const po_data* d, po_result** out);

po_status po_result_status(const po_result* r, int32_t* out);
          /* PO_SOLVE_OPTIMAL, _INFEASIBLE, _UNBOUNDED (supremum not attained), _MAX_ITER */
po_status po_result_asset_count(const po_result* r, size_t* out);  /* N or N+1 */
po_status po_result_weights(const po_result* r, double* out, size_t n);
po_status po_result_iterations(const po_result* r, int32_t* out);  /* total, per method; diagnostics only */
po_status po_result_multiplier(const po_result* r, int32_t constraint, double* out);
          /* PO_CONSTRAINT_BUDGET, _TARGET_RETURN, _VOLATILITY, _LEVERAGE */
po_status po_result_bound_multipliers(const po_result* r, double* out, size_t n);
po_status po_result_residual(const po_result* r, double* out);
          /* risk parity: sum (c_i - b_i)^2 at the point; 0 elsewhere; NaN when the weights are NaN */
void po_result_free(po_result* r);
```

- Erros de chamada (`INVALID_ARGUMENT`) são entradas que não fazem sentido: NaN ou infinito em μ, Σ, R ou b; Σ
  assimétrica ou não definida positiva; tamanhos incompatíveis; setter fora do domínio (δ ≤ 0, λ ∉ [0,1), lo > hi);
  combinação objetivo × setter fora da matriz da seção 3.3. Um setter com valor inválido deixa a spec inalterada.
- Status do resultado (`OPTIMAL`, `INFEASIBLE`, `MAX_ITER`) é sobre o problema bem formado: restrições
  incompatíveis, vol-alvo abaixo da mínima, Sharpe sem excesso positivo, risk parity sem convergir. Um resultado
  com status diferente de `OPTIMAL` existe e pode ser lido. Pesos e multiplicadores são NaN quando o status não
  fornece ponto: sempre em `INFEASIBLE` e `UNBOUNDED`, e em `MAX_ITER` dos objetivos resolvidos por QP (todos da v1
  exceto `RISK_PARITY`), porque o iterado dual do solver não é primal viável; a semântica de `MAX_ITER` do
  `RISK_PARITY` é documentada quando ele for acrescentado.
- A spec e os dados são copiados no `po_solve`; mudá-los depois não muda nenhum resultado. `po_result` é imutável e
  pode ser lido de várias threads. `po_solve` pode ser chamado concorrentemente com handles diferentes ou com os
  mesmos handles somente leitura.
- Multiplicadores: são os de KKT do problema exatamente como escrito na seção 3.3, nas variáveis do usuário (w,
  nunca nas homogeneizadas nem no x do Kelly), com a escala do objetivo escrito (min w'Σw, sem ½). No MAX_SHARPE e
  no MAX_DIVERSIFICATION são os do problema de razão em w, obtidos dos homogeneizados por mudança de escala. No
  RISK_PARITY, os do problema de mínimos quadrados no ponto de convergência. Uma restrição que não existe no
  problema (orçamento com caixa, target_return não definido) devolve NaN, sem erro. A convenção de sinal é do
  builder: uma só para todos os objetivos, documentada no header. Cada multiplicador é a derivada do f ótimo em
  relação ao lado direito da sua restrição (uma restrição `<=` tem multiplicador `<= 0`, como beta num bound
  superior); com restrições dependentes os multiplicadores não são únicos e um conjunto válido é devolvido.

### 5.3 Métricas

Funções sobre dados e pesos (não sobre o resultado), para servirem a qualquer carteira, inclusive uma de pesos
iguais. Os pesos têm o tamanho N ou N+1 dos dados.

```c
po_status po_portfolio_return(const po_data* d, const double* w, size_t n, double* out);
po_status po_portfolio_volatility(const po_data* d, const double* w, size_t n, double* out);
po_status po_risk_contributions(const po_data* d, const double* w, size_t n, double* out); /* soma 1 */
po_status po_sharpe_ratio(const po_data* d, const double* w, size_t n, double* out);
po_status po_diversification_ratio(const po_data* d, const double* w, size_t n, double* out);
po_status po_utility(const po_data* d, const double* w, size_t n, double delta, double* out);
po_status po_tracking_error(const po_data* d, const double* w, size_t n, double* out);
          /* root mean square, not demeaned: sqrt(sum((Rw - b)^2) / T), the quantity TRACKING_ERROR
             minimizes (author's decision, 2026-11-15) */
```

Leverage (Σ|w|) e contribuição de retorno (wᵢμᵢ / retorno) ficam fora: são uma linha em qualquer linguagem.

### 5.4 Lote (depois, sem quebra)

Quando vier: `po_context` (pool de threads, como o `ve_context`) e
`po_solve_batch(ctx, specs[], datas[], count, results[])`, sobre os handles que já existem. Internamente o lote é
um laço de `po_solve` no pool, sem caminho duplicado. Como os dados são copiados no handle, nada no contrato da v1
precisa mudar para isso.

## 6. Python

ctypes + numpy, pacote `portopt`, gerenciado com uv, como o VolEngine (`_native.py` com descoberta da biblioteca
por `PORTOPT_LIB`, checagem de ABI, protótipos e posse dos handles).

```python
import portopt as po

data = po.Data(mu, cov, risk_free=0.10)                  # cópia; returns=, benchmark= para tracking error
spec = po.Spec(po.Objective.MAX_SHARPE, bounds=(0.0, 0.3), leverage=1.5)
res  = po.solve(spec, data)

res.status, res.weights, res.iterations                  # weights: numpy (N ou N+1)
data.volatility(res.weights), data.risk_contributions(res.weights), data.sharpe_ratio(res.weights)

frontier = [po.solve(po.Spec(po.Objective.MIN_VARIANCE, target_return=r), data) for r in targets]
```

- `bounds=` aceita escalares ou arrays; o wrapper escolhe `set_bounds` ou `set_asset_bounds`.
- Status diferente de `OPTIMAL` não levanta exceção; erro de chamada levanta (`ValueError` para
  `INVALID_ARGUMENT`).
- Métricas como métodos de `Data` recebendo pesos.

## 7. C++

`cpp/portopt.hpp`, header-only C++23 sobre a interface C, com `std::expected`, `std::span` e RAII nos handles,
como o `volengine.hpp`.

## 8. Performance

Não é a motivação principal, mas vale o método do VolEngine (benchmark, profile, assembly antes de aceitar
otimização que piora a leitura). Goldfarb-Idnani denso é O(N³) por fatoração e adequado até alguns milhares de
ativos. Não há meta numérica na v1; entra uma suíte Google Benchmark com casos de N = 10, 100 e 1000 para ter
linha de base.

## 9. Validação

Não há oráculo R: as equações são derivações do autor. A validação verifica a matemática diretamente, e as
tolerâncias por grandeza ficam em `docs/guide/validation.md`.

- **Condições de KKT** em todo teste de QP: estacionariedade, viabilidade primal, viabilidade dual (sinal dos
  multiplicadores de desigualdade) e complementaridade, nas variáveis do usuário.
- **Formas fechadas** onde existem: mínima variância Σ⁻¹1 / 1'Σ⁻¹1; tangente Σ⁻¹(μ − rf) normalizada (só quando
  1'Σ⁻¹(μ − rf) > 0; caso contrário o teste espera `UNBOUNDED`); utilidade e
  Kelly sem restrições; fronteira sem bounds (parábola em (retorno, variância)).
- **Consistência entre caminhos**: bounds folgados dão o mesmo resultado da forma fechada; Sharpe homogeneizado
  bate com o máximo da fronteira; vol-alvo bate com o ponto da fronteira de mesma volatilidade; Spinu e
  Lagrangiano aumentado concordam quando ambos se aplicam.
- **Invariantes**: contribuições de risco somam 1 e batem com b no risk parity; nenhuma carteira viável aleatória
  (semente fixa) tem Sharpe ou razão de diversificação maior que o ótimo; invariância de escala de Σ onde a teoria
  diz.
- **Casos degenerados**: Σ quase singular, bounds que fixam ativos, restrições incompatíveis (`INFEASIBLE`),
  leverage ativo e inativo, problema de 1 ativo.
- Matrizes aleatórias definidas positivas com semente fixa, geradas no teste; sem fixtures externas.
- O header compila como C11 (`c_header_check.c`, como no VolEngine).

## 10. Dependências

| Área | Escolha |
|---|---|
| QP, Newton, Lagrangiano aumentado | Próprios |
| Álgebra linear | Própria; sem Eigen |
| Testes / benchmarks | GoogleTest, Google Benchmark (FetchContent) |
| Build | CMake + `CMakePresets.json` (presets como os do VolEngine) |
| Python | ctypes + numpy; uv |

## 11. Plataformas, distribuição, licença

Iguais ao VolEngine: C++23; CI em MSVC, Clang e GCC 15+; wheels para Windows x64, macOS arm64 e x64, Linux glibc e
musl em x64 e arm64; pacote no PyPI, biblioteca nativa nas GitHub Releases; licença MIT; `CITATION.cff`. O build
MinGW liga `libstdc++`, `libgcc` e `winpthread` estaticamente, como no VolEngine.

## 12. Estrutura de pastas

```
PortfolioOptimization/
  include/portopt/         portopt.h (o contrato)
  src/                     core/, linalg/, qp/, problems/, riskparity/, metrics/, capi/
  cpp/                     portopt.hpp + testes
  python/                  pacote portopt (src/, tests/, pyproject.toml)
  tests/  unit/ style/ benchmarks/
  docs/   guide/ (math.md com as equações, validation.md)  performance/  superpowers/
  .github/workflows/
  ARCHITECTURE.md STYLEGUIDE.md CONTRIBUTING.md CHANGELOG.md README.md LICENSE CITATION.cff
  .clang-format .clang-tidy .editorconfig CMakeLists.txt CMakePresets.json
```

STYLEGUIDE.md, `.clang-format`, `.clang-tidy`, `.editorconfig` e `check_style.py` partem dos do VolEngine,
ajustados aos nomes deste projeto (prefixo, pastas). Documentos do repositório em inglês; documentos de trabalho em
`docs/superpowers/` em português. `docs/guide/math.md` registra cada equação da seção 3 com sua derivação, já que
são a referência da validação.

## 13. Deixado ao builder

- Nomes exatos de funções, enums e handles, respeitando a seção 5.1.
- Organização interna do núcleo, nomes de pastas de `src/`.
- Critérios de parada, tolerâncias numéricas, tratamento de degenerescência no Goldfarb-Idnani.
- Busca em δ do `MAX_RETURN` (bissecção ou outro método de raiz).
- Tratamento exato do leverage, se encontrar um melhor que a família implícita.
- Convenção de sinal dos multiplicadores (uma só, documentada).
- Ordem das fases. Sugestão: fundação (repo, build, estilo, CI) e `MIN_VARIANCE` de ponta a ponta (núcleo →
  interface C → Python → validação por KKT); depois os outros objetivos de QP; depois homogeneização (Sharpe,
  diversificação); leverage; vol-alvo; risk parity; wrapper C++ e empacotamento.

## 14. Decisões tomadas pelo Claude (para conferir)

- Nome `portopt`, prefixo `po_`, pacote Python `portopt`.
- Sem `pesos_carteira`, `merton_portfolio`, `portfolios`; `momentum_utility_portfolio` sai junto com o GA.
- Leverage é Σ|w| em todos os objetivos (o `port_vol_target` do R usava norma-2) e conta só os ativos de risco.
- Bounds escalares valem também para o caixa.
- `MAX_RETURN` com vol abaixo da mínima variância devolve `INFEASIBLE` (o R devolvia a mínima variância).
- Caixa sai com peso 0 em `MAX_SHARPE` e `MAX_DIVERSIFICATION` (bounds do caixa que excluem 0 dão `INFEASIBLE`);
  `RISK_PARITY` rejeita caixa.
- Status `UNBOUNDED` para Sharpe e diversificação com supremo não atingido.
- Tracking error com caixa pela eliminação de w₀ (sem coluna constante); `cov` opcional nos dados.
- KELLY e IMPORTANCE implementados como reduções exatas de UTILITY; leverage e bounds do Kelly sobre w.
- Multiplicadores de KKT do problema como escrito, em w; NaN para restrição ausente.
- MAX_RETURN com vol inativa: carteira de retorno máximo, `OPTIMAL`, multiplicador zero.
- Kelly com caixa equivale ao `use_rf = TRUE` do R; não há Kelly com rf como limiar sem caixa investível (se
  `risk_free` está nos dados, o caixa está no universo).
- Risk parity sem bounds devolve a solução de pesos positivos.
- Padrões δ = 5 e λ = 0.5 (os do R).
- Nenhuma anualização; tudo nas unidades do usuário.
- Dados copiados no handle (não emprestados), para não haver bug de tempo de vida entre linguagens e para o lote
  futuro não exigir mudança de contrato.
- Métricas como funções sobre dados e pesos; leverage e contribuição de retorno ficam fora.
- Status de resultado em vez de NULL silencioso; erros de chamada só para entradas sem sentido.
- Sem meta numérica de performance na v1.
