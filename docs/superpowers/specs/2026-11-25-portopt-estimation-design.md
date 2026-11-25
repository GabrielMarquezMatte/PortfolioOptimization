# PortOpt 1.1: utilidades de estimação

Data: 2026-11-25. Origem: as propostas do DataPanel (`DataPanel/docs/bibliotecas/portopt-propostas.md`).

## 1. Intenção

O DataPanel calcula o retorno esperado em dois lugares: no batch em Python e na API em C#. Os dois precisam chegar ao
mesmo número, e esse número é a entrada da otimização (`mu` de `po_data_create`, `MarketInputs.Mu`). Hoje cada lado tem
a sua implementação. Com os estimadores no PortOpt, as duas linguagens chamam o mesmo código C e recebem os mesmos bits.

Os estimadores entram como **utilidades sem efeito colateral**. Leem as entradas (arrays, um `po_data`, um `po_views`)
e escrevem só as saídas. O usuário passa o resultado ao `po_data_create` por conta própria. O core da otimização
(`solve`, `problems`, `qp`, `riskparity`, `metrics`) não chama esse módulo nem muda por causa dele. Nos bindings
Python e C#, as utilidades são funções puras que recebem arrays, e os handles ficam escondidos.

Isso revê o não-objetivo "Estimar μ ou Σ [...] Os dados chegam prontos" de
[2026-11-15-portopt-design.md](2026-11-15-portopt-design.md). Para o solve, os dados continuam chegando prontos. A
estimação passa a existir ao lado do solve, como utilidade. Aquele documento ganha uma nota que aponta para este.

## 2. Objetivos

Funções novas na API C, nos bindings Python, C# e C++, com testes e documentação:

| Função | O que calcula |
|---|---|
| `po_sample_moments` | média e covariância amostral (divisor `t - 1`) de um painel de retornos `t x n` |
| `po_bayes_stein` | a média com o encolhimento de Bayes-Stein de Jorion (1986), a partir de `mu_hat`, `S` e `t` |
| `po_bl_equilibrium` | o retorno de equilíbrio do Black-Litterman, `Pi = delta Sigma w_mkt` |
| `po_views_create`, `po_views_set_omega`, `po_views_free` | o handle das visões: P, Q e Ω |
| `po_bl_omega_he_litterman` | `Omega = diag(tau P Sigma P')` |
| `po_bl_omega_idzorek` | `Omega` de Idzorek na forma fechada de Walters, `omega_k = ((1 - c_k) / c_k) tau p_k Sigma p_k'` |
| `po_bl_posterior` | o μ posterior e, se pedida, a covariância posterior |

Além delas, quando Σ não é positiva definida, a mensagem passa a dar o índice do pivô que falhou (§6.7), e os bindings
ganham a conversão do posterior para a entrada da otimização (§8.1).

O alvo é a 1.1.0. As funções são adicionadas, e `PO_ABI_VERSION` continua em 1 ([versioning.md](../../guide/versioning.md)).
Este trabalho só entra no código: as entradas ficam em `Unreleased` no CHANGELOG e em `PublicAPI.Unshipped.txt`, e a
versão no CMake continua 1.0.0. A release é uma tarefa separada, feita depois que entrar tudo o que for possível das
propostas do DataPanel (decisão do autor).

## 3. Não-objetivos

- Calcular retornos a partir de preços, anualizar, alinhar datas, tratar lacunas. O painel chega alinhado e finito,
  e a anualização é do usuário.
- Valores padrão de `delta`, `tau` e `rf`, e de onde vêm os pesos de mercado. Ficam no DataPanel.
- A matriz positiva definida mais próxima (Higham, 2002) e o menor autovalor na mensagem de erro. Os dois exigem uma
  decomposição em autovalores (um Jacobi simétrico portável), que fica para um design separado, feito antes da release
  1.1.0 (decisão do autor).
- A covariância preditiva de Jorion, Ledoit-Wolf e outros estimadores. Ficam para quando alguém pedir.
- Códigos de status novos.

## 4. Anti-objetivos

- **O core depender de `estimation`.** A dependência vai só de `capi` para `estimation`, e de `estimation` para
  `linalg` e `core`.
- **Mais de 5 parâmetros numa função**, a C incluída ([STYLEGUIDE.md](../../../STYLEGUIDE.md#parameters--maximum-5-per-function)).
- **Configuração opcional em parâmetro.** Cada função faz uma conta fixa. Uma variante nova vira uma função nova, como
  manda a regra de estabilidade.
- **Correção silenciosa.** Pesos de mercado que não somam 1 não são renormalizados, e `w_hat` não é truncado. A entrada
  fora do domínio vira erro.
- **Bits diferentes entre plataformas.** As contas usam só `+`, `-`, `*`, `/` e `sqrt`, em ordem de soma fixa, sob o
  `-ffp-contract=off` que `portopt_core` e `portopt_capi` já usam.

## 5. Abordagens consideradas

- **A. Funções puras com arrays soltos, cada uma com uma conta fixa.** O autor escolheu esta primeiro. Na hora do
  plano, apareceu o conflito com a regra de no máximo 5 parâmetros, que vale também para a interface C e que o
  `clang-tidy` aplica: o posterior teria 10 parâmetros, e outras três funções passariam de 5. Descartada.
- **A'. Σ via `po_data` e as visões via um handle `po_views` (escolhida pelo autor).** O `po_data` já valida Σ
  (simetria, Cholesky) e guarda o fator, e o `po_views` guarda P, Q e Ω. Toda função fica com até 5 parâmetros. As
  funções não mudam os handles, e os bindings mostram funções puras sobre arrays. Os custos: o `po_views` tem estado (Ω
  por setter), e o `mu` do `po_data` significa μ̂ no Bayes-Stein e Π no posterior.
- **A com exceção à regra** (`NOLINT` nas quatro funções). Abriria um precedente justamente na regra que veio dos
  parâmetros posicionais soltos. Descartada.
- **B. Um handle com P, Q, τ e a forma de Ω, mais `po_black_litterman(views, ...)`.** Põe configuração (a forma de Ω)
  no estado. Descartada em favor de A', em que cada forma de Ω é uma função que calcula a matriz.
- **C. Implementar em cada binding.** Perde a garantia de mesmo número entre Python e C#, que é o motivo do trabalho.
  Descartada.

Na entrada do Bayes-Stein, a decisão do autor foi separar os momentos: `po_sample_moments` calcula `mu_hat` e `S` do
painel, e `po_bayes_stein` recebe os momentos (num `po_data`). Os momentos também servem direto ao `po_data_create`,
sem encolhimento.

## 6. Contrato C

Convenções comuns (as mesmas de `portopt.h`):

- As matrizes são row-major. `n` é o número de ativos, `t` o de períodos e `k` o de visões.
- As saídas são arrays do chamador. O tamanho vem dos handles: `n` é `po_data_asset_count`, e `k` é o `k` dado a
  `po_views_create`. Um handle nulo dá `PO_INVALID_HANDLE`, e um ponteiro obrigatório nulo dá `PO_INVALID_ARGUMENT`.
- Toda entrada fora do domínio dá `PO_INVALID_ARGUMENT`, com o motivo em `po_last_error()`. Nenhum código novo.
- Em uma falha, nenhuma saída é escrita. As contas usam temporários, e a cópia para a saída acontece no fim.
- **O `po_data` das utilidades:** tem de ter covariância. Um `po_data` com rf (`po_data_set_risk_free`) é recusado,
  porque o rf não tem papel nessas contas e ignorá-lo seria falhar em silêncio. As séries de retornos
  (`po_data_set_returns`) não são usadas. Σ é a covariância já validada e simetrizada pelo `po_data_create`, com o
  fator de Cholesky guardado.

### 6.1 `po_sample_moments`

```c
PO_API po_status po_sample_moments(size_t t, size_t n, const double* returns, double* out_mu, double* out_cov);
```

- `returns` é `t x n`, no mesmo layout de `po_data_set_returns`. Exige `t >= 2`, `n >= 1`, valores finitos e `t * n`
  sem overflow.
- `out_mu[j]` é a média da coluna `j`. `out_cov` (`n x n`) é a covariância calculada em duas passadas (a média, depois
  os desvios), com divisor `t - 1`. O triângulo superior é calculado e espelhado, então a saída é exatamente simétrica.
- Não confere se a saída é positiva definida. Com `t <= n`, S é singular, o que é uma saída legítima: o
  `po_data_create` recusa S depois, se ela for usada como covariância.

### 6.2 `po_bayes_stein`

```c
PO_API po_status po_bayes_stein(const po_data* d, size_t t, double* out_mu, double* out_grand_mean,
                                double* out_weight);
```

- O `mu` de `d` é `mu_hat`, e a covariância de `d` é a S amostral com divisor `t - 1` (as saídas de
  `po_sample_moments`). Exige `mu` em `d` e `t > n + 2`.
- Cálculo, com `1` o vetor de uns:
  - `Sigma = S (t - 1) / (t - n - 2)`;
  - `mu_0 = (1' Sigma^-1 mu_hat) / (1' Sigma^-1 1)`;
  - `w_hat = (n + 2) / ((n + 2) + t (mu_hat - mu_0 1)' Sigma^-1 (mu_hat - mu_0 1))`;
  - `mu_BS = (1 - w_hat) mu_hat + w_hat mu_0 1`.
- `Sigma^-1` nunca é formada. As resoluções usam o fator de S que o `po_data` guarda, com
  `Sigma^-1 x = S^-1 x (t - n - 2) / (t - 1)`.
- Saídas: `out_mu` (`n`), `out_grand_mean` (`mu_0`) e `out_weight` (`w_hat`), todas obrigatórias.
- **A documentação diz em destaque:** `w_hat` depende da frequência, porque multiplicar `mu_hat` e `S` por 252
  multiplica a forma quadrática por 252. A conta é feita na frequência original dos dados (diária), e só o μ resultante
  é anualizado depois. Mudar a unidade (multiplicar os retornos por uma constante) não muda `w_hat`.

### 6.3 `po_bl_equilibrium`

```c
PO_API po_status po_bl_equilibrium(const po_data* d, const double* w_mkt, size_t n, double delta, double* out_pi);
```

- `Pi = delta Sigma w_mkt`. Exige `n` igual à contagem de ativos de `d`, `w_mkt` finito com `|1'w_mkt - 1| <= 1e-8` e
  `delta` finito e `> 0`. O `mu` de `d`, se houver, não é usado.
- `Pi` é retorno **em excesso do rf**. A documentação diz que, com `po_data_set_risk_free(rate)`, o core usa
  `mu - rate`, então o usuário passa `Pi + rate` (ou o μ posterior `+ rate`) ao `po_data_create` da otimização. A
  conversão de §8.1 faz isso.

### 6.4 `po_views`

```c
typedef struct po_views po_views;
PO_API po_status po_views_create(size_t k, size_t n, const double* p, const double* q, po_views** out);
PO_API po_status po_views_set_omega(po_views* v, const double* omega, size_t k);
PO_API void po_views_free(po_views* v);
```

- `p` é `k x n`, e `q` tem `k` entradas. `Q` está na mesma unidade de `Pi`, em excesso do rf (decisão do Claude, a
  convenção do BL). Visões relativas (+1/−1) não dependem disso, porque o rf se cancela. Visões absolutas dependem, e
  a documentação diz isso.
- `po_views_create` exige `n >= 1`, `k <= n`, `p` e `q` finitos e P com posto de linha cheio: Householder QR de `P'`
  (`n x k`) e `has_full_column_rank` de `linalg/qr.hpp`. Na falha do posto, a mensagem fala em visões redundantes.
- `k = 0` é aceito, com `p` e `q` podendo ser NULL: são as visões vazias, e o posterior devolve `mu = Pi` e
  `Sigma_post = (1 + tau) Sigma`.
- `po_views_set_omega` copia Ω (`k x k`). Exige o mesmo `k` do create, valores finitos, simetria pela tolerância de Σ
  (`1e-10` relativa à maior das duas entradas da diagonal) e diagonal `>= 0`. Uma falha mantém o Ω anterior. Ω é
  esperada positiva semidefinida. Fora a diagonal `>= 0` e o Cholesky de `A` no posterior, isso não é conferido:
  conferir exige autovalores, e os dois helpers sempre produzem Ω válida.
- As convenções de handle valem como nos outros: depois de uma falha no create, `*out` é NULL, e o `free` aceita NULL.

### 6.5 `po_bl_omega_he_litterman` e `po_bl_omega_idzorek`

```c
PO_API po_status po_bl_omega_he_litterman(const po_data* d, const po_views* v, double tau, double* out_omega);
PO_API po_status po_bl_omega_idzorek(const po_data* d, const po_views* v, double tau, const double* confidence,
                                     double* out_omega);
```

- `out_omega` é `k x k`, diagonal, com zeros fora dela. O usuário passa essa matriz a `po_views_set_omega`.
- He e Litterman: `omega_kk = tau p_k Sigma p_k'`. Idzorek (Walters): `omega_kk = ((1 - c_k) / c_k) tau p_k Sigma p_k'`,
  com `confidence` de `k` entradas e `0 < c_k <= 1`. `c_k = 1` dá `omega_kk = 0`, a visão de confiança total.
- Exigem o mesmo `n` em `d` e `v` e `tau` finito e `> 0`. O Ω já guardado em `v` não é lido.

### 6.6 `po_bl_posterior`

```c
PO_API po_status po_bl_posterior(const po_data* d, const po_views* v, double tau, double* out_mu, double* out_cov);
```

- O `mu` de `d` é `Pi`, e a covariância de `d` é Σ. Exige `mu` em `d`, o mesmo `n` em `d` e `v`, `tau` finito e `> 0`
  e, com `k >= 1`, Ω dada por `po_views_set_omega`.
- Com `A = P (tau Sigma) P' + Omega` e `B = tau Sigma P'` (`n x k`):
  - `mu = Pi + B A^-1 (Q - P Pi)`;
  - `Sigma_post = Sigma + tau Sigma - B A^-1 B'`, a forma sem `Omega^-1`, que aceita `Omega = 0`. A saída é
    exatamente simétrica: o termo `B A^-1 B'` é calculado como `X'X`, com `X = L^-1 B'`, só no triângulo superior, e
    espelhado.
- A resolução usa Cholesky de `A`, nunca `A^-1` explícita. Se `A` não for positiva definida, o erro informa o pivô
  (§6.7).
- `out_cov` pode ser NULL, e nesse caso a covariância posterior não é calculada.

### 6.7 O pivô na mensagem de matriz não positiva definida

Decisão do autor: entra neste trabalho. `linalg::cholesky` passa a informar o índice do primeiro pivô recusado, e não
só `nullopt`. A forma da mudança fica com quem constrói. A mensagem de Σ não positiva definida do `po_data_create`
passa a dar esse índice, com base zero como as outras mensagens. Por exemplo:
`covariance is not positive definite: pivot 3 is not positive`. A mensagem de `A` no posterior segue o mesmo padrão.
Mensagens podem mudar em qualquer versão ([versioning.md](../../guide/versioning.md#error-messages)), então isso não
quebra contrato. Os chamadores atuais de `cholesky` (o QP, o Spinu, a descida poliédrica e `Data::create`) continuam
com o mesmo comportamento, e os testes existentes garantem isso.

## 7. Onde fica no código

- `src/estimation/`: `moments`, `bayes_stein`, `black_litterman` e `views` (`.hpp` e `.cpp`). As funções devolvem
  `std::expected<..., Failure>`. O módulo depende de `linalg` e de `core` (`Data`, `Failure`). Entra em
  `portopt_core`.
- `src/capi/estimation/`: `views_api.cpp` (create, set_omega, free) e `estimation_api.cpp` (as seis utilidades), com
  `guard`, `null_handle`, `null_argument` e `set_error`, como em `metrics_api.cpp`. Fica numa subpasta porque
  `src/capi` já tem 8 módulos, o limite por pasta. `struct po_views` vai para `capi/handles.hpp`.
- `include/portopt/portopt.h` e `portopt.abi`: as declarações, numa seção "Estimation utilities" depois das métricas.
- `ARCHITECTURE.md`: o módulo no diagrama de camadas e nas tabelas.

## 8. Bindings (decisões do Claude)

Python e C# expõem funções puras sobre arrays e criam por dentro o `po_data` e o `po_views`. Todos os bindings passam
`out_cov` sempre e devolvem as duas saídas do posterior.

- **Python**, em `_estimation.py`, exportado no topo de `portopt`:
  - `sample_moments(returns) -> tuple[mu, cov]`;
  - `bayes_stein(mu_hat, cov, t) -> BayesStein` (named tuple `mu`, `grand_mean`, `weight`);
  - `bl_equilibrium(cov, w_mkt, delta) -> ndarray`;
  - `bl_omega_he_litterman(cov, p, tau)` e `bl_omega_idzorek(cov, p, tau, confidence)`;
  - `bl_posterior(cov, pi, tau, p, q, omega) -> BlackLittermanPosterior` (§8.1).

  Formas erradas dão `ValueError` antes da chamada nativa, como em `_data.py`. `p` vazio `(0, n)` representa `k = 0`.
- **C#:** a classe estática `Estimation`, com `SampleMoments`, `BayesStein`, `BlackLittermanEquilibrium`,
  `OmegaHeLitterman`, `OmegaIdzorek` e `BlackLittermanPosterior`. Os retornos são os records `Moments(Mu, Covariance)`,
  `BayesSteinEstimate(Mu, GrandMean, Weight)` e `BlackLittermanPosterior(Mu, Covariance)`. As matrizes são `double[,]`,
  como em `MarketInputs`. As entradas vão para `PublicAPI.Unshipped.txt`, as `LibraryImport` para `NativeMethods.cs` e
  o `ViewsHandle` para `Handles.cs`. Os testes de contrato e de cobertura existentes passam a cobrir as funções novas.
- **C++ (`cpp/include/portopt/portopt.hpp`):** a classe `Views` (RAII, como `Data`) e funções livres com os nomes sem
  `po_`, recebendo `const Data&` e `const Views&` e devolvendo `Expected<...>`. O `check_abi_coverage.py` exige que o
  wrapper use toda função `PO_API`.

### 8.1 Do posterior à otimização (proposta 4, decisão do autor)

O posterior do BL ganha, nos bindings, uma conversão para a entrada da otimização. Ela resolve uma armadilha: o μ do
BL está em excesso do rf, e quem usa cash precisa somar o rate antes do `po_data_create`.

- **C#:** `BlackLittermanPosterior.ToMarketInputs(double? riskFree = null)` devolve
  `MarketInputs { Mu = mu + rf, Covariance = covariância posterior, RiskFree = riskFree }`. Sem `riskFree`, o μ fica
  como está e não há cash.
- **Python:** `bl_posterior` devolve o `NamedTuple` `BlackLittermanPosterior(mu, cov)`, com o método
  `to_data(risk_free=None) -> Data`, que segue a mesma regra.
- **C e C++:** nada. Na API C, a conversão é somar um escalar a um vetor (decisão do Claude: o DataPanel usa Python e
  C#, e o wrapper C++ ganha a conversão quando alguém pedir).
- Usa a covariância posterior, como a proposta pede. Quem quiser Σ em vez dela monta a entrada à mão.

## 9. Testes e referências

- **C++, em `tests/unit/estimation/`:**
  - momentos de um painel pequeno com valores à mão; `t = 1` recusado; saída exatamente simétrica;
  - Bayes-Stein com `mu_hat = c 1`, que dá `mu_0 = c`, `w_hat = 1` e `mu_BS = mu_hat`; a mudança de unidade não muda
    `w_hat`; multiplicar `mu` e `S` por 252 muda `w_hat`; `t <= n + 2` recusado; `po_data` com rf recusado;
  - BL com `k = 0`, que dá `mu = Pi` e `Sigma_post = (1 + tau) Sigma`; Ω = 0 (e Idzorek com `c = 1`), em que `P mu`
    reproduz `Q` dentro do arredondamento; com Ω > 0, `Sigma_post` igual à forma
    `Sigma + ((tau Sigma)^-1 + P' Omega^-1 P)^-1`; P redundante recusado no create; `w_mkt` com soma fora da tolerância
    recusado; posterior com `k >= 1` sem Ω recusado.
- **C API, em `tests/unit/capi/estimation/`:** handle nulo, ponteiro nulo, saída intocada na falha, mensagem em
  `po_last_error`, `*out` NULL depois de um create que falha.
- **Pivô:** uma Σ com o pivô 2 não positivo dá a mensagem com o índice 2 no `po_data_create`. `cholesky_test.cpp`
  confere o índice devolvido.
- **Python:** referências em numpy (`np.mean`, `np.cov(ddof=1)`, Bayes-Stein e BL escritos com `np.linalg.solve`),
  incluindo a visão relativa +1/−1 pedida pelo DataPanel.
- **C#:** um teste por função contra os mesmos valores do Python, mais os testes de contrato e de cobertura.
- **Conversão:** `ToMarketInputs` e `to_data` com e sem rf. Com rf, o `Mu` é o posterior `+ rf` e há cash. Sem rf, o
  `Mu` é o posterior e não há cash. Nos dois casos, a covariância é a posterior.
- **Determinismo:** os valores de um caso fixo ficam registrados em hex-float e conferidos em toda a matriz de CI,
  como já se faz com os solves.

## 10. Documentação

- `docs/guide/math.md`, seção nova "Estimation utilities": as fórmulas, a nota de frequência do Bayes-Stein, a nota de
  excesso do rf do BL e as referências (Black e Litterman, 1992; He e Litterman, 1999; Idzorek, 2005; Walters, 2014;
  Jorion, 1986).
- `docs/guide/validation.md`: as regras de domínio e as tolerâncias de §6.
- `CHANGELOG.md`, em Unreleased, Added; os READMEs de Python, C# e C++ com um exemplo de BL até o `solve`, usando uma
  visão relativa.

## 11. Fica para quem constrói

- Os nomes internos em C++ e a divisão exata em arquivos de `src/estimation/`.
- A forma da mudança em `linalg::cholesky` (§6.7).
- As mensagens de erro, fora o padrão do pivô.
- Escolher o caso fixo do teste de determinismo.
