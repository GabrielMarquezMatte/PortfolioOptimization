# PortOpt — Plano 5: multiplicadores precisos no RISK_PARITY com `target_volatility`

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fazer os multiplicadores do `RISK_PARITY` com `target_volatility` terem a precisão que o ponto KKT permite (hoje λ_vol erra 3,5e-8 a cond(Σ) = 1e6 e até 7,4e-4 a 1e8, com os pesos a 1e-16..1e-12), sem mexer nos pesos nem nas faixas.

**Architecture:** O λ de mínimos quadrados (LS) em qualquer w double erra por ≈ ε·|H||w|, porque H = 2J'J é enorme. O par KKT (w, λ) não é mal condicionado: um passo de Newton-KKT no ponto final, em double, com λ + dλ mantido, dá λ a 1e-15..1e-13 (cond 1e4), ~1e-12 (1e6), 1e-12..2e-10 (1e8) e 1e-10..1e-8 (1e10). Não é preciso w estendido. A Fase 1 mantém o dλ do último sistema de Newton em `newton_polish.cpp`, corrige os multiplicadores de bound dos pesos presos por H·dw e protege o caso próximo de σ_min. O portão KKT continua com o λ de LS.

**Tech Stack:** o dos planos 1–4. Sem `long double` em `src/`; sem dependência nova.

**Spec / pesquisa:** `docs/research/2026-11-16-multiplier-accuracy.md` (a fonte das decisões; seções 0 e 6 são o roteiro). Spec: `docs/superpowers/specs/2026-11-15-portopt-design.md` §5.2 (multiplicadores). Plano anterior: `docs/superpowers/plans/2026-11-16-portopt-04-risk-parity.md`.

## Global Constraints

- Valem as do plano 4 (interface C só por acréscimo, convenção de multiplicadores do header, estilo do `STYLEGUIDE.md`, commits **sem** `Co-Authored-By`, gcc-release apenas, documentos de `docs/superpowers/` em português).
- **Nenhuma faixa afrouxada. Tolerância é medida, nunca chutada:** todo teto novo é medido e dito "medido, não derivado".
- Os pesos, o status e o portão `is_kkt_point` não mudam de comportamento; só mudam os multiplicadores devolvidos. A suíte inteira (317 testes) e as baterias de fuzz continuam verdes sem mudar um valor de teste que não seja um teto de multiplicador.
- **Metas de precisão** (se a Task 0 confirmar): λ ≤ 1e-12 relativo até cond(Σ) = 1e6; ≤ 1e-9 a 1e8 (folga 10 sobre o medido); a 1e10 fica "não validado" (a Fase 1 dá ~1e-4).
- **Fora de escopo:** aritmética double-double com w_ext (Fase 2 da pesquisa, só se a meta a 1e8/1e10 exigir; vira plano próprio), o caminho QP (uma medida só, na Task 2), `long double`.

## Review Focus

1. **Salvaguarda perto de σ_min.** s = σ_min(1 + δ), δ ∈ {1e-2, …, 1e-10} (tabela de `validation.md`, seção do `RISK_PARITY`): o dλ pode explodir onde a LICQ quase falha. Teste: `RiskParity.NearMinimumVolatilityKeepsLeastSquaresMultipliers` (Task 1).
2. **Portão × multiplicadores novos.** Com λ\* o resíduo de estacionaridade em w é H·δw inteiro e pode passar de 2× o piso onde o erro progressivo de w passa de ~ε. O portão e `expect_risk_parity_kkt` ficam com o λ de LS; o λ novo é conferido **contra a referência**, não pela estacionaridade. Teste: `expect_risk_parity_kkt` com o λ do resultado e com o λ de LS refeito no teste (Task 1).
3. **Multiplicadores de bound de pesos presos** (`LongOnly`/`Capped` com vol-alvo): β_j sofre o mesmo H·δw; sinal e complementaridade continuam valendo. Teste: `RiskParityFuzz.TargetVolatilityBoundMultipliersByCondition` (Task 2).
4. **Identidade de homogeneidade** λ_orçamento = −λ_retorno·r − λ_vol·s − Σβᵢwᵢ continua valendo com os multiplicadores novos. Teste existente `RiskParity.HomogeneityIdentity` com cond 1e6 (Task 2).
5. **Caminho só com `target_return`** (já a 1e-14 a 1e8) não pode piorar: regressão em `RiskParityFuzz.*` existente (Task 1).

---

### Task 0: Medir o dλ com a Hessiana real (confirma ou refuta antes de implementar)

**Files:**
- Create: `tests/unit/fuzz/risk_parity/risk_parity_multiplier_measure_test.cpp` (teste `DISABLED_MeasureNewtonMultipliers`, padrão dos outros `DISABLED_Measure…`)
- Modify: `CMakeLists.txt`, `docs/guide/validation.md` (resultado)

**Interfaces:**
- Consumes: `riskparity::newton_system` / `evaluate_kkt` (`kkt.hpp`), `forward_error` e a refinação `long double` de `risk_parity_refinement.hpp`, o script de 80 dígitos da pesquisa (`probe6.py` / `exp_mult.py`, que podem ter sido apagados: o `.md` da pesquisa descreve o método, §0).
- Produces: uma tabela medida por cond ∈ {1e4, 1e6, 1e8, 1e10} e casos (n, seed) ∈ {(5,7), (5,11), (8,3)}, sem box, s = 1,2·σ_min, orçamento igual.

- [ ] **Step 1:** Teste desabilitado que, no ponto devolvido de `po::solve`, monta o sistema de Newton-KKT com a `evaluate_kkt` atual, o resolve e compara λ_LS e λ_LS + dλ com a referência `long double` (e, para 1e10, com a de 80 dígitos). Mede também o erro real da Hessiana analítica contra a de 80 dígitos (η) e o resíduo de estacionaridade de (w, λ_LS + dλ) em unidades do piso.
- [ ] **Step 2:** Rodar uma vez. **Critério:** confirma se λ_LS + dλ ≤ ~1e-10 a cond 1e8 (a previsão E3-H da pesquisa); refuta se ficar na ordem do λ atual. **Se refutar, parar e reportar com os números (NEEDS_CONTEXT): as Tasks 1 e 2 dependem disto.**
- [ ] **Step 3:** Registrar a tabela e η em `validation.md`. **Step 4: Commit** — `test: measure the Newton-KKT multipliers of the target-volatility path`

---

### Task 1: Manter o dλ, corrigir β dos presos, salvaguarda

**Files:**
- Modify: `src/riskparity/newton_polish.cpp`, `src/riskparity/newton_polish.hpp`, `src/riskparity/risk_parity.cpp` (devolver λ novo; portão com λ de LS)
- Test: `tests/unit/riskparity/newton_polish_test.cpp`, `tests/unit/riskparity/risk_parity_test.cpp`, `tests/unit/fuzz/risk_parity/risk_parity_reference.hpp` (o checker, ver Review Focus 2)

**Interfaces:**
- Consumes: `newton_system`/`newton_step` (resolvem K por QR Householder; a coluna b é −scale_b·∇h_b, então dλ_b = scale_b·z_b, com scale_b = `hessian_scale` / max|linha|).
- Produces: `PolishResult` passa a carregar `multipliers` (λ_LS + dλ) e `least_squares_multipliers` (para o portão); `po::Result` devolve os primeiros.

Regras (decididas pela pesquisa, §6 item 3):
- No ponto final do polimento, resolver (ou reaproveitar do último `newton_step`) o sistema KKT e somar dλ às igualdades mantidas, **sem mudar w**.
- β_j dos presos: β_j = (∇f − A'λ_novo)_j + (H_L dw)_j, linha j de H_L sobre as livres; verificar sinal como em `settled`.
- Aceitar dλ só se: o QR tem posto completo, ‖dλ‖ ≤ 1e-2·‖λ‖ e uma segunda solução a partir de (w, λ_novo) devolve ‖dλ₂‖ ≤ ‖dλ‖/2. Senão, manter λ_LS (e contar quantas vezes, para o `validation.md`).
- A decisão de aceitar passos de w não usa o dλ final.

- [ ] **Step 1: Testes**
  - `NewtonPolish.KeepsTheNewtonMultiplierCorrection`: nos casos da Task 0 a 1e6 e 1e8, λ_LS + dλ ≤ o teto medido; e λ_LS sozinho continua no nível antigo (prova que o teste discrimina).
  - `NewtonPolish.RejectsACorrectionThatDoesNotContract`: caso construído com ‖dλ₂‖ > ‖dλ‖/2 → mantém λ_LS.
  - `RiskParity.NearMinimumVolatilityKeepsLeastSquaresMultipliers`: s = σ_min(1 + δ), δ ∈ {1e-2, …, 1e-10} continua `OPTIMAL`, passa `expect_risk_parity_kkt`; registra quantas vezes a salvaguarda manteve λ_LS.
  - `RiskParity.TargetReturnOnlyMultipliersUnchanged`: caminho só com `target_return` a 1e8, bits iguais ou erro ≤ 1e-14 (não piorou).
  - `RiskParityReference.KktUsesLeastSquaresMultipliersForTheFloor`: o checker re-ajusta λ por LS no teste para a estacionaridade e compara os multiplicadores do resultado contra a referência.
- [ ] **Step 2:** Rodar, ver falhar. **Step 3:** Implementar. **Step 4:** PASS (`ctest --preset gcc-release` inteiro, mais as baterias de fuzz habilitadas). Se algum teste de valor (que não seja teto de multiplicador) mudar, parar e reportar.
- [ ] **Step 5: Commit** — `feat: keep the Newton multiplier correction of the RISK_PARITY polish`

---

### Task 2: Tetos medidos, caminho com box, QP, documentação

**Files:**
- Modify: `tests/unit/fuzz/risk_parity/risk_parity_fuzz_test.cpp` (`TargetVolatilityMultipliersByCondition`: novos tetos medidos), `docs/guide/validation.md` (substituir a tabela de precisão dos multiplicadores e o "Future work"), `docs/guide/math.md`, `src/riskparity/risk_parity.hpp`, `include/portopt/portopt.h` (a frase de precisão dos multiplicadores: novas faixas; sem mudar assinatura), `CHANGELOG.md`
- Test: o mesmo arquivo de fuzz

- [ ] **Step 1: Testes**
  - `RiskParityFuzz.TargetVolatilityMultipliersByCondition`: tetos novos medidos por cond, cada um com folga 10 sobre o máximo medido (hoje 1e-9/1e-6/1e-1).
  - `RiskParityFuzz.TargetVolatilityBoundMultipliersByCondition`: box `LongOnly` e `Capped` com vol-alvo, β de presos contra a referência `long double`, sinal e complementaridade.
  - `RiskParity.HomogeneityIdentity` estendido a cond 1e6.
  - `DISABLED_MeasureQpMultipliers`: uma medida dos multiplicadores do Goldfarb-Idnani contra a referência estendida a cond 1e8 (item 5 da pesquisa); só registra o número em `validation.md`.
- [ ] **Step 2:** Rodar; medir; fixar os tetos. **Step 3:** Docs: o header diz a nova precisão por cond e mantém "acima de ~1e10 não validado"; `validation.md` troca a tabela e registra a hipótese da Fase 2 como plano futuro, com a condição que a ativaria.
- [ ] **Step 4:** Verificar: `ctest --preset gcc-release` inteiro, pytest, mypy, `check_style.py`.
- [ ] **Step 5: Commit** — `test+docs: measured multiplier accuracy of the target-volatility path`
