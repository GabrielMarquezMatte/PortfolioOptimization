# Binding C# (.NET 10) — design

2026-11-19. Pedido do autor: "o bind de C#, como foi feito no VolEngine", com a mesma motivação (ampliar o alcance e
usar num serviço do autor) e sem mais perguntas. Este design é o do VolEngine
(`VolEngine/docs/superpowers/specs/2026-11-19-volengine-csharp-design.md`) adaptado ao PortOpt; o que difere está
marcado "minha escolha".

## Intenção

Um consumidor .NET da C API, como o Python (ctypes) e o wrapper C++: `include/portopt/portopt.h` é o único contrato.
Nenhuma função nova na C API, nenhuma lógica numérica no C#.

## Metas

- **Cobertura**: as 36 funções do `portopt.abi`, com os conceitos do Python: `Spec` (objetivo, bounds escalares ou por
  peso, target return, risk aversion, importance, max volatility, leverage, target volatility, risk budget),
  `MarketData` (o `Data` do Python: mu, covariância, risk-free, retornos e benchmark; o nome `Data` colide com o
  namespace `System.Data`, CA1724), `Optimizer.Solve` e `Result` (status, pesos, iterações, os quatro
  multiplicadores nomeados, multiplicadores dos bounds, residual), as métricas de `Data` (retorno, volatilidade,
  utilidade, tracking error, Sharpe, diversificação, contribuições de risco) e `NativeAbi` (versão do ABI e da
  biblioteca).
- **Versão do pacote**: a forma SemVer da tag (`release_version.py nuget`: `v0.1.0rc1` e `v0.1.0-rc.1` dão
  `0.1.0-rc.1`); uma execução manual empacota a tag de ensaio (`-rc.0`), nunca uma versão final.
- **Pacote**: `PortOpt` no nuget.org, versão = tag do release, com as 7 bibliotecas que o `wheels.yml` já produz
  (`portopt-native-<rid>.<ext>`) em `runtimes/<rid>/native/libportopt.<ext>`; `[LibraryImport("libportopt")]` resolve
  nos três sistemas. Publicado pelo `release.yml` depois do GitHub Release, por Trusted Publishing (`NuGet/login`,
  `vars.NUGET_USER`); uma falha no NuGet não bloqueia o release nem o PyPI. Antes da primeira tag o autor confere que o
  ID `PortOpt` está livre e registra a policy.
- **Biblioteca local**: `PORTOPT_LIB` (a mesma do Python) via `NativeLibrary.SetDllImportResolver`; sem ela, a
  resolução do NuGet. Na carga, `po_abi_version()` diferente de `PO_ABI_VERSION` (0) lança.
- **AOT/trimming**: `[DisableRuntimeMarshalling]`, tudo blittable, `IsAotCompatible`, warnings como erro; um
  consumidor limpo publicado com Native AOT no CI.
- **API 0.x** (minha escolha): o projeto está em 0.x e quebra o ABI entre versões menores, então a API C# não é
  declarada estável. `PublicApiAnalyzers` fica, com tudo em `PublicAPI.Unshipped.txt`, para que toda mudança da
  superfície apareça no diff; o baseline de `EnablePackageValidation` só entra na 1.0.

## Forma da API

- **Erros**: `PortOptException` com `Status` e a mensagem de `po_last_error`, lida logo após a chamada. Forma errada do
  lado C# (covariância não quadrada, tamanhos que não batem) é `ArgumentException` antes do nativo; uma `Spec` com
  bounds escalares e por peso ao mesmo tempo é `InvalidOperationException` (estado do record, não um argumento).
- **`Spec`**: `sealed record` imutável com propriedades `init`, sem handle; o handle nativo é criado no `Solve` e
  liberado em seguida. `Validate()` roda os setters nativos para quem quer falhar cedo. Bounds como no wrapper C++:
  `Bounds` (um par para todo peso) ou `AssetBounds` (um par por peso, caixa por último); igualdade por valor, também
  dos arrays.
- **`MarketData`**: `IDisposable` sobre um `SafeHandle`; construída de uma vez (como no Python), imutável depois. Covariância
  e retornos como `double[,]` (ordem C, fixados e passados sem cópia). Usar depois do `Dispose` é
  `ObjectDisposedException`.
- **`Result`**: valores copiados (`double[]`), sem handle; válidos para sempre.
- **Threads**: o que o header promete: `po_solve` em paralelo sobre os mesmos spec e data só de leitura.

## Validação

Teste do contrato (cada `function` e `define` do `portopt.abi` contra as declarações e os enums, por reflexão), teste
de cobertura (toda função chamada fora das declarações), testes da camada (erros, `Dispose`, formas, cada objetivo),
paridade bit a bit com chamadas nativas diretas num caso, xUnit v3; `csharp.yml` nos três sistemas, com o pacote
conferido por um consumidor limpo e Native AOT; benchmark BenchmarkDotNet do custo do binding contra chamadas diretas.

## Abordagens

1. Declarações à mão + teste do contrato (a do VolEngine; escolhida).
2. Gerar as declarações do `portopt.abi`: mais um gerador para 36 funções. Descartada.
3. ClangSharp ou C++/CLI: dependência externa ou só Windows. Descartadas.

## Não-metas

`async`, netstandard2.0/.NET Framework, RIDs além dos 7, tipos de terceiros na API, função nova na C API.
