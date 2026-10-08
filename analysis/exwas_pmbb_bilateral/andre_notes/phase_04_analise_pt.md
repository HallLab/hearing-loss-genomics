# Fase 4 — a pergunta, finalmente · análise

**Autor:** Andre Rico · **Data:** 2026-10-08
**Cópia de trabalho pessoal, em português.** Não é página de Confluence.
Resultados: [`phase_5/results/FINDINGS.md`](../phase_5/results/FINDINGS.md) ·
Premissas: [`PREMISES.md`](../PREMISES.md) · Fase anterior:
[`phase_03_analise_pt.md`](phase_03_analise_pt.md)

---

## O que esta fase faz

As três primeiras montaram peças. Nenhuma respondeu nada.

| fase | montou |
|---|---|
| 1 | quem está no estudo — 3.164 casos, 50.746 controles |
| 2 | quais variantes contam — três máscaras |
| 3 | o que descontar — idade, sexo, lote, ancestralidade |

Esta é onde a pergunta é feita:

> **Gene por gene: quem tem perda auditiva bilateral neurossensorial carrega mais mutações raras
> estragadoras nesse gene do que quem não tem?**

É a primeira fase que produz p-valores.

---

## Conceito 1 — por que somar as mutações por gene

Mutação rara é rara mesmo. Uma variante específica pode aparecer em 3 pessoas das 53.910 — e com 3
pessoas não se conclui nada, nunca.

A saída é não testar variante por variante, e sim **somar por gene**. Em vez de "esta mutação está
associada?", a pergunta vira "**qualquer** mutação estragadora neste gene está associada?".

As 3 pessoas de uma variante, mais 5 de outra, mais 2 de outra, viram 40 pessoas com algum dano
naquele gene. E aí dá para testar. Isso é o **burden**: a carga de dano da pessoa naquele gene.

É também por isso que a Fase 2 importava tanto — se a lista de "mutações estragadoras" estiver
errada, você soma as coisas erradas.

## Conceito 2 — por que o SAIGE trabalha em dois passos

O trabalho se separa em duas peças de natureza diferente:

**O que não depende do gene** — idade, sexo, lote, ancestralidade e, sobretudo, **quem é parente de
quem**. Idêntico quer você teste o primeiro gene ou o décimo oitavo milésimo.

**O que depende do gene** — as mutações raras naquele gene.

---

### Passo 1 — a linha de base

Constrói a expectativa de quem tem perda auditiva **sem olhar nenhum gene**. Roda uma vez por
coorte: três vezes no total.

**Não é um modelo por pessoa — é um modelo só, sobre as 53.910.** Ele *produz* um número por pessoa,
mas o modelo é um.

A fórmula que o SAIGE registrou no log:

```
PHENO ~ AGE + AGE2 + SEX + Batch + PC1 + PC2 + PC3 + PC4 + PC5
```

Só que isso é metade. O modelo completo é uma regressão logística **mista** — o `MM` de
`fitNULLGLMM`, de *Mixed Model*:

```
                    ┌──────────── efeitos fixos ─────────────┐   ┌── aleatório ──┐
  logit( P(casoᵢ) ) = β₀ + β₁·AGEᵢ + β₂·AGEᵢ² + β₃·SEXᵢ + ... +         bᵢ

  com    b ~ Normal( 0 , τ · GRM )
```

**O `b` é o parentesco.** Um valor por pessoa, e pessoas aparentadas têm valores **correlacionados**
— a matriz GRM, aquela de 1,45 bilhão de pares, é quem diz o quanto. Sem esse termo, dois irmãos
entrariam no modelo como se fossem dois indivíduos independentes, e a semelhança genética entre eles
viraria sinal falso.

O **`τ`** (tau) é o tamanho desse efeito: quanto do fenótipo o parentesco explica. O SAIGE o aprende
por iteração, e dá para ver acontecendo no log:

```
inital tau is  1 0
tau:  1 0.1040913
tau:  1 0.1544189
tau:  1 0.1553402      ← convergiu
```

É o `tau₂` que aparece na tabela de resultados. **Quando ele converge para zero — como no AFR — o
termo aleatório desaparece e o modelo *vira* a regressão logística comum** que você descreveu.

### O que realmente atravessa para o passo 2: os resíduos

Aqui está a ideia que faz os dois passos serem uma coisa só, e não duas etapas administrativas.

O que interessa não é o ajuste — é o quanto cada pessoa **escapa** dele:

```
resíduoᵢ  =  yᵢ − μ̂ᵢ        y = 1 se caso, 0 se controle
                              μ̂ = probabilidade que o modelo previu
```

| | o modelo previu | é caso? | resíduo | leitura |
|---|---|---|---|---|
| pessoa A | 8% | sim | **+0,92** | muito surpreendente |
| pessoa B | 40% | sim | +0,60 | meio esperado |
| pessoa C | 35% | não | −0,35 | esperado |

A pessoa A é **surpreendente**: idade, sexo, ancestralidade e parentesco dela não explicam o
diagnóstico. Sobrou algo por explicar.

---

### Passo 2 — a pergunta por gene

E a pergunta vira uma só:

> **As variantes raras deste gene estão concentradas nas pessoas surpreendentes?**

Formalmente, para a soma (*burden*), colapsa-se as variantes do gene num número por pessoa e
calcula-se:

```
  Gᵢ = Σⱼ wⱼ · genótipoᵢⱼ          a carga de dano da pessoa i naquele gene

  S  = Σᵢ Gᵢ · ( yᵢ − μ̂ᵢ )          o escore
       └─────────┬──────────┘
         carga × surpresa, somado sobre todo mundo
```

Se quem carrega o gene é justamente quem surpreende, os dois fatores são grandes juntos e `S` fica
longe de zero. Se o gene está espalhado igualmente entre surpreendentes e esperados, os termos
positivos e negativos se cancelam e `S` fica perto de zero.

Para a dispersão (*SKAT*) a conta é a mesma, mas sem colapsar e elevando ao quadrado, o que mata o
sinal de menos:

```
  Sⱼ = Σᵢ genótipoᵢⱼ · ( yᵢ − μ̂ᵢ )      um escore por variante

  Q  = Σⱼ wⱼ² · Sⱼ²                     soma dos quadrados
```

Em ambos, **o `(yᵢ − μ̂ᵢ)` vem pronto do passo 1**. É literalmente a única coisa do modelo nulo que o
passo 2 precisa — e é por isso que ele não tem que reajustar nada.

### O variance ratio, o segundo truque

Calcular a **variância** de `S` exigiria a matriz de parentesco inteira em cada gene — caro de novo.
Então o SAIGE faz uma conta barata, com uma versão esparsa, e corrige:

```
  Var(S)  ≈  r × Var_barata(S)
```

Esse `r` é o **variance ratio** — estimado uma vez no passo 1, a partir de marcadores sorteados, e
reusado em todos os genes. Os nossos: 0,9847, 0,9942 e 1,0000. Perto de 1 significa que a conta
barata quase não erra.

---

### E há uma razão que não é só economia

A linha de base precisa ser construída **sem nunca ter visto o gene**. Se você ajustasse as
covariáveis e o gene ao mesmo tempo, eles competiriam pela mesma variação, e o modelo poderia
absorver parte do sinal real dentro das covariáveis. Separando, a linha de base fica congelada antes
de qualquer gene ser perguntado — e **todos os genes são medidos com a mesma régua**.

Uma analogia: para avaliar se uma escola ensina bem, primeiro você monta a previsão de nota de cada
aluno a partir do que **não** tem a ver com a escola — renda, notas anteriores. A diferença entre a
nota real e a previsão é o resíduo. Depois pergunta, escola por escola, se os alunos dela escapam da
previsão para cima. Montar a previsão uma vez e reusá-la em todas é o que torna a conta viável.

**Consequência prática:** os 66 trabalhos do passo 2 se apoiam nos três modelos nulos do passo 1.
Um erro no passo 1 contamina tudo; não existe gene que escape.

### O que atravessa, em arquivo

```
combined.rda                      372 MB   o modelo ajustado, com os resíduos
combined.varianceRatio.txt         57 B    um número: 0,9847
```

372 megabytes e um número, calculados uma vez em ~15 minutos, lidos pelas 22 tarefas daquela coorte.

---

## Conceito 3 — qual é o teste

Não é razão de verossimilhança. O SAIGE usa um **teste de score**.

| teste | o que faz | ajustes de modelo |
|---|---|---|
| razão de verossimilhança | ajusta sem o gene e com o gene, compara | 2 |
| Wald | ajusta com o gene, olha efeito sobre erro | 1, com o gene |
| **score** | ajusta **só sem** o gene, e pergunta "para que lado ele gostaria de se mover?" | 1, sem o gene |

O score só precisa do modelo do passo 1 — do resíduo `(yᵢ − μ̂ᵢ)`, como no Conceito 2. É literalmente
o que torna a estrutura de dois passos possível. E tem uma segunda vantagem: com dado raro ele
continua **válido** onde os outros quebram.
Um gene com 5 alelos todos em casos é separação perfeita; o Wald e o LRT precisariam estimar um
efeito que vai para infinito. O score nunca estima nada sob a alternativa.

**O "SPA" do nome** é *saddlepoint approximation*. O teste de score padrão supõe curva normal, o que
falha **na ponta** quando casos e controles estão desbalanceados e a variante é rara. Temos os dois:
3.164 casos contra 50.746 controles, com variantes de poucos alelos. Sem corrigir, os p-valores
sairiam otimistas demais, e exatamente nos genes mais raros.

---

## Conceito 4 — os quatro p-valores por gene

As estatísticas estão no Conceito 2; aqui é o que cada uma **responde**, e quando cada uma falha.

**`Pvalue_Burden` — a soma** (`S = Σᵢ Gᵢ · resíduoᵢ`).
Forte se todas as mutações do gene empurram para o mesmo lado. **Frágil se metade é danosa e metade
protetora** — na soma elas se cancelam e o teste não vê nada.

**`Pvalue_SKAT` — a dispersão** (`Q = Σⱼ wⱼ² · Sⱼ²`).
Elevar ao quadrado mata o sinal de menos, então direção deixa de importar: danosas e protetoras
contribuem igual. **Sobrevive a direções misturadas; perde força quando todas apontam para o mesmo
lado**, porque jogou fora a informação de direção.

O "kernel" do nome (*Sequence Kernel Association Test*) é só uma forma de medir **o quanto duas
pessoas se parecem** naquele gene. O SKAT pergunta: quem se parece aqui também se parece no
fenótipo?

**`Pvalue` — o SKAT-O, a mistura.** Os dois são casos extremos de uma mesma família:

```
  Q(ρ)  =  (1 − ρ) · Q_SKAT  +  ρ · Q_Burden            ρ entre 0 e 1

     ρ = 0  →  SKAT puro      (só dispersão)
     ρ = 1  →  Burden puro    (só soma)
```

O SKAT-O varre os valores de `ρ`, pega o melhor, e **corrige o p-valor por ter varrido** — senão
escolher o melhor de vários seria trapaça. Resultado: funciona nos dois cenários, sem precisar
adivinhar antes qual deles é o seu gene.

**`Cauchy` — o omnibus.** Este é o principal, e é o que esta análise reporta. Enquanto o SKAT-O
combina **duas estatísticas dentro de uma célula**, o Cauchy combina **as nove células** — é a
seção seguinte.

---

## O desenho que torna o Cauchy possível

Cada gene é testado **nove vezes**: três máscaras × três cortes de raridade. Saem nove p-valores, e
aí vem o problema: qual é *o* p-valor daquele gene?

Pegar o menor dos nove e comparar com um limiar calculado como se fosse um teste só é **cobrar
barato** — você garimpou nove e pagou por um.

A **combinação de Cauchy** resolve. Ela transforma cada p-valor antes de somar, de um jeito em que
**evidência forte grita e ausência de evidência sussurra**:

```
p = 0,000001  →  318.310        p = 0,5  →      0
p = 0,0001    →    3.183        p = 0,9  →     -3
```

A média fica dominada pelo menor p-valor, e o resultado segue uma regra simples:

```
omnibus ≈ menor p × ( 9 ÷ em quantas das 9 células o sinal aparece )

   sinal em 1 célula  →  9,0×   exatamente o preço de Bonferroni
   sinal nas 9        →  1,0×   nenhuma penalidade
```

Ou seja: **"pegue a melhor célula, e pague pelo garimpo que deu para achá-la."**

### Uma chamada ao SAIGE em vez de nove

O SAIGE faz essa combinação **sozinho** — mas só se as nove combinações estiverem na mesma execução.
Os dois parâmetros aceitam lista:

```
--annotation_in_groupTest="pLOF,pDM,pLOF:pDM"
--maxMAF_in_groupTest=0.0001,0.001,0.01        ← este é até o padrão do SAIGE
```

Então rodamos **uma chamada por coorte e cromossomo**: 3 × 22 = **66 tarefas**. Cada uma testa as
nove combinações internamente e emite **uma linha `Cauchy` por gene**, cobrindo todas.

O número principal vem da ferramenta, não de um cálculo meu rodando depois.

> Esse desenho foi testado antes de ser adotado, comparando uma chamada contra nove no chr21:
> as células compartilhadas batem 100%, e o Cauchy nativo bate com o cálculo manual em 165 de
> 208 genes. Detalhes na premissa **P10** do [`PREMISES.md`](../PREMISES.md).

---

## O que rodou, e o que saiu

### Passo 1 — três modelos nulos

| coorte | amostras | variance ratio | tau₂ |
|---|---:|---|---:|
| combined | 53.910 | 0,9847 | 0,155 |
| EUR | 40.143 | 0,9942 | 0,059 |
| AFR | 10.539 | 1,0000 | **0** |

O **tau₂** é a variância poligênica — quanto do fenótipo o parentesco explica. No AFR ele converge
para **zero**, o que significa que ali o SAIGE vira regressão logística com covariáveis: com 490
casos não há componente poligênico detectável. Não é defeito; é o tamanho da amostra aparecendo.

O **variance ratio** é o fator que corrige o cálculo barato que o SAIGE faz em cada gene para o que
o cálculo caro daria. Perto de 1 é o esperado.

### Passo 2 — 66 tarefas, nenhuma falha

Cada arquivo cobre um cromossomo inteiro. O chr1 da coorte combinada, por exemplo:

```
pLOF       5.795 linhas     ┐
pDM        5.687 linhas     │ as nove células
pLOF;pDM   5.844 linhas     ┘
Cauchy     1.948 linhas     ← uma por gene. É o número principal.
```

Tempo: cerca de 36 minutos por tarefa na coorte combinada.

**O resultado está na Fase 5** — [`phase_5/results/FINDINGS.md`](../phase_5/results/FINDINGS.md).
Adianto o essencial: **nenhum gene atinge significância em nenhuma coorte.**

---

## Em aberto

- **O GRM do combined é magro** — 38.833 marcadores contra 114.694 do AFR, por causa do filtro de
  Hardy-Weinberg aplicado entre ancestralidades. Herdado da Fase 3 e não resolvido.
- **O AFR com 490 casos e tau₂ = 0.** Rodou porque não rodar é viés de reporte próprio, mas tem que
  ser lido como nulo quase sem poder.
- **Sem cromossomo X**, como nas máscaras.
- **Os genótipos do passo 2 vêm de fora desta pasta** — 438 GB que não dá para copiar. É o único
  caminho que aponta para fora, e está registrado no [`PROVENANCE.md`](../PROVENANCE.md).
