# Fase 2 — quais variantes contam · análise

**Autor:** Andre Rico · **Data:** 2026-10-08
**Cópia de trabalho pessoal, em português.** Não é página de Confluence.
Evidência herdada: [`phase_2/results/FINDINGS.md`](../phase_2/results/FINDINGS.md) (documento da replicação, copiado) ·
Premissas: [`PREMISES.md`](../PREMISES.md) · Fase anterior:
[`phase_01_analise_pt.md`](phase_01_analise_pt.md)

---

## O que esta fase faz

A Fase 1 separou **pessoas**. Esta separa **variantes**.

Cada participante tem cerca de 20 mil variantes raras no exoma. A esmagadora maioria não faz nada —
troca uma letra que não muda a proteína, ou cai num pedaço do gene que não é usado. Se você somar
todas, o sinal de verdade afunda no ruído.

Então antes de testar, é preciso decidir **quais variantes contam**. Essa decisão vira uma
**máscara**.

> **Esta fase não foi refeita aqui.** As máscaras vêm da `elena_replication`, copiadas prontas e
> conferidas byte a byte. Elas dependem da anotação do genoma e **não dependem de quem é caso**, então
> trocar o fenótipo não as invalida. O que está descrito abaixo foi o trabalho feito lá.

---

## Conceito 1 — o que é uma máscara

Uma máscara é uma lista: **para cada gene, quais variantes entram no teste**.

```
A1BG   var    19:58347031:T:C   19:58347349:T:C   19:58347352:C:A   ...
A1BG   anno   pLOF              pLOF              pLOF              ...
```

Duas linhas por gene: as variantes, e o rótulo de cada uma. O SAIGE lê isso e soma só o que está ali.

**Mudar a máscara muda o resultado**, porque muda o que está sendo somado. Por isso ela é tão
decisiva quanto o fenótipo.

## Conceito 2 — as duas formas de estragar um gene

As máscaras do estudo se dividem por **mecanismo de dano**:

**`pLOF` — *predicted Loss of Function*, perda de função prevista.**
A variante **quebra** a proteína. Põe um ponto final no meio (`stop_gained`), embaralha a leitura a
partir dali (`frameshift`), ou destrói o ponto de emenda onde o RNA é costurado
(`splice_acceptor`, `splice_donor`). O resultado é uma proteína truncada ou nenhuma proteína.

**`pDM` — *predicted Damaging Missense*, missense danosa prevista.**
A variante troca **um** aminoácido por outro. A proteína continua inteira, mas pode funcionar mal. A
maioria das trocas é inofensiva — então aqui entra a previsão: um algoritmo estima se *aquela* troca
é danosa.

**`pLOF_pDM`** é a união das duas, e **`ALL`** é todas as variantes raras do gene sem filtro nenhum.

## Conceito 3 — os três preditores

Como se decide se uma variante é danosa? Com ferramentas que o campo construiu para isso:

| ferramenta | o que prevê | escala |
|---|---|---|
| **SpliceAI** | se a variante atrapalha a costura do RNA | 0 a 1 · usamos ≥ 0,2 |
| **REVEL** | se uma troca de aminoácido é danosa | 0 a 1 · usamos ≥ 0,5 |
| **AlphaMissense** | idem, do DeepMind | classes: *pathogenic*, *benign* |

O **VEP** (Variant Effect Predictor) é quem junta tudo: recebe a variante, devolve em que gene ela
cai, que tipo de consequência tem, e os escores dessas ferramentas. A release do PMBB já vem com a
anotação do VEP pronta — é dela que as máscaras são construídas.

---

## De onde vêm os dados

```
/static/PMBB/PMBB-Release-2026-4.0/Exome/   anotação VEP da release
```

**165.369.501 linhas** de anotação foram lidas para reconstruir as máscaras. O número é grande
porque **uma variante aparece várias vezes**: um gene tem vários transcritos, e a mesma variante
pode ser `stop_gained` num e `missense` noutro. Esse detalhe é a causa de um dos defeitos abaixo.

---

## As regras, exatamente como estão no código

Quatro regras, uma por máscara. Estão no
[`phase_2/scripts/07_rebuild_masks.py`](../phase_2/scripts/07_rebuild_masks.py) e são curtas o
bastante para caber aqui inteiras.

### `pLOF` — quebra a proteína

Uma variante entra se **qualquer** uma das duas condições vale:

```
1.  a consequência é exatamente um destes sete termos:

      stop_gained            põe um ponto final no meio
      frameshift_variant     embaralha a leitura dali em diante
      splice_acceptor_variant  ┐ destrói o ponto onde o RNA é costurado
      splice_donor_variant     ┘
      start_lost             apaga o códon de início
      stop_lost              apaga o códon de parada
      transcript_ablation    elimina o transcrito inteiro

2.  OU a consequência menciona "splice" de qualquer outra forma
    E o SpliceAI é >= 0,2
```

A segunda condição é o **portão do SpliceAI**, e é o que separa as variantes de splice que
importam das que só ficam perto. Há vários termos de splice de impacto baixo — `splice_donor_region`,
`splice_polypyrimidine_tract` e outros — que descrevem variantes **próximas** ao ponto de emenda.
Elas podem atrapalhar a costura, e na maioria das vezes não atrapalham. O SpliceAI é quem decide
caso a caso, e só passa quem ele aponta.

**Correspondência é exata, não por substring.** Procurar `"splice_donor_variant"` dentro do texto da
consequência também casaria com `splice_donor_region_variant`, que é outra coisa. A comparação é
termo a termo.

### `pDM` — troca um aminoácido de forma provavelmente danosa

```
o AlphaMissense classifica como "pathogenic" ou "likely_pathogenic"
OU
o maior REVEL entre os transcritos é >= 0,5
```

**"o maior entre os transcritos"** é a parte que importa. O VEP devolve o REVEL como uma **lista**,
um valor por transcrito, porque a mesma variante pode ser missense num transcrito e outra coisa
noutro:

```
REVEL_score = ".,0.65,0.712,."
```

A regra lê a lista e pega o maior valor. Tratar esse texto como um número só devolve vazio, e a
variante some da máscara.

### `pLOF_pDM` — a união

As duas juntas. Quando uma variante satisfaz as duas regras — acontece com **11.517** delas, por
serem perda de função num transcrito e missense noutro — ela entra rotulada como `pLOF`. O rótulo
mais forte vence.

### `ALL` — tudo

Toda variante do gene presente no conjunto de genótipos, sem filtro de dano.

### E uma condição que vale para as quatro

```
BIOTYPE == protein_coding
```

O gene precisa codificar proteína. Um teste de perda de função pergunta se **perder a função da
proteína** se associa à doença — num lncRNA ou pseudogene não há proteína para perder, e o teste
roda, dá um número, e o número não significa nada.

Isso retira **1.101 genes de 19.038**. Aplicamos também à máscara `ALL`, que a rigor não faz
afirmação sobre dano: um gene sem proteína não é candidato a nada neste desenho, então deixá-lo em
qualquer máscara só gasta correção múltipla.

---

## Por que essas regras e não outras

Elas não foram inventadas aqui. São as regras do plano de análise original do estudo, **aplicadas
como escritas** — e a Fase 2 da `elena_replication` existe porque a implementação que rodou se
afastava delas em quatro pontos.

A história completa, com os quatro defeitos, o tamanho de cada um e quem achou, está em
[`elena_replication/docs/09_fase_2_resumo.pt.md`](../../elena_replication/docs/09_fase_2_resumo.pt.md)
— e em detalhe no
[`05_phase_2_fechamento.pt.md`](../../elena_replication/docs/05_phase_2_fechamento.pt.md).

Aqui basta saber que as máscaras desta análise seguem as regras acima, e que isso foi verificado.

---

## O resultado

As quatro máscaras, construídas do zero a partir da anotação da release:

| máscara | genes | variantes |
|---|---:|---:|
| `pLOF` | 17.841 | 462.144 |
| `pDM` | 17.623 | 890.332 |
| `pLOF_pDM` | 17.945 | 1.340.936 |
| `ALL` | 18.038 | 15.951.289 |

Uma leitura que ajuda a calibrar a escala: são ~18 mil genes com proteína no exoma humano, e
praticamente todos aparecem nas três primeiras máscaras. O que muda entre elas não é *quais genes*,
e sim **quantas variantes cada gene leva para o teste** — 26 por gene na `pLOF`, 51 na `pDM`, 75 na
união.

A `ALL` tem 884 variantes por gene, e é por isso que ela não é usada. Decisão sua, e o Doug chegou
à mesma por conta própria três dias antes:

> *"why would you ever use the all category... it's going to give you a lot of noise"*

**Esta análise roda as três primeiras.** A `ALL` fica no disco, para a decisão ser reversível sem
reconstruir nada.

---

## Por que isto está copiado e não refeito

As máscaras dependem de **duas coisas**: a anotação do genoma, e as regras de dano. Nenhuma das duas
muda quando o fenótipo muda. Um gene é um gene; uma variante que quebra a proteína quebra a proteína,
independente de quem foi diagnosticado.

Então reconstruir aqui daria **exatamente os mesmos arquivos**, gastando horas de máquina. Foram
copiadas e conferidas byte a byte — está em [`PROVENANCE.md`](../PROVENANCE.md).

O que **não** é reaproveitável, e por isso foi refeito aqui: coorte, covariáveis, modelos nulos e o
teste. Tudo que toca em quem é caso.

---

## Em aberto

- **O portão do SpliceAI em 0,2.** É o valor do plano de análise e é defensável, mas é escolha.
  Mover para 0,5 tornaria a `pLOF` bem mais estrita, e ninguém mediu o quanto.
- **REVEL 0,5 ou 0,6?** Usamos 0,5 como primário. É o valor mais comum na literatura, mas é escolha,
  não lei.
- **O AlphaMissense não cobre tudo.** Para genes sem escore, só o REVEL decide. Quantos genes ficam
  nessa situação não foi medido.
- **Sem cromossomo X.** As máscaras cobrem os cromossomos 1 a 22. Num estudo de surdez isso é
  limitação de verdade — existem genes de surdez ligados ao X.
