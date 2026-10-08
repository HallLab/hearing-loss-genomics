# Fase 2 — quais variantes contam · análise

**Autor:** Andre Rico · **Data:** 2026-10-08
**Cópia de trabalho pessoal, em português.** Não é página de Confluence.
Versão canônica: [`phase_2/results/FINDINGS.md`](../phase_2/results/FINDINGS.md) ·
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

## Os quatro defeitos encontrados

A Fase 2 começou como verificação e virou reconstrução. Quatro defeitos, quatro causas diferentes.

### 1 — três termos de baixo impacto entraram como pLOF

A lista de termos "quebra o gene" do pipeline incluía três que o próprio VEP classifica como
**IMPACT = LOW**:

```
splice_donor_5th_base_variant
splice_donor_region_variant
splice_polypyrimidine_tract_variant
```

São variantes **perto** do ponto de emenda, não nele. Podem atrapalhar, e na maioria das vezes não
atrapalham. Estavam escritas uma a uma no script — não é filtro esquecido, é uma decisão que alguém
tomou.

### 2 — o filtro do SpliceAI nunca foi aplicado

Esse é o grande. O plano de análise dizia: variante de splice só conta como pLOF **se o SpliceAI
passar de 0,2**. O código calculava o escore do SpliceAI e **nunca o usava**.

**63,8% das entradas da máscara pLOF não passavam no próprio critério documentado** —
694.802 de 1.089.876, em todos os 22 cromossomos.

A intenção estava escrita, o número estava calculado, e o `if` não existia.

### 3 — o REVEL foi lido errado · *levantado pela Nikki*

O campo REVEL do VEP não é um número. É uma **lista**, um valor por transcrito:

```
REVEL_score = ".,0.65,0.712,."
```

O código converteu esse texto para número de uma vez só. Texto com vírgulas não vira número — o
resultado é *vazio*, e a variante foi tratada como se não tivesse escore.

```
linhas missense com escore, como foi lido   :   65.806
linhas missense com escore, lendo a lista   :  517.274     (95,3% estavam lá)

variantes passando REVEL ≥ 0,5, como lido   :    5.391
variantes passando REVEL ≥ 0,5, correto     :   20.618
                                    perdidas:   15.227     (73,9%)
```

**Três de cada quatro variantes danosas sumiram da máscara pDM.**

### 4 — genes que não produzem proteína · *levantada pela Nikki*

**1.101 genes de 19.038** na máscara pLOF não codificam proteína — são lncRNA, pseudogenes, RNAs
antisense.

Um teste de perda de função pergunta se **perder a função da proteína** se associa à doença. Esses
genes não têm proteína para perder. O teste roda, dá um número, e o número não pode significar nada.

Lembra do `TMC3-AS1`? Era o gene **nº 1** dos resultados dela, e é um lncRNA. Esse defeito.

---

## O resultado

Todas as quatro máscaras reconstruídas do zero, direto da anotação da release:

| máscara | antes | depois | |
|---|---:|---:|---|
| `pLOF` | 1.002.120 | **462.144** | −54% · saiu o que falhava no SpliceAI |
| `pDM` | 720.983 | **890.332** | **+23%** · voltou o que o REVEL tinha perdido |
| `pLOF_pDM` | 1.717.883 | **1.340.936** | −22% |
| `ALL` | 21.415.507 | **15.951.289** | −26% · saíram os genes não-codificantes |

Repare que o `pDM` **cresceu**. Os outros encolheram porque tiravam coisa que não devia estar; o pDM
cresceu porque recuperou coisa que devia e não estava. Defeitos em direções opostas.

E só três dessas máscaras são usadas: a `ALL` ficou de fora por decisão sua, com o Doug chegando à
mesma posição por conta própria três dias antes — *"why would you ever use the all category... it's
going to give you a lot of noise"*.

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

- **De onde veio a lista de termos de splice?** Os três de baixo impacto estão escritos à mão no
  script. Alguém decidiu incluí-los, e o motivo não está registrado em lugar nenhum. Pergunta para a
  Elena ou a Nikki.
- **REVEL 0,5 ou 0,6?** Usamos 0,5 como primário. É o valor mais comum na literatura, mas é escolha,
  não lei.
- **O AlphaMissense não cobre tudo.** Para genes sem escore, só o REVEL decide. Quantos genes ficam
  nessa situação não foi medido.
- **Sem cromossomo X.** As máscaras cobrem os cromossomos 1 a 22. Num estudo de surdez isso é
  limitação de verdade — existem genes de surdez ligados ao X.
