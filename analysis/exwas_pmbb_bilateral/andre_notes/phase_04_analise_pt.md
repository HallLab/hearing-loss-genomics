# Fase 4 — o teste de associação · análise

Braço: **ExWAS de variante comum, fenótipo bilateral neurossensorial**.
Rodado em 2026-10-09, job LSF `49576771` (66 tarefas) + `49576838` (1 refeita).

---

## O que esta fase faz

Testa cada variante que passou MAF ≥ 1%, uma por uma, contra o desfecho bilateral
neurossensorial. 66 tarefas: 3 coortes × 22 cromossomos.

Sem step 1 — o modelo nulo vem do braço de burden bilateral
(`../gene_burden_pmbb_bilateral/phase_4/results/step1/`). Mesmo fenótipo, mesmas pessoas,
mesmas covariáveis, mesmo parentesco. O modelo nulo não depende do filtro de MAF.

---

## O resultado

| coorte | testadas | barra secundária | menor p | achados |
|---|---:|---:|---:|---:|
| combined | 70.350 | 7,11e-07 | 3,08e-05 | 0 |
| EUR | 74.566 | 6,71e-07 | 9,47e-06 | 0 |
| AFR | 118.049 | 4,24e-07 | 1,07e-06 | 0 |

Nada significativo. λ de 1,0159 / 1,0268 / 1,0130 — bem calibrado.

---

## O ponto que governa a leitura deste braço inteiro

**As duas coortes são aninhadas.** Verifiquei em vez de supor:

- os **3.164** casos bilaterais são todos casos no braço amplo
- os **50.746** controles são **as mesmas pessoas** nos dois
- a restrição remove 3.588 casos e mais nada

Então a análise bilateral é a análise ampla com casos deletados. Três consequências, e
todas mudam como eu posso falar dos números:

### 1. Concordância entre os dois braços não mede nada

Não é replicação. É o mesmo dado duas vezes. Se eu dissesse "nossos dois braços
concordam", estaria dizendo "os dados concordam consigo mesmos".

O Spearman entre eles é 0,388 / 0,398 / 0,302 — e eu reporto só para deixar explícito
que **não** estou usando isso como estatística de replicação.

### 2. Variante que ficou mais fraca não diz nada

Perder 53% dos casos faz isso sozinho. Se uma variante enfraqueceu, não dá para separar
"a associação não era do subtipo bilateral" de "perdemos o poder de ver".

### 3. Variante que ficou MAIS FORTE é a única direção que informa

Essa não pode vir de dado compartilhado. Remover casos remove informação — se o p-valor
melhorou mesmo assim, é porque a associação se concentra no subtipo.

Pelo menos é o que parece. A fase 5 testa isso.

---

## O que ficou mais forte

Critério: mais de uma unidade log mais forte que o braço amplo, **e** p abaixo de 1e-4.

| coorte | quantas | destaque |
|---|---:|---|
| combined | 2 | `chr19:44819634` e `chr19:44819705`, as duas em `BCAM` |
| EUR | 4 | as mesmas duas de `BCAM`, OR 1,17 → 1,35 |
| AFR | 8 | `chr7:158133844`, OR 1,80 → 2,90, o menor p da coorte |

`BCAM` é "basal cell adhesion molecule", no cromossomo 19q13.32 — o grupo sanguíneo
Lutheran. Uma das variantes é missense, a outra sinônima. **Aparece em duas coortes
independentes de análise.**

Na cara, isso é o achado mais interessante do braço. A fase 5 mostra que não é.

---

## Arquivos

| | |
|---|---|
| `phase_4/results/merged_{combined,EUR,AFR}.tsv` | uma linha por variante testada |
| `phase_4/results/MERGE_REPORT.md` | a tabela e os cheques que passaram |
| `phase_4/scripts/01_saige_step2.bsub` | 66 tarefas, step 1 herdado |
| `phase_4/scripts/02_merge_results.py` | recusa com job em voo, cabeçalho divergente, contagem fora do teto do censo, p fora de (0,1] |

---

## Em aberto ao fim desta fase

Uma coisa só, e é grande: **as variantes que fortaleceram são escolhidas entre 70 a 118
mil comparações.** Sem controle, apontar para as que melhoraram é escolher o vencedor
depois da corrida terminar.

É o que a fase 5 resolve.
