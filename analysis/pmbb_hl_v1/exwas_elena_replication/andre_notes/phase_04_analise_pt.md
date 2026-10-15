# Fase 4 — o teste de associação · análise

Documento em português, linguagem de leigo, para eu conseguir explicar em reunião.
Braço: **ExWAS de variante comum, fenótipo amplo (replicação)**.
Rodado em 2026-10-09, jobs LSF `49569375` (66 tarefas) + `49571242` (1 refeita).

---

## O que esta fase faz

Testa, uma por uma, cada variante que passou o filtro de frequência da fase 2, para ver
se ela aparece mais em quem tem perda auditiva do que em quem não tem.

São 66 tarefas: 3 coortes × 22 cromossomos. Cada tarefa roda o SAIGE sobre um cromossomo
de uma coorte e devolve uma linha por variante testada.

---

## Conceito 1 — por que não existe step 1 aqui

Nos braços de variante rara o SAIGE roda em dois passos. O **passo 1** monta o modelo de
base: a expectativa de quem tem perda auditiva a partir de idade, sexo, batch,
ancestralidade e parentesco — **sem olhar gene nenhum**. O **passo 2** compara o
observado contra essa expectativa.

O passo 1 é a parte caríssima, por causa do parentesco: 1,45 bilhão de pares de pessoas
na coorte combinada.

E ele **não depende do filtro de MAF**. São as mesmas pessoas, o mesmo desfecho, as
mesmas covariáveis. Rodar de novo daria exatamente os mesmos números.

Então herdamos:

```
../gene_burden_elena_replication/phase_4/results/step1/{combined,EUR,AFR}/
```

E isso é o que a **Elena** fez, não invenção nossa. Na config dela:

```
step1_outdir = "${base}/HL_only_rarevariant/saige_step1"
```

O ExWAS dela leu o passo 1 do braço de variante rara.

---

## Conceito 2 — o que muda em relação ao burden

| | burden (variante rara) | ExWAS (aqui) |
|---|---|---|
| unidade testada | o **gene** | a **variante** |
| quantos testes | ~18.000 genes × 9 células | ~70 a 118 mil variantes |
| precisa de máscara? | sim — define quem entra no gene | não — não há agrupamento |
| precisa de VEP? | sim, para construir o teste | não, só para interpretar depois |
| combinação de Cauchy | sim, junta as 9 células | não existe |

Na prática, o comando do SAIGE perde quatro parâmetros: `--groupFile`,
`--maxMAF_in_groupTest`, `--annotation_in_groupTest` e as opções de peso. **Tirar o
`--groupFile` é o que coloca o SAIGE em modo de variante única** — nada mais seleciona
isso.

---

## O resultado: nada significativo

| coorte | testadas | barra secundária | menor p | achados em 5e-8 | achados na secundária |
|---|---:|---:|---:|---:|---:|
| combined | 70.366 | 7,11e-07 | 3,02e-06 | 0 | 0 |
| EUR | 74.536 | 6,71e-07 | 2,12e-05 | 0 | 0 |
| AFR | 117.953 | 4,24e-07 | 4,11e-06 | 0 | 0 |

O melhor resultado de todos está **4,3× longe** da barra mais frouxa e **60× longe** da
primária.

---

## Mas a replicação em si é a melhor do projeto

Juntando por `MarkerID` com os arquivos entregues dela:

| coorte | dela | nossa | em comum | correlação de −log10(p) |
|---|---:|---:|---:|---:|
| combined | 70.429 | 70.366 | 70.333 (99,9%) | **0,989** |
| EUR | 74.528 | 74.536 | 74.483 (99,9%) | **0,979** |
| AFR | 117.943 | 117.953 | 117.723 (99,8%) | **0,990** |

Os mínimos dela — 2,57e-06, 1,99e-05, 3,70e-06 — estão a um fator de 1,2 dos nossos. E a
coluna `Significance` dos arquivos dela diz **"Not significant"** em todos.

As variantes que aparecem em só um dos lados são 33 a 230 por coorte, e se explicam pela
diferença de coorte: nossa versão corrigida tem 40 casos a mais e 174 controles a menos,
o que move um punhado de variantes para os dois lados do limiar de 1%.

As maiores discordâncias individuais estão todas no **MHC** (chr6:31–32 Mb), onde o LD é
extremo e nossos PCs diferem dos dela (4/3 contra 9/10 dela).

---

## O pôster diz que há seis resultados significativos

Em `elena_publishes/Horizontal Draft LANS Poster - SUIP-PGSIP.pptx`, sob *"Genome Wide
Significant Results Total: 6"*:

| gene | variante | coorte | p |
|---|---|---|---|
| KNL1 | chr15:40611486:A:G | Combined | 3,51e-12 |
| KNL1 | chr15:40611486:A:G | EUR | 7,21e-12 |
| ZFPM1 | chr16:88533289:A:C | EUR | 2,46e-10 |
| ZFPM1 | chr16:88533290:G:C | EUR | 2,52e-10 |

Três variantes únicas. Nenhuma delas está nos arquivos de genótipo que o pipeline usa,
na saída bruta do SAIGE dela, nas tabelas anotadas dela, nem nos nossos resultados.

### Onde elas estão: no release do PMBB

E aí está a resposta. As três **existem** no `PMBB-Release-2026-4.0`. O que elas não
fazem é passar o filtro que **o próprio pôster descreve**:

> *"variants with a genotype call rate < 0.99 and a Minor Allele Frequency (MAF) < 0.01
> were removed."*

Medido direto no release:

| variante | faltante | call rate | MAF | p no pôster |
|---|---:|---:|---:|---:|
| `chr15_40611486_A_G` | 56.646 / 70.925 | **0,2013** | 0,284 | 3,51e-12 |
| `chr16_88533289_A_C` | 2.717 / 70.925 | 0,9617 | 0,295 | 2,46e-10 |
| `chr16_88533290_G_C` | 3.736 / 70.925 | 0,9473 | 0,283 | 2,52e-10 |

As três estão abaixo de 0,99. E a de `KNL1`, a que carrega o `3,51e-12`, tem genótipo em
**20% da coorte** — 79,9% das pessoas não têm chamada naquela posição.

---

## Conceito 3 — por que call rate baixo fabrica p-valor

Isto é a parte que eu preciso saber explicar.

"Call rate" é a fração de pessoas em que a máquina conseguiu ler o genótipo naquela
posição. Call rate 0,20 quer dizer que em 4 de cada 5 pessoas o dado simplesmente não
existe.

O problema não é ter pouco dado — é que **a falha de leitura não é aleatória**. Ela
depende de qualidade de amostra, de lote de sequenciamento, de quando a pessoa entrou no
biobanco. E essas coisas se correlacionam com quem tem diagnóstico no prontuário: quem
frequenta mais o hospital tem mais código registrado, e pode ter sido sequenciado em
outro lote.

Então se os casos têm 70% de leitura e os controles 15%, a variante vai parecer
diferente entre os dois grupos — e o teste estatístico não tem como saber que a
diferença é da máquina e não da biologia. Ele devolve um p-valor minúsculo com toda
confiança.

**É exatamente para isso que o filtro de call rate existe.** Ele não é burocracia de
pipeline; é o que separa sinal de artefato de genotipagem.

Um jeito de enxergar: um p-valor de `3,51e-12` em cima de 20% de dado é uma afirmação
muito mais forte do que o dado permite sustentar. O resto do exoma inteiro, com dado
completo, não chega perto de `5e-8`.

---

## Onde cada variante caiu, estágio por estágio

| estágio | variantes no chr16 | as do pôster |
|---|---:|---|
| release PMBB v4 | 780.404 | presentes |
| `variant_qc/chr16_callrate.bim` | 670.730 | **ainda presentes** |
| `plink_deduplicated/chr16_allele_ids.bim` | 662.670 | removidas |
| `chr16_deduplicated.bim` | 662.667 | removidas |

E um detalhe que vale guardar: aplicar um filtro **real** de 0,99 ao release deixa
**475.376** variantes no chr16. O arquivo dela que leva o nome `callrate` tem
**670.730**. Ou seja, o estágio batizado de call rate não está filtrando em 0,99 — é bem
mais permissivo. Foi por ali que passaram as de 0,96 e 0,95. A do chr15, com 0,20, foi
removida nesse estágio mesmo; as duas do chr16 caíram só depois.

---

## Como eu levo isso para a Nikki

**O que NÃO dizer:** "nosso resultado é diferente do pôster". Isso coloca o nosso
pipeline contra o dela e convida uma defesa.

**O que dizer:** os arquivos que ela nos entregou **concordam com os nossos a r = 0,99**,
e a coluna `Significance` dos arquivos dela já diz "Not significant" em tudo. Quem
discorda do pôster é o pipeline dela, não o nosso. A gente só foi ver por quê.

**Por onde começar:** pelo número que não depende de pipeline nenhum — `chr15_40611486`
tem 79,9% de genótipo faltante **no release do PMBB**. Isso é medido na fonte. Não é
nosso cálculo, não é resultado de análise, não há o que replicar ou contestar.

**A pergunta concreta:** qual conjunto de variantes gerou a tabela do pôster, já que o
método declarado removia call rate abaixo de 0,99 e essas três estão em 0,20, 0,96 e
0,95?

**Duas coisas que mudam o tom:**

- O arquivo se chama `Horizontal Draft LANS Poster` — é **rascunho**.
- **Eu sou coautor listado**, junto com a Nikki, a Molly e o Doug. Não é reportar o
  trabalho de outra pessoa: é o meu nome também. Isso me dá tanto a legitimidade quanto
  a obrigação de levantar.

**E a cronologia é inocente:** o pôster saiu no verão, de uma rodada que ainda tinha
variantes de baixo call rate. O QC foi apertado depois. Os arquivos entregues são de 2 de
outubro, subidos no dia da reunião, e não têm nada significativo. É a sequência normal de
uma correção, não de um problema de conduta.

---

## Um erro meu nesta fase, registrado

Eu pus na tarefa um cheque que exigia que a contagem de linhas do SAIGE fosse **igual**
ao censo do plink2 da fase 2. Quarenta e uma das 66 células falharam.

O cheque estava errado de desenho. Os dois programas calculam MAF com denominadores
diferentes:

| | denominador |
|---|---|
| plink2 | `ALT / observados` — só quem teve leitura |
| SAIGE | `ALT / 2N` — a coorte inteira |

Para variante com dado faltante, a MAF do plink2 é a maior das duas. Então o plink2
conta como ≥1% variantes que o SAIGE exclui — e o contrário é impossível. O censo é um
**teto**, não uma igualdade.

Verifiquei em AFR/chr1, a célula com a maior diferença (−6): todas as 6 variantes
divergentes atravessavam o limiar exatamente como previsto, e **zero** iam no outro
sentido. Nas 66 células, as 41 divergências apontam todas para o mesmo lado.

Efeito colateral útil: **o piso do SAIGE é mais estrito que 1% nominal** para qualquer
variante com dado faltante. Com 5% de faltante ele exige MAF observada de ~1,053%. O
conjunto testado é um pouco menor e um pouco mais limpo do que o censo sugeria.

O cheque agora testa o teto com margem (recusa se passar do teto, ou se ficar mais de
0,5% abaixo).

---

## Em aberto

- **A pergunta para a Nikki**, acima. É o único item que não é nosso para resolver.
- **A anotação dos achados** (premissa P9) não tem o que anotar: não houve achado. Fica
  registrado que a etapa existia e não foi acionada por falta de alvo, não por escolha.
- **A premissa P7** foi reescrita. Ela nomeava `KNL1` e `ZFPM1` como alvos declarados
  para acertar ou errar; não existe alvo nos dados.
