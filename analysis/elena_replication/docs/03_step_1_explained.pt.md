# Step 1, explicado — quem entra no estudo, e como 427 pessoas ficaram de fora

**Autor:** Andre Rico · **Data:** 2026-09-30 · **Tipo:** explicativo
**Cópia de trabalho pessoal, em português.** Não é página de Confluence — a versão canônica, que vai
para o lab, é [`03_step_1_explained.md`](03_step_1_explained.md). O registro formal está em
[`phase_1/results/FINDINGS.md`](../phase_1/results/FINDINGS.md); a decisão que ele pede está na
[página 01](01_sample_frame_decision.md).

---

## Para que serve o Step 1

O pipeline pergunta: *existe algum gene onde variantes raras e danosas aparecem mais em quem tem
perda auditiva do que em quem não tem?*

Antes de poder perguntar isso, duas coisas precisam estar resolvidas — **quem tem perda auditiva** e
**quem está no estudo**. O Step 1 resolve as duas. Tudo que vem depois herda essa resposta, e é por
isso que ele roda primeiro e que um erro aqui não tem conserto lá na frente.

---

## O que o Step 1 faz, em quatro movimentos

### 1. Puxar todo diagnóstico relacionado a ouvido do prontuário

**Por quê:** o fenótipo precisa vir de algum lugar, e num biobanco hospitalar esse lugar são os
códigos de cobrança lançados durante o atendimento normal.

Duas armadilhas moram aqui. Códigos de diagnóstico são agrupados em *phecodes* — `SO_396` é perda
auditiva, `SO_39x` é a família auricular inteira — porque um código de cobrança cru é fino demais
para servir de fenótipo. E no PMBB v4 os códigos de zumbido **mudaram de tabela**, separados de todos
os outros diagnósticos. Quem não souber disso perde silenciosamente quase toda a evidência de
zumbido. Este pipeline tratou isso corretamente.

### 2. Decidir quem é caso — a regra de 2

Uma pessoa só conta como caso se o diagnóstico aparecer em **duas datas distintas**.

**Por quê:** uma menção isolada pode ser uma suspeita que foi descartada, um código lançado para
justificar um exame, ou erro de digitação. Dois encontros separados são evidência de que alguém
continuou achando aquilo. É uma regra deliberadamente conservadora: ela custa casos reais — 4.007
pessoas aqui — para não admitir falsos.

### 3. Separar o meio ambíguo

Pessoas com evidência auricular que *não* é perda auditiva ficam de fora — nem caso, nem controle.

**Por quê:** são exatamente as pessoas que borrariam a comparação. Alguém com doença de ouvido pode
ter perda auditiva não diagnosticada, e contá-la como controle saudável faz os dois grupos parecerem
mais parecidos do que são, escondendo sinal verdadeiro. São 9.411 pessoas separadas assim.

### 4. Montar a lista de quem entra na análise

**Por quê:** fenótipo sem dado genético para parear não serve de nada. Quem está no prontuário mas
nunca foi sequenciado não tem como contribuir, e deixá-lo na lista quebraria o software.

**É este movimento que dá errado.**

---

## Onde dá errado

O código são três linhas:

```python
fam_file = ".../Imputed/common_snps_LD_pruned/..._genetic_imputed.commonsnps.ldpruned.ALL.fam"
# Keep only samples with genotype data
matched = pheno[pheno["PMBB_ID"].isin(fam["IID"])].copy()
```

O comentário declara a intenção certa: manter só quem tem dado genético.

O problema é **qual** dado genético. O PMBB guarda dois tipos, vindos de dois processos de
laboratório diferentes:

| | o que é |
|---|---|
| **Sequenciamento de exoma** | lê os genes codificadores direto, letra por letra. Acha variante rara, inclusive nunca vista antes. |
| **Array + imputação** | lê um painel fixo de posições comuns e infere o resto estatisticamente. Barato e abrangente; cego a variação genuinamente rara. |

A maioria dos participantes tem os dois. Alguns têm só um.

**Esta análise é construída sobre exoma do começo ao fim.** O modelo nulo é ajustado em variantes de
exoma (`exome_ldpruned/*`), e cada gene é testado em variantes de exoma (`plink_deduplicated/*`). O
dado imputado não é usado para nada.

Mas a linha acima consulta a lista do **imputado**. Ou seja, a pergunta que ela realmente faz é
*"esta pessoa tem dado de array?"* quando a pergunta necessária era *"esta pessoa tem dado de
exoma?"*.

Uma lista de exoma existia — construída por este mesmo pipeline, na mesma árvore de diretórios. Não
foi ela a consultada.

---

## Quanto custou

```
70.925   pessoas com dado de exoma
57.507   analisáveis: caso ou controle para perda auditiva
57.080   sobreviveram ao filtro
   427   removidas -- 40 casos, 387 controles
```

As 427 fazem parte de um grupo de **517** que têm dado de exoma, têm os componentes de ancestralidade
de exoma completos, e não têm dado de array nenhum. Todas poderiam ter sido analisadas. Nenhuma foi.

Vale ser preciso em duas coisas que costumam confundir na primeira leitura:

- **Todas as 517 têm o mesmo problema** — falta de dado de array. A diferença entre 517 e 427 não é
  sobre qual dado faltava; é que 90 delas já tinham sido excluídas por regra de fenótipo, então
  removê-las de novo não mudou nada.
- **Ninguém foi removido por falta de exoma.** Zero pessoas no release estão sem os componentes de
  ancestralidade de exoma. Se o filtro tivesse consultado a lista de exoma, não removeria ninguém.

---

## Por que isso é inconsistência interna, e não divergência de abordagem

O achado seria mais fraco se fosse questão de preferência — um analista optando pela coorte imputada,
outro pela de exoma. Não é.

**Toda etapa analítica deste pipeline roda sobre exoma.** O filtro barra com base num conjunto de
dados que nenhuma etapa posterior consome. O pipeline discorda de si mesmo, e a checagem que teria
pegado isso — *a lista de participantes bate com o dado que está sendo testado?* — é uma que ninguém
faz, porque parece obviamente verdadeira.

---

## O que **não** está sendo afirmado

- **Não** que o fenótipo está errado. Ele está certo. Casos e controles foram re-derivados de forma
  independente a partir do release e batem pessoa a pessoa, 70.925 de 70.925. A regra de 2, a
  exclusão do meio ambíguo e o tratamento da tabela de zumbido reproduzem exatamente.
- **Não** que houve descuido. A intenção está declarada no comentário e é a intenção correta. Os dois
  arquivos `.fam` diferem por uma palavra num caminho longo, e 0,7% de uma coorte está abaixo do que
  qualquer estatística de resumo revelaria.
- **Não** que os resultados são nulos. 40 casos em 6.752 é perda pequena de poder. O efeito é tornar
  associações mais difíceis de detectar, não criar associações falsas. O que ainda **não se sabe** é
  se quem não tem dado de array difere sistematicamente de quem tem — por época de recrutamento, por
  unidade, por ancestralidade. Se diferir, a perda não é neutra. Isso não foi verificado.

---

## A versão de uma frase

O estudo roda sobre dado de exoma, mas a lista de participantes foi filtrada por quem tinha dado de
array — então 517 pessoas que poderiam ter sido analisadas, 40 delas com perda auditiva, ficaram de
fora por lhes faltar algo que a análise nunca usa.
