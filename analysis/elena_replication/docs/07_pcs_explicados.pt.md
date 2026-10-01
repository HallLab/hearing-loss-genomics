# Componentes principais, explicados — e quantos usar

**Autor:** Andre Rico · **Data:** 2026-09-30
**Cópia de trabalho pessoal, em português.** Não é página de Confluence.
**Notebook com os gráficos e o código:** [`phase_3/scripts/01_pc_selection.ipynb`](../phase_3/scripts/01_pc_selection.ipynb)

---

## O problema que os PCs resolvem

Imagine que você quer saber se um gene causa perda auditiva. Você compara portadores com
não-portadores e acha diferença. Parece achado.

Mas suponha que esse gene seja mais comum em pessoas de origem do Leste Asiático, e que essa
população, no seu hospital, tenha sido recrutada numa clínica que atende mais idosos. Aí você achou
uma associação entre o gene e a perda auditiva que é, na verdade, uma associação entre **ancestralidade
e idade**.

Isso se chama **confusão por estrutura populacional**, e é o modo clássico de uma análise genética
produzir achado falso. Grupos populacionais diferem nas frequências de milhões de variantes ao mesmo
tempo, por razões históricas — migração, isolamento, deriva genética. Nada disso tem a ver com
doença.

A correção é **medir a ancestralidade de cada pessoa e descontá-la do modelo.** Os componentes
principais são como se mede.

---

## O que é um componente principal

Suponha que você tenha 70.925 pessoas e meio milhão de variantes genéticas. Isso é uma tabela de
70.925 × 500.000 — impossível de usar como covariável.

Mas a informação ali é muito redundante. Variantes próximas no genoma são herdadas juntas, e
populações inteiras carregam padrões comuns. A **análise de componentes principais (PCA)** comprime
essa redundância: ela procura as **direções de maior variação** na tabela e resume cada pessoa em
poucos números.

Uma analogia: para descrever onde alguém mora no Brasil, você não precisa da lista de todas as ruas
que ele já pisou. Duas coordenadas — latitude e longitude — já colocam a pessoa no mapa. A PCA
encontra as "coordenadas" genéticas equivalentes.

**PC1** é a direção que mais separa as pessoas. **PC2** é a segunda, perpendicular à primeira. E
assim por diante.

### E ninguém diz à PCA o que é ancestralidade

Esta é a parte que impressiona, e o notebook mostra num gráfico. A PCA só vê genótipos — nenhum
rótulo de origem, nenhuma informação demográfica.

Mesmo assim, quando você desenha PC1 contra PC2 e **depois** colore os pontos por ancestralidade, os
grupos aparecem separados. A estrutura emergiu sozinha, porque ancestralidade é, de fato, a maior
fonte de variação genética entre pessoas.

É por isso que os PCs funcionam como medida de ancestralidade: eles não foram feitos para isso, mas é
isso que eles capturam primeiro.

---

## O autovalor, e por que ele importa

Cada PC vem com um **autovalor** — um número que diz **quanta variação aquele componente captura**.

O PC1 tem o maior. O PC2, menor. E a sequência cai. A pergunta prática é: **onde parar?**

```
coorte combinada do PMBB v4:

PC1   762,4   ← captura 33% da variação do top 20
PC2   326,9   ← 14%
PC3   204,9   ← 8,9%
PC4   106,9   ← 4,6%
PC5    81,4   ← 3,5%
PC6    73,7   ← 3,2%
PC7    69,2   ← 3,0%
```

Repare no comportamento: 762 → 327 → 205 → 107 → 81 são quedas grandes. Depois 81 → 74 → 69 é quase
nada.

Do PC6 em diante, **cada componente explica quase o mesmo que o vizinho.** Isso é a assinatura de
ruído: se houvesse estrutura populacional real ali, ela se destacaria.

---

## O scree plot e o "joelho"

O gráfico desses autovalores em ordem se chama **scree plot** — *scree* é o cascalho que se acumula
no pé de uma encosta, e o gráfico tem essa forma: queda íngreme, depois uma base plana.

O ponto onde a curva deixa de cair é o **joelho**. É ali que se corta.

### Como identifico o joelho sem olhar no olho

Olhar o gráfico é subjetivo — depende da escala dos eixos. O critério que usei é a **queda
proporcional** entre PCs consecutivos:

> O joelho é o primeiro PC cuja queda em relação ao anterior fica **abaixo de 10%**. Mantemos até o
> PC anterior a ele.

Aplicado às três coortes do estudo:

| coorte | quedas sucessivas | joelho |
|---|---|---|
| **combined** | 57% · 37% · 48% · 24% · **9,5%** | **5 PCs** |
| **EUR** | 55% · 25% · 12% · **7,1%** | **4 PCs** |
| **AFR** | 44% · 19% · **3,8%** | **3 PCs** |

O limiar de 10% é uma convenção, não uma lei. Mas os três casos aqui não são de fronteira — a queda
cai de 24% para 9,5% no combined, de 12% para 7,1% no EUR, e de 19% para 3,8% no AFR. Mover o limiar
para 8% ou 12% não mudaria nenhuma das três respostas.

---

## Por que não usar PCs de sobra "por segurança"

A intuição diz: se PC é para remover confusão, mais PCs = mais segurança. **É falso**, e por dois
motivos.

**Cada covariável consome graus de liberdade.** O modelo tem um orçamento de informação. Gastar em
covariáveis que não explicam nada reduz o poder de detectar o que você está procurando — você perde
achado real para se proteger de um viés que não existe.

**E PCs do platô podem correlacionar com o fenótipo por acaso.** Um componente que captura ruído
técnico — lote de sequenciamento, data de processamento — pode se correlacionar com qualquer coisa
por coincidência. Incluí-lo pode *introduzir* viés em vez de remover.

Foi exatamente esse o raciocínio da reunião de 2026-07-01 ao rejeitar os 20 PCs do Daniel:

> *"Do not use more than 5–6 PCs for PMBB data; using 20 PCs likely overcorrects for population
> stratification."*

---

## O que achamos no pipeline

| coorte | rodou com | scree sustenta | fonte dos PCs |
|---|---|---|---|
| combined | 5 | **5** ✓ | PCA do release |
| EUR | 9 | **4** | PCA própria, dentro da ancestralidade |
| AFR | 10 | **3** | PCA própria, dentro da ancestralidade |

**O combined estava certo.** EUR e AFR usaram mais que o dobro do que o scree sustenta.

E vale notar: os scree plots **são dela**. Ela os gerou — os arquivos `EUR_variance_explained.tsv` e
`AFR_variance_explained.tsv` estão no diretório dela, com data de 14/jul. O plano de análise dela
dizia *"vou usar um gráfico de variância explicada para determinar quantos PCs — provavelmente 4 ou
5"*.

Ou seja: **ela planejou o método certo, gerou o dado certo, e não aplicou.** Os números usados não
vêm do scree que ela própria produziu.

### Uma assimetria que o notebook expôs

Os PCs do `combined` vêm do release (verifiquei por valor: correlação 0,99965 com `exome_PC1`). Os de
EUR e AFR são PCA que ela rodou **dentro de cada grupo**.

Isso é a prática **correta**: dentro de EUR, a estrutura relevante é fina — gradientes norte-sul da
Europa, por exemplo — e a PCA global não a captura, porque ali o PC1 está ocupado separando
continentes.

E explica por que não existe `variance_explained` para o combined: ela não rodou PCA ali, usou a do
release. Os autovalores do combined vêm do próprio release
(`Exome/PCA/combined/...no1KG.eigenvalues.tsv`).

---

## O que adotamos

```
braço de reprodução :  5 / 9 / 10   ← fixo, é o que rodou
braço corrigido     :  5 / 4 /  3   ← cada número com o scree correspondente atrás
REVEL               :  0,5 primário  +  0,6 sensibilidade
```

O REVEL fica em 0,5 porque é o que a reunião decidiu como primário, com 0,6 como sensibilidade — que
também é o que ela decidiu. Não há motivo para desviar disso.

### E o custo de mudar

Vale saber antes de prometer qualquer coisa numa reunião:

- **mudar REVEL** → refaz máscara (minutos) + SAIGE passo 2 (dias)
- **mudar PCs** → refaz covariáveis (horas) + SAIGE passo 1 **e** 2 (dias)

Não é reversível de graça. Os parâmetros valem ser travados **na** reunião.
