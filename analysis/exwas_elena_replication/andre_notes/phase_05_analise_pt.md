# Fase 5 — resultados, calibração e a auditoria do pôster · análise

Documento em português, linguagem de leigo, para eu conseguir explicar em reunião.
Braço: **ExWAS de variante comum, fenótipo amplo (replicação)**.
Rodado em 2026-10-09.

---

## O que esta fase faz

Três coisas:

1. **Calibração** — o teste está bem comportado, ou está inflado?
2. **Concordância** — nossos números batem com os que a Elena entregou?
3. **Auditoria do pôster** — por que o pôster diz que há seis resultados significativos
   e nenhum dos dois conjuntos de resultados tem algum?

A premissa **P9** (anotar com VEP tudo que passar uma barra) foi registrada antes do
step 2 e **não foi acionada**, porque nada passou. Fica escrito que a etapa existia e
não teve o que anotar — não que ela foi abandonada no caminho.

---

## Conceito 1 — o que é λ (lambda)

O fator de inflação genômica responde a uma pergunta simples: **se não houvesse
associação nenhuma, os p-valores pareceriam com isso?**

Se não há sinal, os p-valores devem ser uniformes entre 0 e 1 — tantos perto de 0,5
quanto perto de 0,9. λ compara o p-valor **mediano** observado com o mediano esperado
nessa situação.

- **λ = 1** → exatamente o esperado sob nenhuma associação
- **λ > 1** → p-valores pequenos demais no geral. Sinal de estrutura populacional não
  corrigida, viés de lote, ou algum problema sistemático
- **λ < 1** → conservador demais

λ **não** detecta achado verdadeiro. Um gene real afeta um punhado de variantes, não
a mediana de 70 mil. λ mede viés de fundo, não descoberta.

### Por que λ é legítimo aqui e precisava de defesa no burden

λ foi inventado para teste de variante única, onde a estatística sob a hipótese nula é
um qui-quadrado com 1 grau de liberdade. **É exatamente o que este braço produz.**

No braço de burden, o p-valor vem do SKAT-O, cuja estatística não é essa. Lá eu
precisei de um parágrafo de justificativa (o argumento é que o que se testa é se os
p-valores são uniformes, e isso vale seja qual for o mecanismo que os gerou). Aqui não
precisa de nada: é o caso de uso original.

Detalhe técnico que mantivemos: uso o `qchisq(0.5, 1) = 0.4549364` exato, e não o
`0.456` arredondado no `step10_1_QQplot.R` dela, que infla λ em 0,23%. Não muda
conclusão nenhuma, mas é de graça acertar.

---

## Conceito 2 — por que reporto duas correlações

Para comparar nossos p-valores com os dela eu uso duas medidas:

| | o que mede |
|---|---|
| **Pearson (r)** | se os valores em si batem |
| **Spearman (ρ)** | se a **ordem** bate — quem é o 1º, o 2º, o 100º |

Por que os dois: quando **nada é significativo**, 70 mil dos 70 mil p-valores são
ruído perto de zero. O Pearson olha essa massa toda e fica alto quase de graça — os
dois lados concordam que "não há nada", e concordar sobre nada é fácil.

O Spearman é mais duro. Ele pergunta se as duas rodadas põem as mesmas variantes nas
mesmas posições do ranking, inclusive no meio da tabela, onde pequenas diferenças
mudam a ordem. Se o ρ cair muito abaixo do r, é sinal de que a concordância era só da
massa de ruído.

Aqui os dois se sustentam, e é por isso que posso dizer que replicamos de verdade.

---

## O resultado: nada significativo, e isso bate com ela

| coorte | menor p nossa | menor p dela | barra primária | barra secundária |
|---|---:|---:|---:|---:|
| combined | 3,02e-06 | 2,57e-06 | 5e-08 | 7,11e-07 |
| EUR | 2,12e-05 | 1,99e-05 | 5e-08 | 6,71e-07 |
| AFR | 4,11e-06 | 3,70e-06 | 5e-08 | 4,24e-07 |

Nenhum dos dois chega a duas ordens de grandeza de `5e-8`. E a coluna `Significance`
dos arquivos dela já diz **"Not significant"** em todos.

### Calibração

| coorte | λ nossa | λ dela |
|---|---:|---:|
| combined | 0,9893 | 0,9919 |
| EUR | 0,9991 | 0,9907 |
| AFR | 1,0207 | 1,0120 |

Os dois bem calibrados. A afirmação do pôster de *"adequate statistical calibration and
minimal evidence of systematic inflation"* **está correta e se reproduz.**

### Concordância

| coorte | dela | nossa | em comum | r | ρ |
|---|---:|---:|---:|---:|---:|
| combined | 70.429 | 70.366 | 70.333 (99,9%) | 0,989 | 0,978 |
| EUR | 74.528 | 74.536 | 74.483 (99,9%) | 0,979 | 0,962 |
| AFR | 117.943 | 117.953 | 117.723 (99,8%) | 0,990 | 0,981 |

**Toda** discordância individual maior que uma unidade log está no MHC — 14 de 14 no
EUR, 16 de 16 no AFR. Não "a maioria": todas.

---

## Um defeito nosso, achado aqui

Separando λ por região em vez de olhar o genoma inteiro:

| coorte | λ fora do MHC (dela / nossa) | λ dentro do MHC (dela / nossa) |
|---|---|---|
| combined | 0,990 / 0,989 | **1,145** / 1,041 |
| EUR | 0,992 / 0,995 | 0,965 / **1,339** |
| AFR | 1,012 / 1,020 | 0,998 / **1,092** |

Fora do MHC somos indistinguíveis dela. Dentro, no EUR e no AFR, **nossos 4 e 3 PCs
deixam inflação que os 9 e 10 dela não deixam.** O EUR é o caso claro: λ = 1,34 contra
0,97 dela, em 1.032 variantes.

### Por que isso acontece

O MHC (cromossomo 6, 28,5 a 33,4 Mb) é a região de maior diversidade do genoma humano.
Ela carrega os genes do sistema imune, e a composição dela varia muito mais entre
ancestralidades do que o resto. É justamente onde a correção por estrutura populacional
pesa mais.

Nossos PCs foram escolhidos por scree plot olhando o genoma inteiro, e essa escolha
continua defensável **para o genoma inteiro** — fora do MHC ela não custa nada. Mas na
janela onde a estrutura é mais forte, os PCs extras dela estão fazendo trabalho real.

### Isso muda alguma conclusão?

Não. λ = 1,34 em 1.032 variantes não fabrica um achado em `5e-8`, e o menor p do EUR é
`2,12e-05`. Mas registro porque é uma limitação **nossa**, encontrada olhando, e porque
qualquer análise futura que se importe com o MHC não deve herdar 4 PCs sem reabrir a
pergunta.

Na coorte combinada a posição se inverte: a inflada é a dela, 1,14 contra nossos 1,04.

---

## A auditoria do pôster

O script `02_poster_audit.py` faz quatro coisas, e a ordem importa:

1. **lê as reivindicações direto do `.pptx`** — não redigitadas por mim
2. **mede cada variante no release do PMBB** com plink2
3. **rastreia** em que estágio do pipeline cada uma desaparece
4. **recusa** escrever o relatório se alguma variante não puder ser medida

Seis resultados reivindicados, três variantes únicas:

| gene | variante | coorte | p no pôster |
|---|---|---|---|
| KNL1 | chr15:40611486:A:G | Combined | 3,51e-12 |
| KNL1 | chr15:40611486:A:G | EUR | 7,21e-12 |
| ZFPM1 | chr16:88533289:A:C | EUR | 2,46e-10 |
| ZFPM1 | chr16:88533290:G:C | EUR | 2,52e-10 |
| ZFPM1 | chr16:88533290:G:C | Combined | 4,54e-09 |
| ZFPM1 | chr16:88533289:A:C | Combined | 6,18e-09 |

Medido no release — **não no nosso pipeline nem no dela**:

| variante | call rate | MAF | passa o QC do próprio pôster? |
|---|---:|---:|---|
| chr15:40611486:A:G | **0,2013** | 0,284 | não |
| chr16:88533289:A:C | 0,9617 | 0,295 | não |
| chr16:88533290:G:C | 0,9473 | 0,283 | não |

O método do pôster diz: *"variants with a genotype call rate < 0.99 and a Minor Allele
Frequency (MAF) < 0.01 were removed."* Nenhuma das três passa de 0,99.

**Para comparar:** o call rate de tudo que este braço testou vai de 0,9534 a 1,0, com
mediana 0,9999. A pior variante do pôster não está na cauda disso — está em outro
regime.

---

## As figuras

Todas em `figures/`, na mesma paleta do braço de burden para poderem conviver num
documento só.

### `qq.png` — calibração

Eixo x: o p-valor que a gente esperaria se não houvesse associação nenhuma. Eixo y: o
observado. Se não há sinal nem viés, os pontos seguem a diagonal.

As duas rodadas estão sobrepostas e quase indistinguíveis — essa é a mensagem. A linha
vermelha pontilhada é `5e-8`, e **nada chega perto dela** em nenhuma das duas.

### `manhattan.png` — o exoma inteiro

Cada ponto é uma variante, na ordem do genoma, com as cores alternando por cromossomo.
Quanto mais alto, menor o p-valor.

As duas linhas vermelhas são as barras: a contínua é `5e-8` (a do pôster), a tracejada
é `0,05/n` (Bonferroni sobre os testes que realmente fizemos). **A faixa vazia entre os
dados e as linhas é o resultado** — o ponto mais alto de qualquer coorte está bem
abaixo das duas.

### `concordance.png` — nós contra ela, variante por variante

Eixo x: o −log10(p) dela. Eixo y: o nosso. Em cima da diagonal = concordância perfeita.

Os pontos vermelhos são o MHC. Repare que **eles são os únicos que se espalham**, e que
se espalham sistematicamente **para cima** — nós achamos mais sinal ali do que ela. É a
mesma coisa que o λ por região mostrou, vista de outro jeito.

### `poster_variants.png` — a figura que carrega o achado

À esquerda, a distribuição de call rate de tudo que foi testado, em escala log, com a
linha de 0,99. À direita, as três variantes do pôster, na mesma escala.

A faixa sombreada à direita marca onde fica tudo que foi testado. As três barras do
pôster estão fora dela, e a do `KNL1` está em 0,201.

---

## Dois erros meus nesta fase

### O grave: uma falha de script virando afirmação falsa

Na primeira rodada, o script de auditoria imprimiu **"absent from the release"** para as
três variantes — e logo abaixo, no mesmo output, a trilha de estágios dizia
**"present"**. Contradição dentro do próprio relatório.

A causa: em Python, `Path("...exome.gl.chr15").with_suffix(".bim")` troca o **último**
sufixo, que ali é `.chr15`, não `.gl`. Ele foi procurar um arquivo
`...exome.gl.bim`, que não existe. O teste de existência falhou, o código caiu num ramo
de erro — e meu `print` renderizava "não consegui medir" como "não existe".

Isso é o pior tipo de bug: ele não quebra, ele **mente com confiança**. E ia mentir
num documento que acusa outro documento de estar errado.

Correções: o caminho é montado por concatenação de string; falha do plink2 agora é
capturada e reportada; e o script **recusa escrever o relatório** se alguma variante
não puder ser medida. Variante não medida nunca mais pode sair impressa como variante
ausente.

### O menor: a figura escondia a cauda

A primeira versão de `poster_variants.png` cortava o eixo em 0,98, o que escondia o fato
de que o AFR tem variantes testadas até 0,9534. Eu teria apresentado "tudo acima de
0,99" quando o correto é "tudo acima de 0,95". O eixo agora vai até 0,95 e os mínimos
de cada coorte estão escritos na legenda.

---

## O memo para a Nikki

Está em `confluence/memo_poster_qc.md`. A ordem dele é deliberada:

1. **Abre pela concordância.** r de 0,98 a 0,99, 99,9% das variantes em comum, λ a 0,01
   de distância, e os arquivos dela já dizendo "Not significant". Assim a pergunta de QC
   chega como problema compartilhado, não como desafio.
2. **O número que destaca é medido no release**, com plink2 — independente dos dois
   pipelines. Call em 20% dos participantes. Não há o que replicar ou contestar nisso.
3. **Inclui a limitação nossa do MHC.** Isso é honesto, e também deixa claro que foi
   auditoria de tudo, não caça.
4. **Fecha com pergunta concreta** — qual conjunto de variantes gerou a tabela — oferece
   a cronologia inocente como resposta provável, registra que eu sou coautor, e dá duas
   opções em vez de exigência.

---

## Em aberto

- **A pergunta para a Nikki.** É o único item que não é nosso para resolver.
- **Os PCs no MHC.** Se alguma análise futura focar aquela região, 4 PCs não serve.
  Não mexi agora porque mudaria o modelo nulo de todos os braços e nenhuma conclusão
  depende disso.
- **A P9 nunca foi acionada.** Se um dia algum braço produzir achado, o passo de
  anotação existe e está pré-registrado.
