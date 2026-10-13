# Guia do CONSOLIDATED — em português

Para eu ler **antes** da reunião. O `CONSOLIDATED.md` é o documento com os números; este
aqui é o que eu preciso ter na cabeça para apresentá-los e para aguentar as perguntas.

---

## A mensagem, em uma frase

> Rodamos seis análises. Nenhuma achou nada significativo — mas sabemos **por quê**, e o
> pipeline está funcionando. O que falta é tamanho de coorte, não método e não fenótipo.

Se eu só conseguir dizer uma coisa, é essa.

---

## Os seis blocos

Dois métodos × três rodadas:

|  | Elena (original) | Replicação (nossa) | Bilateral (nossa) |
|---|---|---|---|
| **Variante rara** (burden) | entregue | refizemos e auditamos | fenótipo novo |
| **Variante comum** (ExWAS) | entregue | refizemos e auditamos | fenótipo novo |

A Elena nunca rodou bilateral. Aquela coluna é trabalho novo, não replicação — **não
posso chamar de "replicação do bilateral"**.

---

## Os números para saber de cor

| | |
|---|---|
| **Nada significativo** | em nenhuma das 6 análises, nenhuma coorte |
| **O mais perto que chegamos** | 1,9× aquém da barra (burden, Elena, combined) |
| **O mais longe** | 616× (ExWAS, bilateral, combined) |
| **Variante rara chega ~10× mais perto** que variante comum | burden 1,9–25× · ExWAS 21–616× |
| **A replicação do ExWAS bate com a entrega dela** | r = 0,98–0,99 · 99,9% das variantes em comum |
| **O único resultado positivo** | genes de surdez do ClinGen no topo do burden: 5 no top-50 contra 0,26 esperado, **p = 5,9e-06** |
| **Coorte** | 6.752 casos (amplo) · 3.164 (bilateral) · ~50.746 controles nos dois |

---

## Os três conceitos que eu vou ter que explicar

### 1. Por que um resultado nulo vale apresentação

Porque o nulo vem com uma **medida de quanto poderíamos ter visto**.

Com 6.752 casos, só detectaríamos razão de chances a partir de ~1,2 nas variantes comuns
mais frequentes, e ~1,8 nas de 1%. Efeitos de variante comum em traços complexos ligados
à idade ficam em **1,05 a 1,15**.

Então a frase certa é: *"nenhum efeito maior do que esta coorte conseguiria enxergar"* —
e **não** *"não existe associação"*. A diferença entre as duas frases é o documento
inteiro.

### 2. O enriquecimento ClinGen

Pegamos a lista de genes de surdez já estabelecidos (ClinGen, classificação Definitive ou
Strong) e perguntamos: eles aparecem no topo da nossa lista mais do que o acaso explica?

Sim. Cinco deles nos 50 primeiros, quando o esperado era 0,26. `p = 5,9e-06`.

**O que isso significa:** o pipeline encontra biologia real. O sinal existe, mas está
**abaixo do piso de detecção**.

**O que isso NÃO significa:** nenhum daqueles cinco genes é significativo sozinho. Não
posso falar "achamos COCH" — não achamos.

### 3. Por que as duas barras

- **5e-8** — convenção, e a que o pôster usou
- **0,05 ÷ número de testes** — Bonferroni sobre o que realmente rodamos

Fixamos as duas **antes** de rodar. Isso importa porque, se uma variante caísse entre as
duas, a escolha de qual usar decidiria o resultado — e escolher depois de ver o número
não é escolher, é torcer.

Nenhuma chegou perto de nenhuma das duas.

---

## PONTOS DE ATENÇÃO

### A. Onde nós estamos frágeis

**1. O MHC.** Nossos PCs (4 e 3) deixam inflação naquela região que os dela (9 e 10) não
deixam — λ de 1,34 no EUR contra 0,97 dela.

- *Se perguntarem:* é real e eu reporto por iniciativa própria, não por ser pego. Fora do
  MHC as duas rodadas são indistinguíveis (λ 0,99 contra 0,99). E λ de 1,34 em 1.032
  variantes não fabrica um resultado em 5e-8 — o menor p do EUR é 2,12e-05.
- *O que não fazer:* minimizar. É uma limitação nossa e está no documento com o mesmo
  peso dos defeitos dela.

**2. Os braços de ExWAS herdam coorte, covariáveis e modelo nulo dos de burden.**
Se aquelas pastas mudarem, os ExWAS ficam desatualizados e nada avisa.

- *Está registrado* nas premissas P1–P3 de cada braço.

**3. Nossos PCs diferem dos dela** (4/3 contra 9/10). É desvio documentado e defendido
por scree plot, mas é desvio — não dá para dizer "rodamos igual".

**4. O enriquecimento ClinGen não distingue os braços entre si.** Fisher exato entre cada
par: p = 0,436 / 0,056 / 0,678 / 0,495.

- *Ou seja:* **não posso dizer** "nossas correções encontraram a biologia que ela
  perdeu". 5 contra 2 parece muito, mas não é estatisticamente distinguível. E o teste
  já é otimista, porque os braços compartilham pessoas e genes.

### B. O que eu NÃO posso dizer

| Não dizer | Por quê | Dizer no lugar |
|---|---|---|
| "Nosso resultado difere do dela" | concordam a r = 0,99 | "reproduzimos a entrega dela quase exatamente" |
| "Nossas correções acharam a biologia" | Fisher não sustenta | "o enriquecimento aparece nos dois, e a diferença não é distinguível" |
| "A restrição bilateral destruiu o sinal" | o controle diz poder perdido | "com 3.164 casos a pergunta fica sem resposta" |
| "O pôster está errado" | é pergunta de procedência | "aquelas variantes não passam o QC que o próprio pôster descreve" |
| "Não há associação com perda auditiva" | nunca testamos isso | "nenhum efeito do tamanho que esta coorte veria" |
| "Achamos COCH / MYO6 / SIX1" | nenhum é significativo | "genes conhecidos se agrupam no topo mais que o acaso" |

### C. A questão do pôster — a mais delicada

**Os fatos:** o pôster reporta seis resultados significativos. As três variantes por trás
deles existem no release do PMBB, mas **as três falham o filtro de call rate que o próprio
pôster declara**. A do `KNL1`, que carrega o `3,51e-12`, tem genótipo em **20% das
pessoas**.

**Como eu conduzo:**

1. **Começo pela concordância.** Os arquivos que a Nikki entregou batem com os nossos a
   r = 0,99, e a coluna `Significance` dela já diz "Not significant". Quem discorda do
   pôster é o pipeline dela, não o nosso.
2. **O número que eu destaco é medido no release**, com plink2 — não depende de nenhum
   dos dois pipelines. 20% de chamada. Não há o que replicar ou contestar nisso.
3. **A cronologia é inocente.** O pôster saiu no verão; o QC foi apertado depois; os
   arquivos entregues são de 2 de outubro. É sequência normal de correção.
4. **Eu sou coautor.** É o meu nome também. Isso me dá a legitimidade e a obrigação — e
   muda o tom de "reportar o trabalho de alguém" para "corrigir o nosso".

**O que eu não faço:** transformar isso no assunto principal da reunião. É um item, não a
pauta.

---

## As perguntas que podem vir

**"Se nada é significativo, o que vocês fizeram?"**
Mostramos que o pipeline funciona (ClinGen, p = 5,9e-06, e r = 0,99 na replicação),
medimos o que faltaria para achar algo, encontramos seis defeitos no pipeline original, e
testamos duas hipóteses que pareciam promissoras e caíram. Um nulo medido vale mais que um
positivo não verificado.

**"Por que vocês estão certos e ela errada?"**
Não é isso. Concordamos com os resultados dela a r = 0,99. O que não fecha é o pôster
contra os arquivos do pipeline que ele descreve.

**"Por que rodaram o ExWAS se o Doug disse que não valia?"**
Porque era o único braço com **alvo declarado** — `KNL1` e `ZFPM1` do pôster. Todos os
outros não tinham hipótese prévia para confirmar ou refutar. E a ressalva biológica dele
continua de pé: está escrita nas premissas dos dois braços.

**"Esses 5 genes do ClinGen não são um achado?"**
São um achado sobre o **pipeline**, não sobre os genes. Nenhum dos cinco passa a barra
sozinho. O que o enriquecimento diz é que o sinal existe e está abaixo do piso.

**"Por que os PCs diferentes?"**
Escolha por scree plot, documentada. Tem um custo que eu mesmo reporto: no MHC a escolha
dela é melhor que a nossa.

**"O bilateral não deu certo então?"**
Deu certo como análise — o que não deu foi resposta. Com 3.164 casos a restrição não
consegue demonstrar especificidade de subtipo **em nenhuma direção**. Não é que bilateral
seja fenótipo errado; é que a pergunta não é respondível nesse tamanho.

---

## Os dois controles — a parte de que eu tenho mais orgulho

Em ambos, o controle **destruiu um resultado que já estava escrito**.

| | o que parecia | o que o controle disse |
|---|---|---|
| **Burden** | o bilateral perdeu o enriquecimento ClinGen → a restrição matou o sinal | cinco sorteios aleatórios do mesmo tamanho também perdem. **Poder, não sinal.** |
| **ExWAS** | `BCAM` fortaleceu em duas coortes ao restringir o fenótipo → especificidade de subtipo | sorteios aleatórios produzem **mais** fortalecimento que a restrição real. **Seleção.** |

As duas regras de decisão — incluindo o que contaria como inconclusivo — foram escritas
**antes** dos números existirem.

**Por que isso vale apresentar:** é a diferença entre um resultado e uma coincidência
bem-vestida. Sem o controle do ExWAS, `BCAM` estaria no slide.

---

## Como usar os dois documentos

| | |
|---|---|
| `CONSOLIDATED.md` | os números. Levo aberto, projeto as tabelas |
| `GUIA_PT.md` (este) | leio antes. Não projeto |
| `figures/how_close.png` | o slide de abertura |

Se alguém pedir o detalhe de um braço, cada pasta tem `confluence/` com o fechamento em
inglês e `andre_notes/` com a versão em português.

---

## Uma coisa que eu levo da experiência inteira

Três números meus estavam errados só na montagem deste documento — a sobreposição top-50
(escrevi de memória), a comparação do burden na figura (comparei quantidades diferentes),
e quase citei 699.784 linhas como defeito quando só 6.345 eram.

Nenhum dos três quebrou nada. Todos passariam despercebidos numa leitura rápida. Por isso
existe o `03_verify.py`, que rededuz 51 números dos arquivos de origem e falha se o
documento não os carrega.

**Número errado num documento de apresentação é número errado dito em voz alta.**
