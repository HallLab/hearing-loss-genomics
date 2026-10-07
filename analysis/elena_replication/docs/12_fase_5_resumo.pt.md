# Fase 5 — lendo o resultado · resumo

**Autor:** Andre Rico · **Data:** 2026-10-07
**Cópia de trabalho pessoal, em português.** Não é página de Confluence.
Detalhamento completo: [`phase_5/results/FINDINGS.md`](../phase_5/results/FINDINGS.md) ·
Notebook com os gráficos: [`phase_5/notebooks/01_results.pt.ipynb`](../phase_5/notebooks/01_results.pt.ipynb) ·
A Fase 4 explicada do zero: [`11_fase_4_explicada.pt.md`](11_fase_4_explicada.pt.md)

---

## Por que esta fase existe

A Fase 4 produziu **472.453 números** — um p-valor para cada combinação de gene, máscara, corte de
raridade e coorte. Ninguém consegue olhar para isso.

A Fase 5 transforma esses números em figura e tabela. E faz uma coisa que a Fase 4 não conseguia
fazer sozinha: **verificar se o "não encontramos nada" é de verdade.**

Essa distinção é o coração da fase. "Nenhum gene passou a barra" e "não tem nada aí" parecem a mesma
frase, mas não são. Um teste quebrado ou excessivamente cauteloso produz a primeira sem sustentar a
segunda — ele não acha nada porque não **consegue** achar nada, não porque não há nada. Precisávamos
saber em qual dos dois casos estávamos.

## O que fizemos

Três conjuntos de figuras, quatro tabelas, e um caderno de código em duas línguas rodando o mesmo
código — a versão em inglês é a oficial para o laboratório, a em português é minha.

---

## Resultado 1 — o gráfico de Manhattan: nada cruza a linha

### O que é um Manhattan

Imagine os 22 cromossomos esticados lado a lado no eixo horizontal, formando o genoma inteiro da
esquerda para a direita. Cada gene vira um ponto. **Quanto mais alto o ponto, mais forte a
associação com perda auditiva.**

O nome vem da aparência: quando há achado de verdade, surgem torres de pontos subindo acima do resto,
e o gráfico parece o skyline de Manhattan.

Desenhamos uma **linha vermelha** na altura do limiar de significância — o ponto a partir do qual o
resultado conta. Com quase 18 mil genes testados, muitos parecem interessantes por puro azar, e a
linha é onde se separa o que resiste a essa correção.

### O que apareceu

**Nenhum ponto cruza a linha vermelha. Em nenhum dos nove painéis. Em nenhuma das três coortes.**

27 gráficos, e a linha vermelha está sozinha lá em cima em todos. Não há skyline — há uma cidade
plana.

### Qual correção a linha usa, exatamente

A linha é **Bonferroni**, não FDR: `0,05 dividido pelo número de genes daquele painel`. Ou seja, são
nove correções separadas por coorte — cada painel tratado como se fosse a análise inteira.

Isso é de propósito, e vale entender por quê. Essa é a barra **mais frouxa** das três que dá para
defender:

| barra | divide 0,05 por | combined |
|---|---|---|
| **a linha do gráfico** | genes daquele painel (~17.800) | 2,81 × 10⁻⁶ |
| por gene da coorte | 17.943 genes | 2,79 × 10⁻⁶ |
| por teste | todos os 159.960 testes | 3,13 × 10⁻⁷ |

A linha do gráfico não cobra por termos olhado nove painéis. Se você varre nove painéis e pega o
melhor, você fez nove vezes mais buscas do que um painel só — então a barra honesta para "passou
alguma coisa nesta coorte?" é uma das duas de baixo.

Desenhei a mais frouxa justamente porque **nada cruza nem ela**. Assim a conclusão não depende de
qual correção alguém prefere defender. Se algo tivesse cruzado, essa distinção precisaria ser feita
com todas as letras.

---

## Resultado 2 — e o FDR também não acha nada

### Bonferroni contra FDR, em uma linha cada

São duas formas de corrigir o fato de você estar testando 18 mil genes de uma vez.

- **Bonferroni** é o rigoroso: divide o limiar pelo número de testes. Protege muito contra falso
  positivo, e por isso perde achado fraco porém real.
- **FDR** (taxa de falsa descoberta) é o tolerante: em vez de evitar *qualquer* falso positivo,
  aceita uma proporção deles — tipicamente 5% da lista final. Acha mais coisa.

O FDR é **onde um sinal fraco mas verdadeiro apareceria primeiro.** Se o Bonferroni não acha nada
mas o FDR acha, isso é informação de verdade.

### O que apareceu

Rodamos o FDR nos dois níveis — dentro de cada painel e sobre a coorte agregada:

| coorte | menor q, no painel | menor q, na coorte | genes com q < 0,05 |
|---|---|---|---|
| combined | 0,075 | 0,374 | **0** |
| EUR | 0,271 | 0,785 | **0** |
| AFR | 0,117 | 0,260 | **0** |

Zero em todas, nos dois níveis.

E os **próprios q-valores** dizem mais que a contagem. O melhor gene de cada coorte carrega entre
**26% e 79% de chance de ser falso positivo** quando se considera a coorte inteira. Isso não é
"passou perto".

Isso fecha a porta de um jeito que o Bonferroni sozinho não fecha. Um estudo com efeito real, mesmo
fraco, normalmente mostra *alguma coisa* sob FDR. Aqui não mostra nada, e não é por pouco.

*(Não desenhei linha de FDR nos gráficos: o corte do FDR é o maior p com q < 0,05, e não existe
nenhum — a linha não teria onde ficar. Cada painel traz o menor q anotado no canto, que diz a mesma
coisa sem inventar um limiar que não existe.)*

---

## Resultado 3 — o gráfico QQ: a ausência é real

### O que é um QQ

Este é o gráfico que responde "o teste está funcionando?", e é por isso que ele existe.

A ideia: se **nada** no estudo tivesse efeito nenhum, os p-valores ainda assim teriam uma
distribuição previsível, conhecida de antemão. O QQ compara o que **observamos** com o que seria
**esperado** nesse mundo sem efeito algum.

- pontos **em cima da diagonal** → o teste se comporta como deveria
- pontos **acima da diagonal** → inflação; os p-valores estão otimistas demais e não servem
- pontos **abaixo da diagonal** → o teste está conservador demais e esconderia achado real

Junto vai um número, o **λ** (lambda), que resume o gráfico: 1,0 é perfeito.

### O que apareceu

Nos 27 painéis, o λ vai de **0,927 a 1,108**, com mediana **1,002**. Nenhum fora da faixa saudável.
Os pontos ficam colados na diagonal, dentro da faixa cinza de tolerância, e levemente abaixo dela no
meio.

**Tradução:** o teste está calibrado e um pouquinho conservador. Não está inflado, o que invalidaria
tudo. Não está deflacionado a ponto de esconder sinal.

Isso é o que nos autoriza a dizer **"não tem nada aí"** em vez de só "nada passou". E bate com uma
verificação que a Fase 4 tinha feito por um caminho completamente diferente — dois métodos
independentes, mesma resposta.

---

## Resultado 4 — a maior parte do topo é frágil

Mesmo sem ninguém significativo, existe uma lista dos "mais próximos". Era de esperar que alguém
pegasse essa lista e começasse a investigar os primeiros.

Então marcamos cada um com um sinal de **frágil**, por três motivos possíveis:

- o sinal vem todo da dispersão e nenhum da soma (sem direção consistente)
- o gene inteiro colapsou numa unidade só, porque todas as mutações eram ultra-raras
- o gene tem menos de 20 alelos no estudo inteiro

| coorte | marcados, dos 30 primeiros |
|---|---:|
| combined | 16 |
| EUR | 18 |
| **AFR** | **26** |

O AFR é o caso extremo, e a razão é o tamanho: 11.334 pessoas, 1.285 casos. Nessa escala, os menores
p-valores são alcançados por genes com um punhado de alelos. O gene que lidera aquela coorte, o
`CFH`, repousa em **cinco alelos** — tire uma pessoa do estudo e ele evapora.

**A lição prática:** a fração marcada deve ser lida **antes** dos nomes dos genes. Uma tabela de
top-30 do AFR não é uma lista de 30 candidatos.

---

## Resultado 5 — a figura pela qual tudo isso existiu

Esta é a única figura que a Elena não tem equivalente, e é a que resume a replicação inteira.

Cada ponto é um gene. O p-valor **dela** no eixo horizontal, o **nosso** no vertical. Genes na
diagonal são os que as correções deixaram em paz. Quanto mais longe da diagonal, mais a correção
mexeu naquele gene.

E embaixo, uma **faixa vermelha**: os genes que ela testou e que nossa Fase 2 removeu por não
produzirem proteína. Eles não têm posição vertical, porque nunca os testamos — são
**1.099 / 1.011 / 897** genes, conforme a coorte.

Deixá-los fora da figura teria escondido o achado mais importante. Porque **a marca mais à direita
dessa faixa, no painel principal, é o `TMC3-AS1` — o gene em primeiro lugar nos resultados dela.**

Ele aparece como o ponto mais extremo do eixo dela e simplesmente não existe no nosso.

### E o resto da lista

Dos 50 genes do topo de cada coorte, compartilhamos entre **15 e 24** com ela. Metade a dois terços
de cada lista é diferente.

---

## O que concluímos

**1. Nenhum achado.** Nada atinge significância, em nenhuma coorte, em nenhum dos dois braços — nem
o nosso nem o dela. Nem por Bonferroni, nem por FDR, em nenhum nível de agregação.

**2. A ausência é real.** O teste está calibrado. Não encontramos nada porque não há nada a
encontrar nesta coorte, com este tamanho. 6.752 casos é pouco para análise de variante rara.

**3. As correções ainda importaram muito.** Metade de cada lista de topo muda, e o gene nº 1 dela é
um que nunca testamos porque ele não poderia estar ali.

### Por que o ponto 3 não é consolo

Pode parecer que, se ninguém é significativo, errar o ranking não importa. Importa, por dois motivos.

Primeiro: tabela de top genes **serve** para escolher o que investigar depois. Se o laboratório
pegasse aquela lista para desenhar o próximo experimento, o `TMC3-AS1` estaria no topo da fila — um
gene sem proteína, com duas mutações sustentando tudo.

Segundo: numa coorte maior, esses mesmos defeitos estariam operando com os mesmos efeitos, só que
com genes efetivamente cruzando a linha. O tamanho pequeno desta amostra é o que impediu o problema
de virar um achado publicado. Isso é sorte, não método.

**O resultado é negativo, com um achado metodológico sólido dentro.** É diferente de — e mais útil
que — um estudo que não acha nada e não aprende nada.

---

## Em aberto

- **Teste sistemático contra a lista ClinGen de genes de surdez.** Diria se o braço corrigido
  recupera biologia conhecida melhor que o dela. Não fiz porque essa lista vive no `cycle_2`, e a
  regra de isolamento desta pasta mantém o `cycle_2` fora. É endpoint público, dá para buscar de
  forma independente se o time quiser.
- **`GJB3`** (conexina 31, gene de surdez estabelecido) ficou em rank 14 de 17.943 e com o
  comportamento esperado de um gene real. Continua **anotação, não achado** — está duas ordens de
  magnitude longe do limiar.
- **Falta o README** do `elena_replication`, deixado para o fim de propósito. Agora é o fim.
- O braço de reprodução (57.632) **não foi rodado** na Fase 4, por decisão sua. As entradas seguem no
  disco, então ele continua executável se o time quiser o lado a lado completo.
