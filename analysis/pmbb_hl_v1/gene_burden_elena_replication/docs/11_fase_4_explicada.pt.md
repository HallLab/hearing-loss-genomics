# Fase 4, explicada — onde a pergunta finalmente é feita

**Autor:** Andre Rico · **Data:** 2026-10-05
**Cópia de trabalho pessoal, em português.** Não é página de Confluence.
**Resultados e achados:** [`phase_4/results/FINDINGS.md`](../phase_4/results/FINDINGS.md)
**Desenho e parâmetros:** [`phase_4/PLAN.md`](../phase_4/PLAN.md)

---

## Onde estávamos

As três primeiras fases montaram peças. Nenhuma delas respondeu nada.

| fase | montou | a pergunta que responde |
|---|---|---|
| 1 | quem está no estudo | quem tem perda auditiva e quem não tem? |
| 2 | quais variantes contam | quais mutações têm chance de estragar um gene? |
| 3 | o que descontar | o que mais, além do gene, explica quem fica surdo? |

É como ter separado os participantes, a lista de mutações relevantes e a ficha de cada pessoa — e
ainda não ter olhado se alguma coisa se relaciona com alguma coisa.

## A pergunta da Fase 4

**Gene por gene: as pessoas com perda auditiva carregam mais mutações raras estragadoras nesse gene
do que as pessoas sem?**

É a primeira fase que produz p-valores. A primeira que pode mudar aquilo em que o laboratório
acredita.

---

## Por que somar as mutações por gene

Mutação rara é rara mesmo. Uma variante específica pode aparecer em 3 pessoas das 57.498. Com 3
pessoas não se conclui nada, nunca.

A saída é não testar variante por variante, e sim **somar por gene**. Em vez de "esta mutação está
associada à surdez?", a pergunta vira "**qualquer** mutação estragadora neste gene está associada à
surdez?".

As 3 pessoas de uma variante, mais 5 de outra, mais 2 de outra, viram 40 pessoas com algum dano
naquele gene. E aí dá para testar.

Isso é o **burden**: a carga de dano que a pessoa tem naquele gene.

É também por isso que a Fase 2 importava tanto. Se a lista de "mutações estragadoras" está errada,
você soma as coisas erradas e a conta toda não significa nada.

---

## Por que o SAIGE trabalha em dois passos

O programa que faz a análise se chama SAIGE, e ele divide o trabalho em dois. A divisão é esperta, e
vale entender porque explica quase tudo o que vem depois.

**O que não depende do gene** — quem são essas pessoas, que idade têm, de onde vem a ancestralidade
delas, e sobretudo **quem é parente de quem**. Isso é idêntico quer você esteja testando o primeiro
gene ou o décimo oitavo milésimo.

**O que depende do gene** — as mutações raras naquele gene específico.

### Passo 1: a linha de base

Constrói a expectativa de quem tem perda auditiva **sem olhar nenhum gene**. Só idade, sexo, lote de
sequenciamento, ancestralidade e parentesco. Roda uma vez por coorte: três vezes no total.

O parentesco é a parte caríssima. Com 57.498 pessoas, para saber quem é parente de quem você compara
**cada par** de participantes: 1,65 bilhão de pares. E isso não tem nada a ver com o gene que você
está testando.

Na prática, esse passo levou **15 minutos**. Se fosse refeito para cada gene, seriam uns **6 meses**
de processamento por coorte — perto de um ano e meio para as três — para chegar na mesma resposta.

### Passo 2: a pergunta por gene

Para cada gene, verifica se as mutações raras explicam algo **além** daquela expectativa.

Uma analogia: imagine avaliar se uma escola ensina bem. Primeiro você monta a previsão de nota de
cada aluno a partir do que não tem nada a ver com a escola — renda da família, notas anteriores.
Depois você pergunta, escola por escola, se os alunos dela vão melhor do que a previsão dizia.

O passo 1 é a previsão. O passo 2 é a pergunta sobre cada escola. Montar a previsão uma vez e reusá-la
em todas é o que torna a conta viável.

### E tem uma razão que não é só economia

Mesmo se fosse barato, separar seria o certo.

A linha de base precisa ser construída **sem nunca ter visto o gene**. Se você ajustasse idade, sexo,
ancestralidade, parentesco e o gene todos ao mesmo tempo, eles competiriam pela mesma variação — e o
modelo poderia absorver parte do sinal real do gene dentro das covariáveis, ou o contrário. Você não
saberia dizer qual.

Separando, a linha de base fica congelada antes de qualquer gene ser perguntado. Resultado: **todos
os genes são medidos com a mesma régua**.

### Uma consequência prática

Como os 594 trabalhos do passo 2 se apoiam nos três modelos nulos do passo 1, **um erro no passo 1
contamina tudo**. Não existe um gene que escape.

Foi por isso que comparei os números do passo 1 com os da Elena antes de submeter o passo 2, em vez
de só verificar que os trabalhos terminaram sem erro. Terminar sem erro e estar certo são coisas
diferentes.

---

## Qual é o teste estatístico

Não é razão de verossimilhança (LRT). O SAIGE usa um **teste de score**, e a diferença importa.

| teste | o que faz | quantos ajustes de modelo |
|---|---|---|
| LRT | ajusta sem o gene e com o gene, compara | 2 |
| Wald | ajusta com o gene, olha o efeito sobre o erro | 1, com o gene |
| **score** | ajusta **só sem** o gene, e pergunta: "para que lado esse modelo gostaria de se mover?" | 1, sem o gene |

O teste de score só precisa do modelo do passo 1. É literalmente o que torna a estrutura de dois
passos possível.

E tem uma segunda vantagem, menos óbvia. Com dado raro o teste de score continua **válido** onde os
outros quebram. Teve um gene nos nossos resultados, o `CFH`, com 5 alelos no total, todos em casos e
nenhum em controles. Isso é separação perfeita, e o LRT e o Wald precisam estimar o efeito nessa
situação — a estimativa vai para infinito. O teste de score nunca precisa estimar nada nessa
situação, então não tem o problema.

### O "SPA" do nome

O arquivo que roda o passo 2 se chama `step2_SPAtests.R`. SPA é *saddlepoint approximation*,
aproximação de ponto de sela.

O teste de score padrão supõe que a estatística segue uma curva normal. Essa suposição funciona bem
no meio da distribuição e **falha na ponta** justamente quando casos e controles estão
desbalanceados e a variante é muito rara. Nós temos os dois: 6.752 casos contra 50.746 controles, com
variantes de 5 alelos.

Sem corrigir, os p-valores sairiam **otimistas demais** — você produziria significância que não
existe, e precisamente nos genes mais raros, que são os que você está estudando.

É isso que o SAIGE resolve e que uma regressão comum não resolveria.

---

## As três colunas de p-valor

Para cada gene saem três p-valores. Eles respondem perguntas diferentes.

**`Pvalue_Burden` — a soma.** Transforma o gene num número por pessoa: "quantos alelos danosos essa
pessoa carrega aqui?" E testa esse número. É forte **se todas as mutações empurram para o mesmo
lado**. E é frágil do jeito oposto: se metade é danosa e metade é protetora, elas se cancelam na soma
e o teste não vê nada.

**`Pvalue_SKAT` — a dispersão.** Em vez de somar, soma os efeitos **ao quadrado**. Elevar ao quadrado
mata o sinal de menos, então direção deixa de importar: danosas e protetoras contribuem igual.
Sobrevive a direções misturadas. Em troca, perde força quando todas realmente apontam para o mesmo
lado, porque jogou fora a informação de direção.

Sobre a palavra "kernel", que está no nome (*Sequence Kernel Association Test*) e soa opaca: kernel
aqui é só **uma forma de medir o quanto duas pessoas se parecem** naquele gene. O SKAT pergunta:
pessoas que se parecem neste gene também se parecem no fenótipo? É isso.

**`Pvalue` — o SKAT-O, a mistura.** Combina os dois, com o peso escolhido a partir do próprio dado.
É a coluna principal, a que usei em todos os rankings. As outras duas são seus componentes.

### Um exemplo real, e por que ele importa

O gene `PSG8`, na coorte combinada, foi um dos nossos melhores p-valores:

```
Pvalue (SKAT-O) = 0,000004     <- parece achado
Pvalue_Burden   = 0,40         <- a soma não vê nada
Pvalue_SKAT     = 0,000003     <- a dispersão vê tudo
```

O efeito direcional é praticamente zero. A soma das mutações não vai a lugar nenhum. Todo o p pequeno
vem do lado SKAT — ou seja, as mutações desse gene têm efeitos espalhados, sem direção consistente.

Isso é **evidência mais fraca** para uma história de burden do que o 0,000004 sugere olhando sozinho.
Olhar só a coluna principal faria esse gene parecer melhor do que ele é.

### Uma coisa que encontrei checando isso

Em **64% dos nossos 472 mil testes, a soma e a dispersão dão exatamente o mesmo p-valor.** Parece bug,
e não é.

O SAIGE junta as mutações ultra-raras numa só antes de testar, porque mutação com 2 ou 3 alelos não
sustenta estimativa própria. Quando *todas* as mutações do gene são ultra-raras, sobra uma única
unidade — e com uma só unidade, somar e elevar ao quadrado dão o mesmo teste.

É o caso do `CFH`: as quatro mutações foram juntadas numa, e o "gene" na prática é uma única unidade
com 5 alelos. Mais uma razão para não confiar nele.

---

## Por que rodamos 594 vezes e não uma

Porque ninguém sabe de antemão sob qual definição o sinal real vive:

- **3 máscaras** — "só mutações que quebram o gene" (pLOF), "só as que trocam um aminoácido de forma
  provavelmente danosa" (pDM), ou as duas juntas
- **3 cortes de raridade** — mutações em até 1%, 0,1% ou 0,01% das pessoas
- **3 coortes** — todos juntos, só europeus, só africanos
- **22 cromossomos** — divisão só de engenharia, para caber na memória

Testar as combinações de verdade e reportar todas é honesto. Testar várias e mostrar só a melhor é o
que gera resultado que não replica.

### Sobre os cortes de raridade, uma surpresa

Eu esperava que os três cortes fossem três testes diferentes. Em grande parte **não são**.

Em quase um terço dos casos, os três cortes dão p-valor **idêntico** — é o mesmo teste rodado três
vezes. O motivo: em exoma, mutação que quebra gene é quase sempre rarissíma. Afrouxar o corte de
0,01% para 1% acrescenta, em média, **duas mutações** ao conjunto do gene.

Mas nos casos em que os cortes diferem, a diferença pode ser enorme. O gene `MLLT6` tem p de
0,000008 no corte estrito e 0,053 no frouxo — seis mil vezes de diferença. Esse é o padrão de sinal
concentrado nas mutações mais raras e **diluído** quando entram as um pouco mais comuns. Faz sentido:
mutação que realmente destrói função é mantida rara pela seleção natural.

---

## O que saiu

**Nenhum gene atingiu significância, em nenhuma coorte.** O melhor p que temos é 4,2 × 10⁻⁶, contra
um limiar de 2,8 × 10⁻⁶ no cálculo mais generoso que se consegue defender.

**E o braço da Elena também não.** O mais perto que qualquer um dos dois chega é um gene dela que
errou o limiar por 1,5%.

### Por que "nada" é um resultado, e não uma falha

"Nenhum gene passou a barra" e "não tem nada aí" são afirmações diferentes. Uma análise mal calibrada
produz a primeira sem sustentar a segunda. Então verifiquei duas coisas.

**A análise está calibrada.** Se tudo fosse ruído, a fração de testes abaixo de um p-valor deveria ser
exatamente aquele p-valor. Medimos entre 0,86 e 1,15 do esperado — nem inflada (que invalidaria tudo)
nem deflacionada (que esconderia sinal real).

**Nada replica entre ancestralidades.** Europeus e africanos são amostras independentes. Se houvesse
efeito real, genes no topo de uma deveriam estar no topo da outra mais do que o azar permite. Não
estão: 22 observados contra 18,8 esperados.

Então a ausência é de verdade. **Nesta coorte, com este tamanho, não há achado de burden de variante
rara para perda auditiva.** 6.752 casos é pouco para esse tipo de análise.

---

## Onde a replicação valeu, então

Não ficou em "achamos um gene que ela perdeu". Ficou em duas coisas.

### O gene número 1 dela não existe

O `TMC3-AS1` é o gene em **primeiro lugar** nos resultados dela, em duas das três coortes. Ele é um
`lncRNA` — um gene que não produz proteína.

A máscara pLOF pergunta se **perder a função da proteína** se associa à doença. Esse gene não tem
proteína para perder. O resultado todo repousa em **duas** mutações.

Nas nossas máscaras ele não existe: caiu na Fase 2 junto com os 1.101 genes não-codificantes que a
Nikki levantou.

E o nome importa para o quanto isso é perigoso. `TMC3-AS1` é o antisense do `TMC3`, e o `TMC1` é gene
de surdez bem estabelecido. Quem passa o olho numa tabela de genes vê "TMC" e lê plausibilidade. Um
artefato que **parece** achado é mais perigoso que um que parece ruído.

### Metade da lista de topo muda

Dos 50 genes do topo de cada coorte, compartilhamos entre 15 e 24 com ela. Metade a dois terços de
cada lista é diferente.

Numa coorte maior — ou num estudo que escolhe genes para investigar a partir dessa tabela, que é para
isso que tabela de top genes serve — **os defeitos decidiriam o que seria investigado.**

---

## Uma coisa que anotei, com ressalva

O `GJB3` (conexina 31) é gene de surdez estabelecido, e apareceu no **rank 14** de 17.943 na coorte
combinada. E com o comportamento certo: no corte estrito de raridade o sinal está lá, e ao afrouxar
ele desaba enquanto nove mutações mais comuns entram.

**Isso é uma anotação, não um achado.** O p-valor dele é duas ordens de magnitude longe do limiar. E
com 17.943 genes testados, um gene de surdez conhecido cair no rank 14 por azar não é nada
surpreendente — é o tipo de coincidência que acontece todo dia.

O que o torna digno de nota é só que identidade e comportamento apontam para o mesmo lado. Fica como
item para a Fase 5.

### E dois genes conhecidos que andaram em direções opostas

O `SIX1`, também gene de surdez, **melhorou** bastante no nosso braço: rank 38 para rank 7. O `COCH`
**piorou**: rank 12 para rank 22.

Dois genes estabelecidos indo para lados opostos **não** é evidência de que o braço corrigido
recupera biologia conhecida melhor. É só mais evidência de que as listas diferem, o que já sabíamos.
Registro assim porque seria fácil contar só a metade que favorece a correção.

---

## Resumo em cinco linhas

1. A Fase 4 é onde a pergunta finalmente é feita: gene por gene, quem tem perda auditiva carrega mais
   mutações estragadoras?
2. O SAIGE separa em dois passos porque a parte caríssima — parentesco — não depende do gene.
3. O teste é de score com correção de ponto de sela, e saem três p-valores por gene: soma, dispersão,
   e a mistura dos dois.
4. **Resultado: nada significativo**, e verifiquei que a ausência é real, não análise quebrada.
5. O valor da replicação foi metodológico: o gene nº 1 dela é um gene sem proteína, e metade da lista
   de topo muda quando os defeitos são corrigidos.
