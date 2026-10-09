# Fase 5 — resultados e o controle que derrubou o achado · análise

Braço: **ExWAS de variante comum, fenótipo bilateral neurossensorial**.
Rodado em 2026-10-09.

Este documento tem um único assunto central: **por que `BCAM` não é um achado**, e como
eu sei disso.

---

## Onde a fase 4 parou

Restringir o fenótipo de "qualquer perda auditiva" para "bilateral neurossensorial"
remove 3.588 dos 6.752 casos. Mesmo assim, algumas variantes ficaram **mais fortes**:

| coorte | quantas fortaleceram | destaque |
|---|---:|---|
| combined | 2 | `BCAM`, chr19 |
| EUR | 4 | `BCAM` de novo, OR 1,17 → 1,35 |
| AFR | 8 | `chr7:158133844`, OR 1,80 → 2,90 |

A leitura atraente: a associação vive no subtipo bilateral, então afinar o fenótipo
afina o sinal mesmo com menos gente.

A leitura chata: são 70 a 118 mil variantes comparadas por coorte. Remova metade dos
casos ao acaso e **todo** p-valor se mexe. Alguns se mexem para o lado certo por sorte.
Apontar para esses é escolher o vencedor depois da corrida.

**Nada nos números separa as duas.**

---

## Conceito — como se testa isso

Segura o número de casos e tira o fenótipo.

Em vez de remover os 3.588 que não têm bilateral neurossensorial, remove **3.588 ao
acaso**. O sorteio não tem significado nenhum de subtipo. Então o que ele produzir de
"fortalecimento" é o **piso** — a quantidade que a simples remoção de casos compra.

| coorte | casos no amplo | sorteados para | controles |
|---|---:|---:|---:|
| combined | 6.752 | 3.164 | 50.746 |
| EUR | 5.183 | 2.547 | 37.596 |
| AFR | 1.285 | 490 | 10.049 |

**Cinco sorteios, não um.** Um sorteio sozinho pode dar sorte para qualquer lado, e com
um número só eu teria escolha sobre qual resposta reportar.

Isso custou 15 modelos nulos (step 1) e 330 tarefas de step 2. O step 1 **não** pode ser
herdado aqui: o modelo nulo é ajustado ao vetor de fenótipo, e cada sorteio tem casos
diferentes.

---

## A regra foi escrita antes

Está em `phase_5/power_control/README.md`, fixada antes de qualquer número existir:

- **observado ≤ máximo dos sorteios** → é seleção, não subtipo. Não reporta.
- **observado > máximo, e Poisson p < 0,05** → sugestivo de concentração no subtipo.
- **entre os dois** → inconclusivo, com os dois números à mostra.

Escrever a terceira opção importa tanto quanto as outras duas. Sem ela, qualquer
resultado ambíguo viraria o que eu quisesse que fosse.

---

## O resultado

| coorte | restrição bilateral | cinco sorteios | média | Poisson p | veredito |
|---|---:|---|---:|---:|---|
| combined | **2** | 8, 5, 7, 3, 5 | 5,6 | 0,976 | seleção |
| EUR | **4** | 5, 6, 3, 7, 3 | 4,8 | 0,706 | seleção |
| AFR | **8** | 18, 13, 10, 15, 15 | 14,2 | 0,972 | seleção |

A regra dispara na primeira cláusula nas três coortes.

### E não é por pouco — é ao contrário

Em **todas** as coortes o observado está **abaixo da média** dos sorteios. No combined e
no AFR está abaixo ou igual a cada sorteio individual.

Ou seja: tirar 3.588 casos **ao acaso** produz *mais* variantes "mais fortes" do que
tirar os 3.588 que não têm bilateral neurossensorial.

A restrição não está só falhando em concentrar sinal. Ela é levemente **menos** produtiva
de fortalecimento aparente do que a sorte pura.

**`BCAM` não é nada.** E eu estava a um passo de escrever que era.

---

## Por que os sorteios produzem tantos

O AFR mostra o mecanismo mais claro: média de 14,2, a maior das três. Com 490 casos a
estimativa de efeito é a mais ruidosa que temos, então mais variantes passeiam mais longe
quando os casos mudam.

É o mesmo fato que aparece na curva de poder, visto de outro ângulo.

---

## O custo de poder, em números

Razão de chances detectável a 80% de poder, α = 5e-8:

| | casos | MAF 1% | MAF 5% | MAF 10% | MAF 50% |
|---|---:|---:|---:|---:|---:|
| bilateral combined | 3.164 | 2,27 | 1,45 | 1,31 | 1,18 |
| **amplo** combined | 6.752 | 1,78 | 1,30 | 1,21 | 1,12 |
| bilateral EUR | 2.547 | 2,50 | 1,52 | 1,35 | 1,20 |
| bilateral AFR | **490** | **7,92** | 2,57 | 1,99 | 1,51 |

Efeitos de variante comum em traços complexos ligados à idade ficam em OR 1,05–1,15.

Este braço está abaixo dessa faixa em todo lugar, exceto nas variantes mais frequentes
das duas coortes maiores. **O AFR, com 490 casos, é praticamente não informativo para
variante comum** — precisaria de OR 7,92 em MAF 1%.

---

## O que este resultado NÃO diz

Esta distinção importa e eu preciso saber fazer.

**Não diz** que bilateral neurossensorial é o fenótipo errado. Não diz que a restrição
foi má ideia.

**Diz** que com 3.164 casos a restrição **não consegue demonstrar especificidade de
subtipo em nenhuma direção** — porque o ruído que ela introduz ao cortar metade dos casos
é maior que qualquer afinamento que ela possa produzir. A pergunta fica sem resposta, e
ficaria sem resposta mesmo se a biologia estivesse do nosso lado.

É a mesma conclusão a que o controle do braço de burden chegou, por outro caminho — lá
foram os genes do ClinGen, aqui foram as variantes que fortaleceram.

---

## Auto-exame de qualidade

O mesmo teste que aplicamos aos achados do pôster, aplicado a nós:

| coorte | ρ (call rate vs −log₁₀p) | p |
|---|---:|---:|
| combined | −0,0001 | 0,98 |
| EUR | −0,0034 | 0,35 |
| AFR | −0,0031 | 0,29 |

Nenhuma relação. Call rate pior não compra significância aqui.

---

## As figuras

| | o que mostra |
|---|---|
| `qq.png` | calibração, bilateral contra amplo |
| `manhattan.png` | o exoma inteiro, com as duas barras |
| `crossarm.png` | o efeito da restrição, com as que fortaleceram circuladas |
| `power.png` | o custo de cortar os casos pela metade |
| `power_control.png` | **o observado contra os cinco sorteios** |

A última é a que carrega a conclusão. A linha vermelha (observado) fica em cima ou
abaixo de todas as barras azuis (sorteios).

---

## Um detalhe operacional que vale guardar

Uma das 330 tarefas saiu sem erro nenhum: `.err` vazio, sem mensagem, SAIGE nunca
começou. Mesma falha transitória de nó que já tinha acontecido no outro braço.
Resubmetida sozinha, rodou normal.

E eu só notei porque o `04_readout.py` **recusa** emitir veredito com sorteio incompleto.
Sem essa recusa, o controle teria rodado com 4,97 sorteios e eu não saberia.

Antes disso, no step 1, aconteceu a outra versão do mesmo problema: contei os `.rda` e
deu 15 de 15 — **nove deles tinham 0 bytes**, porque o SAIGE cria o arquivo quando a
tarefa começa, não quando termina. Contar arquivo não é contar trabalho feito.

---

## Em aberto

- **Nada deste braço precisa de ação.** O achado foi testado e caiu.
- **A P9 nunca foi acionada** — não houve o que anotar.
- **A pergunta aberta com a Nikki** é do outro braço, não deste.
