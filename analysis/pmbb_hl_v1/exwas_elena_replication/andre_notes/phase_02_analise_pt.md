# Fase 2 — quais variantes contam · análise

Documento em português, linguagem de leigo, para eu conseguir explicar em reunião.
Braço: **ExWAS de variante comum, fenótipo amplo (replicação)**.
Rodado em 2026-10-09, job LSF `49564495`, 66 de 66 tarefas.

Este braço é o de cima do diagrama da Elena:

```
ExWAS                     MAF >= 0,01                   -> SAIGE          <- ESTE
Rare Variant Gene Burden  MAF <= 0,01; 0,001; 0,0001    -> SAIGE-GENE
```

---

## O que esta fase faz

Quase nada, e isso é de propósito.

Nos braços de variante rara, a fase 2 construía as **máscaras**: decidia quais
variantes entram no teste de cada gene, usando anotação funcional. Aqui não existe
máscara. Cada variante é testada sozinha. Então não há o que construir.

O que esta fase faz é **contar**. Uma coisa só: quantas variantes o SAIGE vai
efetivamente testar em cada coorte. Esse número é o denominador de uma das duas
barras de significância, e precisa estar fixado **antes** de existir qualquer
p-valor — senão a barra seria escolhida depois de ver o resultado, o que não é barra.

---

## Conceito 1 — por que contar, se o número está no arquivo?

O arquivo de genótipos tem 13.223.743 variantes. Mas não é esse o número.

**MAF é uma propriedade do grupo de pessoas, não da variante.** Uma variante que
aparece em 0,9% da coorte combinada pode aparecer em 1,2% só entre os de
ancestralidade africana. No primeiro caso ela fica de fora; no segundo, entra. A
mesma variante, o mesmo arquivo, resultado diferente — porque mudou quem está sendo
olhado.

Então o censo tem que ser feito três vezes, uma por coorte, nas pessoas exatas que o
SAIGE vai usar. É por isso que são 66 tarefas: 3 coortes × 22 cromossomos.

---

## Conceito 2 — ALT1 não é MAF

O plink2 reporta a frequência do alelo **ALT** (o alternativo, o que difere da
referência). MAF é a frequência do alelo **menor** — o mais raro dos dois.

Normalmente são a mesma coisa, porque o ALT costuma ser o raro. Mas não precisa ser:
se o ALT estiver em 99,5% da população, o alelo menor é o REF, em 0,5%. A variante
teria MAF de 0,5% e seria descartada, mesmo com ALT "alto".

Então `MAF ≥ 1%` quer dizer **frequência de ALT entre 1% e 99%**, nas duas pontas.

Eu pedi as faixas espelhadas (95%, 99%, 99,9%) de propósito, para **medir** isso em
vez de supor. Resultado: **zero variantes acima de 95%** em qualquer coorte. O REF é
sempre o alelo maior nesses arquivos. A suposição estava certa — mas agora é um
número, não uma torcida.

---

## De onde vêm os dados

| o que | de onde |
|---|---|
| genótipos | `analysis/elena/rarevariant_geneburden/plink_deduplicated/chr{1..22}_deduplicated` |
| quem está na coorte | `../gene_burden_elena_replication/phase_3/results/covariates/` |

Os genótipos são **os mesmos arquivos** que os nossos dois braços de variante rara já
leem, e os mesmos que o ExWAS da própria Elena consumiu. O filtro de call rate > 0,99
foi aplicado quando esses arquivos foram construídos, lá atrás.

Isso tem uma consequência boa: não há como este braço divergir dos outros nos
genótipos, porque não há nada sendo reconstruído.

---

## O resultado

De 13.223.743 variantes no arquivo, passam `MAF ≥ 1%`:

| coorte | amostras | testadas | % do arquivo | barra secundária |
|---|---:|---:|---:|---:|
| combined | 57.498 | 70.383 | 0,53% | 7,10e-07 |
| EUR | 42.779 | 74.553 | 0,56% | 6,71e-07 |
| AFR | 11.334 | 117.996 | 0,89% | 4,24e-07 |

Meio por cento. O exoma é dominado por variantes raríssimas: 12 milhões das 13,2
aparecem em menos de 0,01% das pessoas.

### As duas barras

- **Primária: `5e-8`** — a que a Elena usou no pôster. Convenção de genoma inteiro.
- **Secundária: `0,05 / n_testadas`** — Bonferroni sobre os testes que realmente
  fizemos.

A primária é cerca de **10× mais estrita** que a secundária, nas três coortes. Isso
não é defeito, é o motivo de reportar as duas: uma variante em `p = 2e-7` passaria a
correção sobre os testes reais e falharia no limiar do pôster. Com as duas fixadas
antes, esse caso vira descrição; se eu escolhesse depois, viraria discussão.

---

## Espectro de frequência

| frequência ALT | combined | EUR | AFR |
|---|---:|---:|---:|
| [0; 0,01%) | 12.080.815 | 12.448.698 | 11.869.076 |
| [0,01%; 0,1%) | 826.073 | 572.134 | 944.450 |
| [0,1%; 1%) | 246.472 | 128.358 | 292.221 |
| **[1%; 5%)** | **35.251** | **37.825** | **86.371** |
| **[5%; 95%)** | **35.132** | **36.728** | **31.625** |
| ≥ 95% | 0 | 0 | 0 |

As duas linhas em negrito são o que é testado. Soma dá as colunas da tabela anterior.

---

## Dois pontos de atenção que saíram daqui

### 1. A coorte combinada testa MENOS variantes que a EUR sozinha

70.383 contra 74.553 — com 14.719 pessoas **a mais**.

Parece erro, não é. Juntar ancestralidades empurra para baixo de 1% variantes que são
comuns em um grupo só. Uma variante em 3% entre EUR e 0,1% nos outros grupos pode
fechar em 0,9% na combinada, e sai.

**Por que isso importa:** a combinada **não é um superconjunto** da EUR em termos do
que é testado. Uma variante pode aparecer no resultado da EUR e estar ausente do
resultado da combinada sem que nada esteja errado. Eu ia interpretar isso como bug se
não tivesse medido antes.

### 2. O excesso do AFR está todo colado no limiar

O AFR testa 43.443 variantes mais que a EUR — com um quarto das pessoas. E o excesso
inteiro está na faixa de 1% a 5% (+48.546). Na faixa realmente comum, de 5% a 95%, o
AFR testa **menos** que a EUR (−5.103).

Duas causas, e este censo **não consegue separar** as duas:

- **Diversidade real.** Populações de ancestralidade africana carregam mais variação
  de baixa frequência. Isso aparece também nas faixas mais raras da tabela, o que é
  coerente.
- **Ruído de estimativa.** Com 11.334 pessoas, a MAF é medida com menos precisão do
  que com 42.779. Numa MAF verdadeira de 1%, o erro padrão é de cerca de 0,07 ponto
  percentual no AFR contra 0,03 na EUR. Variantes que estão de verdade logo abaixo
  do corte escapam para cima mais vezes.

**O que vale nas duas hipóteses:** o AFR ganha a barra mais frouxa (4,24e-07)
justamente porque testa mais variantes — mas metade desses testes extras está na
faixa de 1–5% da coorte menor, onde o poder é pior. **A barra é mais permissiva
exatamente onde os testes individuais são mais fracos.**

É um argumento para usar a `5e-8` como manchete e a secundária como contexto. Era
como o pré-registro já estava escrito, mas agora por um motivo medido.

---

## Não vamos rodar a fase 3?

Não, e não por economia.

A fase 3 produz as covariáveis e o modelo nulo — idade, sexo, batch, componentes
principais, parentesco. **Nada disso depende do filtro de MAF.** São as mesmas
pessoas, o mesmo fenótipo, as mesmas covariáveis que o braço de variante rara usou.
Refazer daria exatamente os mesmos números, ao custo da parte caríssima (o parentesco,
1,45 bilhão de pares na coorte combinada).

Então este braço **herda**:

| o que | de onde |
|---|---|
| coorte (fase 1) | `../gene_burden_elena_replication/phase_1/` |
| covariáveis (fase 3) | `../gene_burden_elena_replication/phase_3/results/covariates/` |
| modelo nulo (step 1) | `../gene_burden_elena_replication/phase_4/results/step1/` |

E isso não é invenção nossa: **a Elena fez o mesmo.** Está na config dela,
`HL_only_ExWAS/saige_nextflow/saige_hearing_pipeline.config`:

```
step1_outdir = "${base}/HL_only_rarevariant/saige_step1"
```

O ExWAS dela leu o step 1 do braço de variante rara. Confirmação, não suposição.

Herdar tem um preço que está registrado no `PREMISES.md`: se aquela pasta mudar, este
braço fica desatualizado e nada avisa automaticamente.

---

## Aqui não entra filtro de VEP?

**Para rodar o teste, não entra.** E essa é justamente a diferença que define os dois
braços.

No braço de variante rara, a anotação é o que **constrói** o teste. O SAIGE-GENE não
testa variantes, testa genes — e para juntar variantes num gene é preciso decidir
quais entram. É aí que `pLOF` e `pDM` existem: são máscaras feitas de anotação
funcional. Sem VEP, não há o que agrupar, e não há teste.

Aqui cada variante é testada sozinha. Não há agrupamento, então não há nada para a
anotação decidir. O SAIGE testa as 70.383 (ou 74.553, ou 117.996) uma por uma, seja
ela sinônima, intrônica ou uma parada prematura.

**Mas a anotação entra depois, na interpretação — e é exatamente aí que está a
objeção do Doug.** Na reunião de 2026-10-02 ele disse:

> *"if the variants that associated were not pathogenic, then I just… I didn't
> understand the value of that […] from a biological standpoint, it's unclear what
> that would mean."*

E a Nikki, mais cedo na mesma reunião:

> *"for the ExWAS, these all mostly focused on common variants. And for the most part,
> we don't really see much of the damaging there anyway."*

Os dois estão falando da **mesma lacuna**: um ExWAS de variante comum encontra
associação estatística sem dizer se a variante faz alguma coisa. O teste não filtra
por função, então o resultado não vem com função.

Então eu vou fechar essa lacuna, e **de propósito já na fase 5, não como reação a
quem reclamar**: qualquer variante que passar uma das duas barras vai ser anotada com
VEP e AlphaMissense, e o resultado reportado junto do p-valor — consequência,
pontuação de dano, gene. Se as variantes associadas não forem danosas, isso vira um
número no documento em vez de uma pergunta sem resposta na reunião.

Já temos as ferramentas na LPC (VEP e ANNOVAR com dbNSFP 4.5a, e o Biofilter 4), e os
mesmos pipelines foram usados nos braços de variante rara. Vou acrescentar isso ao
`PREMISES.md` como etapa pré-registrada antes de rodar o step 2 — porque anotar só os
achados **depois** de ver quais foram os achados é o tipo de coisa que, feita sem
registro prévio, parece escolha conveniente.

---

## Arquivos

| | |
|---|---|
| `phase_2/data/{combined,EUR,AFR}.keep` | listas de amostras, conferidas contra o `.fam` |
| `phase_2/results/census/<coorte>.chr<N>.afreq.alt1.bins` | os 66 histogramas |
| `phase_2/results/census_summary.json` | o agregado, com o espectro completo |
| `phase_2/results/FINDINGS.md` | o mesmo conteúdo, em inglês, formato do repositório |

Os limites dos histogramas foram colocados em 0,0001, 0,001, 0,01, 0,05 e nos
espelhos. Se algum dia a gente quiser contar em outro limiar, é recontar esses
arquivos — não reler 600 GB de genótipos.

---

## Em aberto

- **Herança silenciosa.** Se `../gene_burden_elena_replication/phase_3/` mudar, este
  braço fica velho sem avisar. Não resolvi; está registrado em `PREMISES.md` P1–P3.
- **A anotação dos achados** entra na fase 5, pré-registrada. Ainda não escrita.
