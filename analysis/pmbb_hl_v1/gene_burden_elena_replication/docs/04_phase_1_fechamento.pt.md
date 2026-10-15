# O que fechamos na Fase 1

**Autor:** Andre Rico · **Data:** 2026-09-30
**Cópia de trabalho pessoal, em português.** Não é página de Confluence. O registro formal está em
[`phase_1/results/FINDINGS.md`](../phase_1/results/FINDINGS.md) e a referência completa em
[`02_phase_1_reference.md`](02_phase_1_reference.md).

> **Nota de terminologia.** Uso *tinnitus* em vez de *zumbido* para bater com o nome que aparece no
> código, nas colunas e nos diretórios do projeto (`tinnitus_PMBBv4_SAIGE.txt`, phecode `SO_397`).
> São a mesma coisa.

---

## A pergunta

Antes de perguntar *qual gene causa perda auditiva*, é preciso saber **quem tem perda auditiva** e
**quem está no estudo**. A Fase 1 verifica essas duas coisas. Tudo que vem depois herda a resposta —
e um erro aqui não tem conserto lá na frente.

---

## A tabela `observation` — a raiz de quase tudo

Isto explica três coisas diferentes mais adiante, então vale entender primeiro.

O PMBB entrega o prontuário eletrônico num formato padronizado chamado **OMOP**, criado para que
hospitais diferentes rodem a mesma análise sem reescrever código. Ele separa os dados por tipo:

| tabela | o que guarda | tamanho |
|---|---|---|
| `condition_occurrence` | **diagnósticos** | 5,6 GB |
| `observation` | **o resto** | 1,2 GB |
| `procedure_occurrence`, `drug_exposure`, `measurement_*` | procedimentos, remédios, exames | — |

A `observation` é a gaveta do "resto". Olhando o que mais aparece lá, são quase todos códigos **Z**
da CID — *fatores que influenciam o estado de saúde*, não doenças:

```
Z23      vacinação
Z79.899  uso prolongado de medicação
Z00.00   exame médico geral
Z94.0    status de transplante renal
Z12.31   rastreamento de câncer de mama
```

São **circunstâncias**, não diagnósticos. Faz sentido estarem separadas.

### Onde mora o problema

Os **phecodes** — o agrupamento de códigos que usamos como fenótipo — são construídos **só** a partir
da `condition_occurrence`. E no PMBB v4, por alguma decisão de mapeamento institucional, os códigos
padrão de tinnitus (388.3x, H93.1x) foram parar na `observation`:

```
tinnitus no phecode (condition_occurrence) :  3.016 eventos   ← só o pulsátil
tinnitus na observation                    : 25.094 eventos   ← 89% do total
```

**Quem monta o fenótipo pelo caminho normal enxerga 11% do tinnitus que existe.** E nada avisa — não
há erro, não há mensagem, o phecode de tinnitus existe e retorna dados. São só os dados errados.

Isso não está documentado em lugar nenhum do release. Quem descobriu foi a **Nikki**.

---

## O que está certo

**A regra que decide quem é caso reproduz exatamente.**

Reconstruí casos e controles do zero, a partir dos dados brutos, sem olhar o código dela. Deu igual:
**70.925 de 70.925 pessoas, uma a uma.** Não é só o total bater — são as mesmas pessoas.

Cobre a regra de 2 (diagnóstico em duas datas), a exclusão do "meio ambíguo", e o tratamento do
tinnitus da `observation`.

**A estratificação por ancestralidade também confere.** Os arquivos EUR e AFR contêm só quem o
release classifica assim, e o `combined` é exatamente a soma com os grupos menores.

### Onde eu errei e o pipeline acertou

Minha primeira reconstrução discordou em **558 pessoas**, porque eu montei a evidência auricular só
com o `condition_occurrence` e esqueci a `observation`. Todas as 558 tinham tinnitus lá.

Guardar esse número — ele reaparece adiante de um jeito que fecha a história.

### Uma regra que costuma ser mal-entendida

**Ter tinnitus não exclui ninguém.** 1.815 dos nossos 6.752 casos têm tinnitus também, e são casos
normalmente.

O que exclui é ter doença de ouvido **e nenhum diagnóstico de perda auditiva**. Esse grupo — o "meio
ambíguo", 9.411 pessoas — é majoritariamente outra coisa:

```
9.411  meio ambíguo
   558    casos de tinnitus
 8.853    outra evidência auricular (otite, cerume, vertigem…)
```

---

## O que está errado

### Achado 1 — 431 pessoas analisáveis ficaram de fora

O estudo roda sobre **exoma**. Mas a lista de participantes foi filtrada por quem tinha dado de
**array/imputado** — outro exame, que a análise não usa em lugar nenhum.

```python
fam_file = ".../Imputed/common_snps_LD_pruned/...ldpruned.ALL.fam"
# Keep only samples with genotype data
matched = pheno[pheno["PMBB_ID"].isin(fam["IID"])]
```

O comentário tem a intenção certa. O arquivo consultado é o errado — e o LD-pruned de **exoma**
existia, construído pelo mesmo pipeline.

Resultado: **431 pessoas** com exoma completo fora, **40 delas com perda auditiva**.

Não é divergência de método. É **inconsistência interna**: o pipeline barra com base num dado que
nenhuma etapa dele consome.

### Achado 2 — 556 pessoas que a regra exclui entraram como controles

E aqui a `observation` volta. Os 556:

```
evidência auricular NÃO-tinnitus    :    0
tinnitus no phecode (SO_397)        :    0
tinnitus na tabela OBSERVATION      :  556  ← 100%
```

A **única** evidência auricular deles está naquela tabela.

E aqui está a parte contraintuitiva: **não foi por consultar a `observation` que eles entraram — foi
por não consultar.** A versão do fenótipo usada não lia aquela tabela, então o tinnitus dessas 556
pessoas simplesmente não existia do ponto de vista do código. Elas apareciam sem nenhuma doença de
ouvido, e a regra que separa o meio ambíguo não tinha do que se agarrar.

Ninguém as viu e decidiu mantê-las como controles. **Ninguém as viu.**

A correção das 20:57 foi justamente adicionar essa fonte — e aí elas passaram a aparecer como
excluídas, corretamente.

E o fecho: as 556 são **exatamente** as pessoas que meu check 02 classificou errado pelo mesmo
motivo. Interseção de 100%. Cometi o mesmo erro, três meses depois, sozinho.

### Achado 3 — a correção existiu, e um dos ramos não a pegou

```
31/jul 20:21   covariáveis de perda auditiva construídas        (57.632)
31/jul 20:57   fenótipo REGENERADO — a correção da observation   (57.080)
01/ago 14:27   covariáveis de TINNITUS refeitas da versão corrigida  ✓
03/ago 14:38   modelos ajustados, nas covariáveis VELHAS de HL
03/ago 16:22   o teste roda
```

A correção veio **primeiro**. A análise rodou **três dias depois**, sobre um arquivo construído antes
dela.

E dá para ver que a regeneração foi percebida: as covariáveis de tinnitus foram refeitas no dia
seguinte, da versão corrigida. O ramo de perda auditiva simplesmente nunca foi refeito.

Então a análise de perda auditiva rodou sobre um fenótipo **anterior à descoberta do tinnitus no v4**
— a mesma descoberta pela qual o pipeline merece crédito.

Não é descuido. É o modo de falha mais comum em pipeline com etapas manuais: dois ramos paralelos, um
atualizado e o outro não, sem nada declarando que um depende do outro.

### Achado 4 — quem ficou de fora não é amostra aleatória

```
taxa de exclusão na coorte inteira     0,73%

Leste Asiático   208 / 1.333  =  15,60%   ← 21× a taxa geral
Europeu          233 / 51.867 =   0,45%
Africano          49 / 14.927 =   0,33%
```

**Uma em cada seis pessoas do Leste Asiático foi excluída.** Também há desvio em sexo (mulheres 1,01%
contra homens 0,43%) e forte em lote (92% dos excluídos vêm do lote 1); os excluídos são mais jovens
e recrutados mais cedo.

O padrão sugere causa técnica — sequenciados por exoma cedo, array nunca concluído para uma parte.
Nada sobre as pessoas.

Mas o efeito é de **representatividade**, não de poder estatístico: um grupo que já era 1,6% da
coorte perdeu um sexto dos seus. Isso não melhora aumentando o N.

---

## O fluxo dos números

```
70.925   coorte de exoma no PMBB v4
           −  4.007  excluídos: evidência em uma data só
           −  9.411  excluídos: meio ambíguo (doença de ouvido, sem perda auditiva)
57.507   ANALISÁVEIS — o que o fenótipo sustenta          ← nosso número
           −    431  cortados pelo filtro do .fam imputado (40 casos)
           +    556  pessoas cujo tinnitus a versão usada NÃO ENXERGAVA,
                      e que por isso entraram como controles
57.632   O QUE A ESTATÍSTICA RODOU                        ← número da Elena
```

O sinal `+556` engana à primeira vista: ninguém *acrescentou* essas pessoas. Elas **nunca foram
retiradas**, porque a regra que as retiraria depende de uma fonte de dados que aquela versão do
fenótipo ainda não lia.

### Duas armadilhas de leitura

**57.080 não é 57.632 menos nada.**

```
57.632 − 556  =  57.076
57.080 −   4  =  57.076   ← o núcleo comum
```

São **duas derivações independentes**, de versões diferentes do fenótipo. Por isso o 57.080 tem 4
pessoas que o 57.632 não tem — se fosse subtração, seria impossível. E o 57.080 **nunca foi
consumido por nada**: serve só como prova de que o arquivo foi regenerado.

**Os dois conjuntos não se contêm.** Nenhum é subconjunto do outro. A rodada é menor num sentido
(−431) e maior no outro (+556).

---

## O que fechamos: dois cenários, em arquivo

| | combined | EUR | AFR |
|---|---|---|---|
| **reprodução** (o que rodou) | 57.632 · 6.712 casos | 43.016 · 5.160 | 11.387 · 1.279 |
| **corrigido** (o que deveria) | 57.507 · 6.752 casos | 42.786 · 5.183 | 11.334 · 1.285 |

Arquivos em `phase_1/results/cohorts/`, gerados por `scripts/05_emit_cohorts.py`.

**Por que dois e não só o corrigido.** Se rodássemos só a versão corrigida, qualquer diferença nos
resultados seria ambígua: veio do defeito, ou da coorte que *nós* mudamos? Perderíamos a capacidade
de dizer "reproduzimos os números dela" — que é o produto principal deste trabalho.

**A diferença entre os dois é o resultado.** Transforma "existe um defeito" em "o defeito move este
gene nesta magnitude" — ou em "o defeito é real e não muda nada", que é uma recomendação diferente.

⚠️ **Um detalhe que precisa ficar guardado.** O braço de reprodução foi **copiado**, não
reconstruído, porque o arquivo de onde ele veio foi sobrescrito. O arquivo de covariáveis
sobrevivente é o **único registro do que foi analisado**. Se ele for apagado, esse cenário se perde
para sempre.

O corrigido tem um número bom de contar: ele **restaura 165 participantes do Leste Asiático**,
fazendo aquele grupo crescer 17,6% enquanto todos os outros mudam meio por cento.

---

## O que fica em aberto

| | |
|---|---|
| Por que o filtro usou o `.fam` do array | pergunta para Nikki e Elena |
| Por que o ramo de perda auditiva não foi refeito, se o de tinnitus foi | idem |
| Os 4 controles legítimos que somem sem explicação | resíduo pequeno, não perseguido |
| Tinnitus como fenótipo | fora de escopo — verificamos só perda auditiva |
| O fenótipo combinado HL-ou-tinnitus da decisão de 2026-07-01 | **nunca foi construído** neste pipeline |
| Se o grupo restaura os 431 | decisão deles, não nossa |

Sobre o fenótipo combinado: quando for construído, as **934 pessoas** que hoje são casos de tinnitus
sem perda auditiva mudam de lado — deixam de ser excluídas e viram casos. É cerca de **+14%** no
número de casos.

Mas note que isso **não resolve o meio ambíguo**: 8.853 das 9.411 pessoas ali não têm nada a ver com
tinnitus, e continuariam excluídas.

---

## Em duas frases

O jeito de decidir quem tem perda auditiva está **certo** — reproduzido pessoa a pessoa. O jeito de
decidir quem entra no estudo está **errado** em três pontos independentes, todos pequenos, todos
empurrando o resultado para o lado conservador, e **nenhum detectado por nada** no pipeline.

O problema não é 0,59% dos casos ter sumido. É que sumiu e o pipeline reportou sucesso.
