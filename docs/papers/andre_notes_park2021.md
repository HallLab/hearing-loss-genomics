# Notas de Estudo — Métodos do Park et al. 2021 (ExoPheWAS)

> **Para quê:** referência pessoal do Andre para fixar os conceitos de **seleção de variantes e desenho do rare-variant burden test**, partindo dos Métodos do Park et al. 2021 (Nat Med). Escrito em linguagem didática, com analogias de software/PRS.
>
> **Paper:** Park J et al. (2021). *Exome-wide evaluation of rare coding variants using EHR identifies new gene-phenotype associations.* Nat Med 27(1):66–72. doi:10.1038/s41591-020-1133-8. PDF: `docs/papers/nihms-1768013.pdf`. Resumo formal: `docs/papers/paper_summary_park2021.md`.
>
> **Por que este paper importa pro projeto:** é a **fonte publicada do achado ZNF175 → tinnitus** (p=3,24×10⁻¹⁰, replicado em BioMe/DiscovEHR/UKB; ortólogo de camundongo *Zfp719* é surdo). Não confundir com o Hui 2023 (PLOS Genetics), que **não** menciona ZNF175.

---

## Índice
1. [pLOF — o que "desliga" um gene](#1-plof)
2. [ANNOVAR e RefSeq — quem anota e com que mapa](#2-annovar-refseq)
3. [Ferramentas de hoje: VEP, LOFTEE, gnomAD v4](#3-ferramentas-modernas)
4. [O que é o gnomAD (e o que ele NÃO é)](#4-gnomad)
5. [Os DOIS filtros: frequência + função](#5-dois-filtros)
6. [MAF ≤ 0,1% — o filtro de raridade](#6-maf)
7. [Por que RARAS e não comuns](#7-raras)
8. [REVEL ≥ 0,5 — a missense ambígua (e por que fica fora do discovery)](#8-revel)
9. [≥25 portadores → 1.518 genes — o funil do lado-gene](#9-portadores)
10. [1.000 phecodes (≥20 casos) — o funil do lado-fenótipo](#10-phecodes)
11. [Modelo estatístico — a regressão por trás de cada célula](#11-modelo)
12. [O desenho inteiro numa imagem](#12-panorama)
13. [Glossário rápido](#13-glossario)

---

<a name="1-plof"></a>
## 1. pLOF — o que "desliga" um gene

**pLOF = predicted Loss-of-Function** = variantes que, **segundo a anotação**, provavelmente quebram o gene (proteína não produzida ou truncada/inútil). "Predicted" = é **inferência por regra**, não confirmação de laboratório.

**Pré-requisito mental — o frame de leitura:** o DNA é lido em **trincas (códons)** de 3 letras; cada códon → 1 aminoácido; há códons de **STOP**. Analogia: um parser que lê um stream em blocos de 3 bytes, com um byte terminador.

Os **4 tipos de pLOF**:

| Tipo | O que faz | Analogia de código |
|---|---|---|
| **Frameshift indel** | inserção/deleção **não múltipla de 3** → desalinha TODOS os códons a jusante | off-by-one num parser de blocos: o resto vira lixo |
| **Stop gain** (nonsense) | troca de letra cria um **STOP prematuro** → proteína truncada | `return`/`exit()` injetado no meio da função |
| **Stop loss** | mutação **destrói o STOP** → tradução continua além do fim | apagaram o `break`/sentinela do loop |
| **Splice-site (canônico)** | atinge os **2 nucleotídeos canônicos** (GT…AG) das bordas do íntron → splicing erra (pula éxon / mantém íntron) | corromper o delimitador que o parser usa para cortar segmentos |

**Por que splicing importa:** genes têm **éxons** (viram proteína) e **íntrons** (removidos). A maquinaria reconhece as bordas dos íntrons por sinais quase invariantes — os mais conservados são 2 nucleotídeos em cada ponta. Mutou ali → corte errado → proteína corrompida.

> **Em uma frase:** pLOF = variante que, pela regra biológica, **desliga o gene** — frameshift, stop gain/loss, ou quebra do splice canônico.

---

<a name="2-annovar-refseq"></a>
## 2. ANNOVAR e RefSeq — quem anota e com que mapa

Não são tipos de variante — são **ferramenta + referência**:

- **ANNOVAR** = software de **anotação**. Recebe a variante e responde "isto é frameshift / stop gain / …". Como um **linter/static-analyzer** que rotula cada mudança. (Aplica as regras da seção 1 automaticamente em milhões de variantes.)
- **RefSeq** = **catálogo de referência** (NCBI) com a "planta baixa" dos genes: onde começam/terminam, éxons, íntrons. O ANNOVAR precisa dela para saber **onde** a variante cai. Alternativas: Ensembl/GENCODE, UCSC — a escolha pode mudar levemente os rótulos (por isso registra-se qual foi usada).

---

<a name="3-ferramentas-modernas"></a>
## 3. Ferramentas de hoje: VEP, LOFTEE, gnomAD v4

As **definições conceituais** dos 4 tipos continuam idênticas. O que mudou desde 2021:

| Papel | Park 2021 | Equivalente moderno |
|---|---|---|
| Anotar consequência | ANNOVAR | **VEP** (termos padronizados Sequence Ontology: `frameshift_variant`, `stop_gained`, …) |
| Planta baixa do gene | RefSeq | RefSeq **ou** GENCODE/Ensembl (gnomAD usa GENCODE) |
| Filtrar LoF "de verdade" | *(não tinha — pegava todos)* | **LOFTEE** (plugin do VEP) |
| Frequência populacional (MAF) | gnomAD **v2** (141k, GRCh37) | **gnomAD v4** (~807k, GRCh38) + constraint (LOEUF) |

**LOFTEE** é o pulo do gato moderno: o VEP sozinho ainda é "burro" (diz a categoria, não a confiança). O LOFTEE pega as variantes rotuladas como LoF e separa em:
- **HC (high-confidence)** — provável LoF de verdade;
- **LC (low-confidence)** — descarta com flag (ex.: stop nos últimos éxons, splice em região fraca).

→ **pLOF "sério" de hoje = VEP + LOFTEE, ficando só com HC.** É **mais conservador** que o critério do Park ("qualquer frameshift/stopgain/splice"). Atenção ao comparar contagens de portadores: menos variantes sobrevivem.

**Decisão de projeto em aberto:** re-anotar com VEP+LOFTEE/gnomAD v4 (GRCh38) **ou** manter ANNOVAR/gnomAD v2 para ficar fiel à replicação do Daniel? Trade-off entre estado-da-arte × fidelidade. *(Lembrete: build do genoma — v4 é GRCh38; v2 é GRCh37; precisa de liftover se misturar.)*

---

<a name="4-gnomad"></a>
## 4. O que é o gnomAD (e o que ele NÃO é)

O gnomAD é um **catálogo agregado de variantes** que serve a dois papéis distintos no nosso fluxo. **Ponto crítico: ele NÃO expõe os indivíduos** — entrega **estatísticas agregadas**, igual ao individual-level protegido do PMBB.

Por variante, ele guarda:
- **AC** (allele count) — quantas cópias do alelo apareceram;
- **AN** (allele number) — total de alelos observados (≈ 2× nº de pessoas com cobertura);
- **AF** = AC/AN → **a frequência** (de onde sai o MAF);
- nº de **homozigotos**.

Duas coisas guardadas lado a lado, mas de naturezas diferentes:
- **Frequência** → **varia por população** (AC/AN/AF por grupo de ancestralidade: AFR, NFE, Latino, EAS, SAS, …, + global e grpmax/popmax). É por isso que o método usa **limiar de MAF por ancestralidade**.
- **Anotação funcional** (consequência VEP, flag LOFTEE, aminoácido) → **igual para todos**, porque depende só de *onde* a variante cai e *o que* muda, não de quem a carrega.
- (Por gene, ainda guarda **constraint**: LOEUF/pLI = "o quanto a população tolera perder esse gene".)

> **Em uma frase:** gnomAD = "para cada variante, **quão rara em cada ancestralidade** (frequência) e **o que provavelmente faz** (anotação igual p/ todos)", sem nunca expor indivíduos.

---

<a name="5-dois-filtros"></a>
## 5. Os DOIS filtros: frequência + função

Erro comum: achar que o burden filtra **só** por frequência. São **dois filtros agindo juntos** (a interseção):

```
milhões de variantes por indivíduo
        │
        ├── filtro 1: rara?      (MAF ≤ 0,1%)
        │
        ├── filtro 2: funcional? (pLOF, ou missense REVEL≥0,5)
        │
        └──► sobra um punhado por gene  ──►  agrega (burden) por gene
```

- Rara **mas inócua** → não entra.
- pLOF **mas comum** → não entra no burden (vai p/ análise de variante única).
- Entra só: **rara E provavelmente quebra o gene.**

---

<a name="6-maf"></a>
## 6. MAF ≤ 0,1% — o filtro de raridade

**MAF = frequência do ALELO, não de pessoas.** Como cada um tem 2 cópias e raras quase nunca são homozigotas:

> **% de portadores ≈ 2 × MAF.** Logo, MAF 0,1% → ~**0,2% das pessoas** carregam.

Concreto (coorte Park, 10.900 pessoas ≈ 21.800 alelos): `MAF ≤ 0,1%` → no máximo ~**22 cópias** na coorte inteira. **Por isso** existe o piso de ≥25 portadores e a necessidade de **agregar** (seção 9): variante rara sozinha = zero poder.

**Comuns não são descartadas** — o 0,1% é um **divisor de águas**:
- **Raras (≤0,1%)** → agregadas em **burden por gene**;
- **Comuns (>0,1%)** → testadas **uma a uma** (univariate), pois já têm gente suficiente.

**Implementação (conferir no runbook do Daniel):** costuma haver **dois** cortes de frequência — o externo (gnomAD, por ancestralidade) **e** o interno da coorte. O Hui 2023 usou `gnomAD MAF < 0,001` **E** `cohort MAF < 0,01`.

---

<a name="7-raras"></a>
## 7. Por que RARAS e não comuns

**A relação efeito × frequência** (decorre da **seleção natural**):

```
efeito
grande │  ● raras (pLOF, mendelianas)
       │    ● ●
       │       ● ●
pequeno│           ● ● ● ● ●  comuns (GWAS/PRS)
       └─────────────────────────────► frequência
```

**Seleção purificadora:** variante que quebra gene importante e causa doença séria → atrapalha o portador → seleção **remove** ao longo das gerações → **não consegue ficar comum**. Então:

> Se uma variante é **rara E quebra o gene**, a raridade já é pista de que ela importa. Filtrar `MAF ≤ 0,1% + pLOF` **enriquece em efeito biológico forte**.

**Contraste com o seu mundo (PRS/GWAS):** são as duas pontas opostas do espectro.

| | GWAS / PRS | Rare-variant burden |
|---|---|---|
| Frequência | comuns (>1–5%) | raras (≤0,1%) |
| Efeito/variante | pequeno (OR ~1,01–1,2) | grande (OR 2–10+) |
| Onde caem | maioria **não-codificante** | **codificante** (quebram proteína) |
| Estratégia | somar **milhares** → escore | agregar **poucas** por gene → burden |
| Captura | risco **poligênico** difuso | efeito quase-**mendeliano** de 1 gene |

São lentes **complementares** (por isso o Hui 2023 fez os dois). E ambos os escores (PRS e burden) tinham **baixo poder preditivo individual** no Hui — coerente: HL adulta é majoritariamente poligênica/ambiental, o componente raro forte só aparece num subconjunto.

**Bônus: interpretabilidade.** Hit comum de GWAS cai em região não-codificante → "tem algo por aqui", mas não se sabe o gene nem o mecanismo (precisa de fine-mapping). Um **pLOF dá a hipótese causal de graça**: "essa proteína não foi feita → o gene X é necessário". É o caso do **ZNF175**: sinal de variantes que **desligam o gene** + biologia confirmando (camundongo *Zfp719* surdo, expressão em células ciliadas).

> **Em uma frase:** queremos as raras porque a seleção mantém raras justamente as de **efeito forte que quebram genes** — sinal grande e **mecanicamente legível**, invisível ao GWAS/PRS. Agregar por gene resolve o problema de poder da raridade.

---

<a name="8-revel"></a>
## 8. REVEL ≥ 0,5 — a missense ambígua (e por que fica fora do discovery)

**Missense** = troca de 1 letra que **troca 1 aminoácido**. A proteína é feita, do tamanho certo, mas com "uma letra errada". É o **meio ambíguo** do espectro: algumas devastadoras, a maioria inócua, e **não dá pra saber pela classe**.

**REVEL** (*Rare Exonic Variant Ensemble Learner*) = escore **0 a 1** de patogenicidade **só para missense**. É um **ensemble** (ML) que combina ~13 preditores (SIFT, PolyPhen-2, conservação…). **≥0,5** é só o **limiar de corte** escolhido ("considero deletéria se ≥0,5").

| | pLOF | Missense + REVEL |
|---|---|---|
| Natureza | **determinística** (regra biológica) | **probabilística** (modelo *chuta*) |
| Confiança | alta, homogênea | variável; erra em variantes individuais |
| Efeito | todos "desligam" (mesma direção) | heterogêneo (forte/fraco/às vezes oposto) |

### Por que o REVEL NÃO está no Discovery burden (é de propósito)

1. **Validação ortogonal (o motivo mais elegante):** o Park **descobre com pLOF** e **confirma com missense REVEL≥0,5**, no mesmo PMBB mas em **portadores não-sobrepostos** (variantes/pessoas diferentes). Se o gene é real, **duas evidências independentes** apontam junto. Misturar tudo no discovery destruiria essa independência. → filosofia **DiCE** (Diverse Convergent Evidence).
2. **Discovery limpo × escala:** ~1,5M de testes (seção 10); o REVEL é ruidoso (erra em variantes individuais) e inflaria falso-positivo. Descobre-se com o sinal mais nítido (pLOF).
3. **Premissa do burden = efeito homogêneo:** burden assume que todas as variantes colapsadas têm efeito parecido e na mesma direção. pLOF satisfaz; missense é heterogênea e **dilui** o sinal.

### ⚠️ Diferença que importa pro projeto
Não é regra universal — é escolha do Park. **O Hui 2023 fez DIFERENTE:**

| | Park 2021 | Hui 2023 |
|---|---|---|
| Papel da missense | **só robustness** (separada) | **dentro do burden primário** |
| Corte REVEL | ≥ **0,5** | > **0,6** |

Como a replicação de vocês segue o **Hui** (runbook do Daniel), vale o critério **REVEL > 0,6 junto com pLOF**, não o esquema do Park. **Conferir no runbook** o corte e o papel reais.

---

<a name="9-portadores"></a>
## 9. ≥25 portadores → 1.518 genes — o funil do lado-gene

**Heterozigoto** = variante em **1** das 2 cópias do gene (a outra normal). Para raras pLOF, **quase todo portador é heterozigoto** (homozigoto raro é raríssimo/severo). Então "≥25 heterozygous carriers" = **≥25 pessoas** na coorte carregam ≥1 pLOF rara naquele gene.

**Por que 25?** Piso de **poder estatístico**: com poucos portadores, o grupo "exposto" é minúsculo e o intervalo de confiança explode — teste cego. O Park fez análise formal de poder (Ext. Data Fig. 2).

**Coincidência elegante** (10.900 pessoas ≈ 21.800 alelos):
$$\frac{25}{21.800} \approx 0{,}00115 \approx 0{,}1\%$$
O piso de 25 portadores ≈ burden do gene a **~0,1%** — mesmo patamar do MAF. O método inteiro é calibrado nesse nível de raridade.

**O grande recado — o funil:**
```
~20.000 genes codificantes
        │  ≥25 portadores de pLOF NA COORTE de 10.900
        ▼
   1.518 genes testáveis   ◄── só ~8% do exoma!
```
**>90% dos genes não foram testáveis**, por dois motivos:
1. **Coorte pequena** (10.900) → muitos genes ficam abaixo de 25 só por falta de gente.
2. **Constraint (a ironia):** genes **essenciais/restritos** (LOEUF baixo) → seleção remove pLOFs → **quase ninguém carrega** → <25 portadores → **caem fora**. Genes que **toleram** perder função acumulam pLOFs → testáveis.

> Paradoxo: o burden enxerga melhor os genes **menos restritos** e fica **cego** para muitos dos mais críticos. Os 1.518 são um **subconjunto enviesado**.

**ZNF175 sob essa luz:** para ter sido **testado** (≥25 portadores), é um gene que **tolera razoavelmente** perda de função → combina com penetrância incompleta / "modificador + segundo hit".

**Por que isso justifica a v3:** coorte maior → mais portadores/gene → **mais genes cruzam o piso** (inclusive restritos antes invisíveis) + mais poder no ZNF175. 10.900 (Park) → ~40k (Hui) → **PMBB v3** (vocês). Esperar contagens de portadores mudarem na migração v2→v3.

*Nota:* 25 é escolha, não lei (outros usam 10, 20). Trade-off cobertura × confiabilidade. **Conferir o piso no runbook.**

---

<a name="10-phecodes"></a>
## 10. 1.000 phecodes (≥20 casos) — o funil do lado-fenótipo

**O problema:** o EHR registra **códigos ICD** (faturamento) — hiper-granulares, ruidosos, aos milhares.

**Phecode** = **agrupamento curado** de ICD em uma doença clinicamente significativa (Phecode Map 1.2; ~1.800 phecodes). Analogia: ICD = enums/log crus inconsistentes; **phecode = camada de normalização** que agrupa em categorias úteis (tipo normalizar mil user-agents em "Chrome/Firefox/Safari"). É um `GROUP BY doença` por cima do faturamento.

**Caso/controle (regras anti-ruído):**
- **Caso** = código em **≥2 datas distintas** ("regra do 2" — 1 vez só costuma ser *rule-out*/erro).
- **Controle** = **nunca** teve o código.
- **1 data só** = **excluído**.
- **Exclusões de controle:** remove dos controles quem tem condições *relacionadas* (ex.: outras doenças do ouvido fora do controle de surdez).

**Piso `≥20 casos`** = gêmeo (lado-fenótipo) do `≥25 portadores`. Sobraram **1.000 phecodes** (de ~1.800).

**Os dois pisos do funil:**
| Eixo | Filtro de poder | Sobrou |
|---|---|---|
| Gene (X) | ≥25 portadores pLOF | 1.518 genes |
| Fenótipo (Y) | ≥20 casos | 1.000 phecodes |

**A grade e múltiplos testes:**
$$1.518 \times 1.000 \approx 1{,}5 \text{ milhão de testes (regressões logísticas)}$$
Com 1,5M testes, `p<0,05` daria ~75 mil falsos-positivos. Daí o limiar **`p < 10⁻⁶`** — escolhido **empiricamente** no ponto em que o QQ-plot se descola do esperado (não Bonferroni ingênuo ≈3×10⁻⁸). Eles **compensam a leniência com replicação** (DiCE): validade vem "tanto da replicação quanto do limiar".

**Amarração com o projeto:**
- **Hearing loss é phecode** (Hui: **phecode 389**; Park: HL N=579, 5,3%). **Tinnitus é phecode separado.**
- O ZNF175 deu sinal **no phecode de tinnitus**, com **hearing loss logo abaixo** do 10⁻⁶ → são **células diferentes da grade**, não contradição.
- **Controle vazado** (Hui: 27% dos "controles" por phecode tinham HL no audiograma) = falha clássica da definição ICD: gente com HL real sem o código vai parar no controle. **Por isso** o Hui cruzou com **audiograma (PTA>25 dB)** e o projeto prioriza os ~4.000 com audiograma+exoma — fenótipo quantitativo de verdade.

> **Em uma frase:** phecodes normalizam ICD cru em doenças usáveis; `≥20 casos` é o piso do lado-fenótipo; a grade de ~1,5M testes força `p<10⁻⁶` + replicação; e a imperfeição do caso/controle por ICD é o que motiva os audiogramas no projeto.

---

<a name="11-modelo"></a>
## 11. Modelo estatístico — a regressão por trás de cada célula

Cada célula da grade (1 gene × 1 phecode) é **um** modelo. A frase dos Métodos empacota **5 decisões**.

### A equação
$$\text{logit}\,P(\text{caso}) = \beta_0 + \underbrace{\beta_{b}\cdot G_{\text{burden}}}_{\text{o que interessa}} + \beta_1\,\text{idade} + \beta_2\,\text{idade}^2 + \beta_3\,\text{sexo} + \sum_{i=1}^{10}\gamma_i\,\text{PC}_i$$

Tudo gira em torno de **um coeficiente: $\beta_b$**. O p-valor testa $H_0:\beta_b=0$; $\exp(\beta_b)=\text{OR}$. O resto são covariáveis de ajuste.

### 1. Por que logística
Desfecho **binário** (caso/controle) → modela o **log-odds**; $\exp(\beta)$ = odds ratio. *(Para desfecho quantitativo — ex.: o **audiograma/PTA** do projeto — troca-se por **regressão linear**, mesma estrutura de covariáveis.)*

### 2. O preditor: "additive, fixed threshold"
- **Fixed threshold:** variante "conta" se passa nos cortes fixos (MAF≤0,1% + pLOF); burden = **soma das qualificadas** no gene (método CAST/burden clássico).
- **Additive:** cada alelo a mais soma **linear** no log-odds (0/1/2); como raras quase nunca são homozigotas, na prática é ~**0/1**.
- ⚠️ O burden **assume efeito homogêneo** (todas mesma direção) → é **por isso** que só pLOF entra (missense heterogênea diluiria). Contraste: **SKAT** (variance-component) tolera direções opostas; burden simples **não** (soma → opostos se cancelam).

### 3. Covariáveis — controle de confundidores
| Covariável | Por quê |
|---|---|
| **idade + idade²** | doença idade-dependente (HL ↑ com idade); o **quadrático** ajusta uma **curva** (efeito acelera). Também: idoso = mais consultas → mais ICD → confundidor de *ascertainment* |
| **sexo** | prevalência/ascertainment diferem por sexo |
| **10 PCs** | **estratificação populacional** (clássico do PRS): se freq. alélica E prevalência variam por subpopulação → associação falsa. PCs capturam estrutura fina e removem o confundimento |

### 4. EUR + AFR separados → meta-análise IVW
Por que não pôr todos juntos só com PCs? EUR (8.198) e AFR (2.172) têm **frequências alélicas muito distintas** → estratificar é mais limpo que confiar só nos PCs.
- **Estratifica:** roda o mesmo modelo **dentro** de EUR e de AFR (cada um com seus PCs) → pares $(\beta_k, SE_k)$.
- **Combina (IVW, efeito fixo):**
$$\beta_{\text{meta}}=\frac{\sum_k \beta_k/SE_k^2}{\sum_k 1/SE_k^2},\qquad SE_{\text{meta}}=\sqrt{\frac{1}{\sum_k 1/SE_k^2}}$$
peso = **precisão** ($1/SE^2$); o estrato maior/mais preciso domina (geralmente EUR).
- Mnemônico: **estratificar = controle grosso da ancestralidade; PCs = controle fino dentro do estrato.**
- As coortes de replicação (PMBB2, BioMe, DiscovEHR, UKB) entram como **mais um $(\beta_k,SE_k)$** nesse mesmo somatório.

### 5. Firth penalized likelihood
- **Problema:** raras + binário → **separação** ("0 casos entre portadores" ou "todos casos") → ML padrão explode ($\hat\beta\to\pm\infty$, SE estoura, p-valor inflado).
- **Conserto:** Firth adiciona penalização (prior de Jeffreys) que **encolhe** $\hat\beta$ → estimativa **finita e estável**. Padrão em baixa contagem (embutido em SAIGE/REGENIE/PLINK).
- **Validação:** $\beta$/p batem com **regressão logística exata** (padrão-ouro caro) → penalização não distorceu.

### Saídas
$\beta_b$ (log-OR) · $\exp(\beta_b)$ (OR) · p-valor (vs. $10^{-6}$) · **DOE** = sinal de $\beta_b$ (triângulos ↑/↓ da Fig. 2).

### Contexto moderno
Park (R 3.3.1) = logística simples + **remoção** de aparentados + Firth. Hoje: **REGENIE / SAIGE** resolvem de uma vez (1) desequilíbrio caso-controle (saddlepoint/SPA), (2) **parentesco** (modelos mistos, em vez de remover), (3) escala. *Radar p/ v3: se crescer muito e tiver aparentados, "remover 3º grau + logística simples" pode virar modelo misto. Conferir no runbook do Daniel.*

> **Em uma frase:** cada gene×phecode é uma **logística** com um único coeficiente de interesse (o **burden aditivo**, premissa de efeito homogêneo); ajusta **idade/idade²/sexo** (confundidores clínicos e de ascertainment) e **10 PCs** (estrutura fina); roda **EUR e AFR separados → meta-análise IVW** (controle grosso de ancestralidade); usa **Firth** para não explodir com poucos portadores.

---

<a name="12-panorama"></a>
## 12. O desenho inteiro numa imagem

```
                        DISCOVERY (PMBB, 10.900)
  ┌───────────────────────────────────────────────────────────────┐
  │  VARIANTES                          FENÓTIPOS                    │
  │  milhões                            ICD cru                      │
  │     │ filtro freq: MAF≤0,1%(gnomAD)    │ → phecodes (Map 1.2)     │
  │     │ filtro função: pLOF              │ caso=≥2 datas            │
  │     ▼                                  ▼ controle=nunca           │
  │  agrega por gene                    ≥20 casos                     │
  │  ≥25 portadores → 1.518 genes  ×  1.000 phecodes                  │
  │            └──────────── ~1,5M regressões logísticas ────────────│
  │            ajuste: idade, idade², sexo, 10 PCs                    │
  │            EUR e AFR separados → meta-análise                     │
  │            limiar: p < 10⁻⁶  → 97 genes                           │
  └───────────────────────────────────────────────────────────────┘
                              │
        ROBUSTNESS / REPLICATION (validação ortogonal — DiCE)
        • mesmo PMBB: missense REVEL≥0,5 + variantes únicas (portadores não-sobrepostos)
        • outras coortes: PMBB2(AFR), BioMe, DiscovEHR, UKB; BioVU (single-variant)
                              │
                              ▼
                   26 genes robustos (5 controles + 21 novos)
                   inclui  ZNF175 → tinnitus
```

---

<a name="13-glossario"></a>
## 13. Glossário rápido

| Termo | Definição em 1 linha |
|---|---|
| **pLOF** | variante que (por regra) desliga o gene: frameshift, stop gain/loss, splice canônico |
| **Frameshift** | indel não-múltiplo de 3 → desalinha o frame de leitura |
| **Stop gain / loss** | cria STOP prematuro / destrói o STOP |
| **Splice site canônico** | os 2 nucleotídeos (GT…AG) das bordas do íntron |
| **Missense** | troca de 1 aminoácido (proteína inteira, "1 letra errada") |
| **ANNOVAR / VEP** | softwares que **anotam** a consequência da variante |
| **RefSeq / GENCODE** | catálogos ("planta baixa") dos genes usados na anotação |
| **LOFTEE** | plugin do VEP que separa LoF **HC** (confiável) de **LC** (descarta) |
| **REVEL** | escore ML 0–1 de patogenicidade de **missense** (corte ≥0,5 no Park, >0,6 no Hui) |
| **gnomAD** | catálogo **agregado** de frequência por ancestralidade + anotação (v2=141k/GRCh37; v4=807k/GRCh38) |
| **MAF** | frequência do **alelo** (≈ ½ da % de portadores p/ raras) |
| **Burden test** | agrega variantes raras de um gene e testa "gene quebrado vs não" |
| **LOEUF / pLI** | métrica de **constraint**: o quanto a população tolera perder o gene |
| **ICD / phecode** | código de faturamento cru / agrupamento curado em doença de pesquisa |
| **Heterozigoto** | variante em 1 das 2 cópias do gene |
| **DiCE** | Diverse Convergent Evidence: empilhar evidências independentes p/ confirmar |
| **Firth** | regressão logística penalizada, robusta p/ poucos casos (conserta separação) |
| **PC (principal component)** | eixo que captura estrutura de ancestralidade (ajuste do modelo) |
| **Separação (separation)** | célula caso/controle vazia → ML explode ($\hat\beta\to\pm\infty$); Firth conserta |
| **logit / OR** | log-odds; $\exp(\beta)$ = odds ratio (efeito do burden) |
| **IVW meta-análise** | combina $(\beta_k,SE_k)$ de coortes ponderando por **precisão** ($1/SE^2$) |
| **SKAT** | teste variance-component; tolera variantes de direções opostas (≠ burden) |
| **REGENIE / SAIGE** | ferramentas modernas: caso-controle desbalanceado + parentesco + escala |

---

### Pendências / próximos tópicos a estudar
- [x] **Modelo estatístico** (regressão logística, burden aditivo, idade/idade²/sexo/10 PCs, meta-análise IVW EUR+AFR, Firth) → **seção 11** ✅ 2026-06-04.
- [ ] **QQ-plot / λ (lambda) e inflação:** como o limiar `p<10⁻⁶` saiu do ponto em que o QQ observado "descola" do esperado (Ext. Data Fig. 3; λ∆95 = 1,558 global / 1,09 AFR / 1,251 EUR). É o controle de qualidade da estatística.
- [ ] **Conferir no runbook do Daniel:** corte de REVEL e papel (burden primário?), piso de portadores, cortes de MAF (externo gnomAD + interno coorte), build do genoma.
- [ ] **Decisão:** re-anotar com VEP+LOFTEE/gnomAD v4 (GRCh38) vs. manter ANNOVAR/gnomAD v2 (fidelidade ao Daniel).

*Notas criadas em 2026-06-04 a partir de sessão de estudo dos Métodos. Pareia com `paper_summary_park2021.md` (resumo formal) e `paper_summary_hui2023.md` (paper de HL que a replicação segue).*
