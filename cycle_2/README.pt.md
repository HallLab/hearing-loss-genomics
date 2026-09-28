# Ciclo 2 — Variação Rara em Genes de Perda Auditiva Mendeliana e Perda Auditiva de Início Adulto

**Status:** RASCUNHO de charter — para revisão de Molly Hall, Douglas Epstein e Nikki Palmiero antes de rodar qualquer análise
**Autor:** Andre Rico · **Aberto em:** 2026-08-26
**Coorte:** PMBB Release 2026-4.0 (v4)

> **Nota:** esta é a versão em português, para uso interno do Andre. A versão canônica, que circula com o time,
> é [`README.md`](README.md). Se as duas divergirem, vale o inglês.

---

## 0. Por que um ciclo novo

Tudo que está fora desta pasta é o **Ciclo 1** (`analysis/`, `data/`, `docs/`, `results/`, `scripts/`). O Ciclo 1
foi organizado em torno de um *gene* — ZNF175 — e cumpriu o seu papel: respondeu a pergunta que recebeu, e a
resposta foi negativa. O Ciclo 2 é organizado em torno de uma *pergunta*, para que o projeto continue tendo
sentido qualquer que seja o resultado.

O Ciclo 1 não é movido nem arquivado — segue vivo (os pipelines de v4 da Elena, dois blockers de dados em
aberto). O Ciclo 2 fica ao lado dele e reaproveita sua infraestrutura.

### O que o Ciclo 1 estabeleceu (fechado — não reabrir)

| Achado | Onde |
|---|---|
| Hui et al. 2023 replicado em PMBB v2; os 6 genes do paper recuperados a FDR<0,05 | [`results/chapter1_paper_replication/`](../results/chapter1_paper_replication/) |
| O sinal ZNF175→tinnitus em 11K reproduz (OR≈14,6) e decai em 44K | [`analysis/chapter_2/findings_znf175_11k_vs_44k.md`](../analysis/chapter_2/findings_znf175_11k_vs_44k.md) |
| O decaimento é **winner's curse**, não artefato de pipeline — os mesmos 4 indivíduos dirigem o sinal nos dois freezes | [`analysis/chapter_2_v2/preconclusion_znf175_4_vs_8.md`](../analysis/chapter_2_v2/preconclusion_znf175_4_vs_8.md) |
| A discrepância "8 vs 4" carriers resolve como `8 = 6 + 2` (fenótipo mais amplo + carriers sem linkage) | idem |
| **ZNF175 é nulo em PMBB v4** em todas as masks × MAF × ancestralidade, com beta negativo | [`analysis/elena/HL_only_rarevariant/`](../analysis/elena/HL_only_rarevariant/) |

**Consequência:** ZNF175 deixa de ser um assunto. No Ciclo 2 ele é uma linha de uma tabela de resultados, nada além.

---

## 1. Pergunta de pesquisa

### Primária — Q1

> **Variação rara danosa nos genes de perda auditiva mendeliana conhecida contribui para perda auditiva de
> início adulto numa biobank hospitalar — e, se contribui, em quais genes e com que tamanho de efeito?**

### Secundária — Q2 (condicional a Q1 ter definido um conjunto de portadores)

> **Entre os portadores dessas variantes, o que separa quem desenvolve perda auditiva de quem não desenvolve?**
> (penetrância; idade, sexo, ancestralidade; um segundo hit em outro gene do conjunto; `CX3CR1` como o único
> modificador nomeado *a priori*, vindo do trabalho em camundongo)

Q1 e Q2 são aninhadas, não paralelas: **Q1 constrói o conjunto de portadores, Q2 pergunta o que acontece dentro
dele.** Q1 vale a pena independentemente do resultado; Q2 só é interpretável depois que Q1 existe.

### Por que essa pergunta e não outra

- **Tem um prior real.** ~100 genes de HL congênita não-sindrômica bem caracterizados. Um conjunto
  pré-especificado de ~100–200 genes substitui a barra exome-wide (~2,5×10⁻⁶) por FDR sobre ~10² testes — um
  ganho grande de poder, de graça.
- **A resposta é genuinamente aberta.** Esses genes são adjudicados para HL *congênita, majoritariamente
  recessiva e majoritariamente severa*. Se o estado de portador heterozigoto eleva risco de HL de *início adulto*
  é uma questão não testada, não uma conclusão dada. **Não** esperamos sinal forte — se esperássemos, o estudo não
  valeria a pena. O desfecho esperado é nulo.
- **O negativo é publicável.** "Portadores de variantes raras danosas em genes de HL mendeliana não estão
  enriquecidos para HL de início adulto no EHR" contradiz diretamente a premissa translacional enunciada no
  kickoff (*ligar HL adulta a genes de HL congênita bem caracterizados*). Isso é informação, não fracasso.
- **Sobrevive ao próprio nulo.** Diferente de um projeto de gene único, nenhum desfecho nos deixa sem nada a reportar.

---

## 2. O que temos (contado, não presumido)

### Coorte de phecodes — PMBB v4

| Fenótipo | N com exoma | casos | prevalência |
|---|---|---|---|
| Hearing impairment | 57.632 | **6.712** | 11,7% |
| Tinnitus | 53.096 | 2.732 | 5,1% |

*(contado a partir dos arquivos de covariáveis SAIGE da Elena, `analysis/elena/HL_only_rarevariant/new_SAIGE_covariates/`)*

Tamanho de amostra deixou de ser a restrição. **A qualidade do fenótipo é.** O phecode `hearing impairment`
junta perda condutiva, súbita, induzida por ruído e relacionada à idade num rótulo só, e a maior parte da HL
adulta é poligênica e etária. É por isso que o audiograma importa — ver abaixo.

### Coorte audiométrica — o ativo subutilizado

A base de audiometria do Brant (`audbase`) já está em disco desde o Ciclo 1, junto com uma **linkagem ao PMBB
que já foi construída uma vez, em fevereiro de 2021**:

| Arquivo | Conteúdo |
|---|---|
| [`data/PMBB_Exome/brant audbase 1.7.21.TXT.gz`](../data/PMBB_Exome/) | audbase cru, 100.470 sujeitos / 216.542 registros — **contém PHI** |
| [`data/PMBB_Exome/audbase_feb252021/RGC21_45k_aud_1.csv.gz`](../data/PMBB_Exome/audbase_feb252021/) | fenótipos derivados: PTA aéreo/ósseo por orelha, `PTA`, `Bilateral_HL`, `Worse_ear`, `Degree_HL` (0–4), `BL_SNHL` |
| [`data/PMBB_Exome/audbase_feb252021/degree_HL_aud.txt.gz`](../data/PMBB_Exome/audbase_feb252021/) | 3.328 `PMBB_ID → Degree_HL_Aud`, desidentificado |

**O número que importa:** esse arquivo tem **45.012 sujeitos de audiometria identificados por MRN, dos quais
apenas 3.328 (7,4%) carregam um PMBB_ID.** Os ~41.700 restantes estão sem linkage. O match foi feito em fev/2021
por MRN + data de nascimento contra o PMBB daquela época; o PMBB cresceu bastante desde então. Os dois lados do
join cresceram.

Ou seja, o item registrado como blocker na ata de 2026-07-01 — *"audiogram data not yet matched to PMBB IDs"* —
é, mais precisamente: **existe uma linkagem, está cinco anos desatualizada, e cobre 7,4% da coorte
audiométrica.** Estamos atualizando e estendendo um exemplo que já funcionou, não construindo do zero.

> ⚠ **Governança de dados.** O audbase cru e o `RGC21_45k_aud_1.csv.gz` contêm identificadores (nome, endereço,
> telefone, e-mail, DOB, MRN). Qualquer re-linkagem é ação de honest broker / IRB sob o projeto *Audiometric
> Phenotyping of PMBB Enrollees* — não é algo para rodarmos por conta própria. A ata de 2026-07-01 registra a
> Elena vinculada ao **IRB errado** (Callback, e não Audiometric Phenotyping); até a correção da Nikki sair,
> esses arquivos estão fora do escopo dela. Ratchet (lab da Ritchie) é o responsável designado pelo refresh
> para v4 — coordenar, não duplicar.

---

## 3. Desenho

Dois níveis, cada um fazendo o trabalho para o qual serve:

| Nível | Coorte | Fenótipo | Papel |
|---|---|---|---|
| **A — Descoberta** | 57.632 exomas | phecode HL (binário), n=6.712 | tem poder; rótulo ruidoso |
| **B — Validação** | ~3.328 exoma + audiograma (cresce com o refresh) | `Degree_HL` / PTA (quantitativo), `BL_SNHL` | rótulo limpo; sem poder sozinho |

O Nível B **não substitui** o Nível A — 3,3K é pequeno demais para descoberta de variante rara. Ele valida e
refina: um gene que sobrevive ao Nível A e mostra deslocamento dose-consistente de PTA no Nível B é uma
afirmação qualitativamente diferente de um que só passa num FDR de phecode. É também o desenho que faz o
refresh da linkagem valer cientificamente, e não apenas administrativamente.

### Escolhas analíticas pré-especificadas

Herdadas das decisões de 2026-07-01, mais as lacunas que aquela reunião deixou em aberto:

| Escolha | Valor | Origem |
|---|---|---|
| Fenótipo | HL **e/ou** tinnitus, combinado | decidido em 2026-07-01 — *nunca construído; o Ciclo 2 constrói* |
| Conjunto de genes (primário) | **ClinGen HL GCEP (40007), Definitive+Strong — n=100 genes**, snapshot 2026-08-26 | **provisório** — Andre 2026-08-26; levar para ratificação |
| Conjunto de genes (sensibilidade) | mesmo painel, +Moderate — n=120 genes | provisório, mesma decisão |
| Masks | pLOF; pLOF+AlphaMissense; pLOF+REVEL; pLOF+AM+REVEL — **mantidas separadas** | decidido em 2026-07-01; o Ciclo 1 colapsou tudo num único `pDM` |
| Limiar REVEL | 0,5 primário, 0,6 como sensibilidade | decidido em 2026-07-01 (não resolvido) |
| MAF | 0,01 / 0,001 / 0,0001 | decidido em 2026-07-01 |
| PCs | 5–6 (nunca 20) | decidido em 2026-07-01 |
| Correção múltipla | **FDR dentro do grupo de MAF**, não entre modelos | decidido em 2026-07-01 — *nunca aplicado a nenhum resultado* |
| Método | SAIGE-GENE+ (burden + SKAT + Cauchy) | pipeline v4 da Elena, reaproveitado |
| Controles técnicos | pares gene–fenótipo fora da audição com efeito grande e conhecido no PMBB (`BRCA1`/câncer de mama, `TTN`/cardiomiopatia, `CFTR`/FC — o conjunto que o Park 2021 usou para validar), mais λ_GC e QQ na camada exome-wide | novo no Ciclo 2 |

> **Não são controle:** `GJB2`, `SLC26A4`, `MYO7A` e o resto do conjunto primário são a *hipótese sob teste*, não
> uma verificação de validade dela. Tratá-los como controle positivo seria assumir a resposta. A validade do
> pipeline se estabelece em pares gene–fenótipo fora da audição.

### O que conta como resposta

- **Q1 positivo:** ≥1 gene do conjunto primário a FDR<0,05 dentro de um grupo de MAF, com os controles técnicos
  se comportando e λ_GC em faixa.
- **Q1 negativo:** controles técnicos se comportam, calibração limpa, nenhum gene sobrevive → variação em estado de
  portador nos genes de HL mendeliana não contribui de forma detectável para HL adulta definida por EHR. **Este é o
  desfecho esperado**, e é reportável: delimita a premissa translacional enunciada no kickoff.
- **Q1 não informativo:** controles técnicos falham ou a calibração está ruim → o achado é sobre o fenótipo ou o
  pipeline, e o Nível B vira prioridade.

---

## 4. Guarda de escopo — o que o Ciclo 2 *não* é

- Não é um projeto sobre ZNF175.
- Não é teste formal de interação gene×gene (reconhecido sem poder desde o kickoff; Q2 fica descritivo).
- Não é Menière (poucos casos; coordenar com Bogdan/Ian antes de entrar).
- Não é UK Biobank (acesso congelado desde maio de 2026).
- Não é um paper metodológico. O trabalho de winner's curse do Ciclo 1 é material de discussão, não a
  afirmação primária.

---

## 5. Dependências em aberto

| # | Item | Responsável | Bloqueia |
|---|---|---|---|
| 1 | Qual lista adjudicada de genes de HL é a autoritativa | Doug / Andre | ~~bloqueia Q1~~ — **resolvido provisoriamente em 2026-08-26** (ver §5.1 no README.md); na pauta da próxima reunião |
| 2 | Re-linkagem audbase ↔ PMBB v4 | Ratchet (lab da Ritchie); escalado pelo Doug | Nível B |
| 3 | Correção de IRB da Elena (Audiometric Phenotyping) | Nikki | Elena tocar no Nível B |
| 4 | Crosswalk mestre GENO_ID ↔ PT_ID (7 carriers do Ciclo 1) | Nikki / curadores do PMBB | resíduo do Ciclo 1 apenas |
| 5 | Limiar final de REVEL (0,5 vs 0,6) | Molly / Nikki | construção das masks |

---

## 6. Próximos passos

1. Revisar este charter com Molly, Doug e Nikki — **antes de rodar qualquer coisa** (a cadência definida no kickoff).
2. ~~Resolver a dependência #1 (a lista de genes).~~ Resolvido provisoriamente — ver §5.1 no [`README.md`](README.md).
   Levar à próxima reunião; se a sala discordar, basta re-rodar o fetch com outra camada — o downstream não muda.
3. Construir o fenótipo combinado HL-e/ou-tinnitus para v4 — a única decisão de 2026-07-01 que nunca foi executada.
4. Preparar a especificação técnica do refresh audbase↔PMBB (chaves de join, rendimento esperado, QC) para que
   Doug e Ratchet recebam um pedido concreto em vez de um pedido em aberto.
