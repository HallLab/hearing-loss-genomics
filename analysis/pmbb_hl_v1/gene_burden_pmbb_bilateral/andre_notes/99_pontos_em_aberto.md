# Pontos em aberto — o que é bloqueio e o que é lembrete

**Autor:** Andre Rico · **Data:** 2026-10-08
**Cópia de trabalho pessoal, em português.** Não é página de Confluence.

Os quatro documentos de fase terminam com "Em aberto", e eles misturam coisas de natureza muito
diferente — um item que mudaria a conclusão está na mesma lista que um parâmetro que poderia ser
outro. Esta página separa, para a pergunta *"o que precisa ser resolvido?"* ter resposta.

**A resposta curta: nada bloqueia a conclusão atual.** Os 14 itens se distribuem em quatro grupos, e
nenhum está no grupo que mudaria o que reportamos.

---

## A · Bloqueiam a conclusão

**Nenhum.**

A conclusão desta análise é: *nenhum gene atinge significância, e o nulo é falta de poder, não falta
de sinal* — esta última parte medida por dois controles. Nenhum dos itens abaixo mudaria essa frase.

Vale afirmar explicitamente porque o inverso seria grave: se houvesse um item aqui, o resultado não
deveria ser apresentado.

---

## B · Escolhas de parâmetro que mudariam o ranking

São decisões defensáveis que poderiam ter sido outras. **Não mudam a conclusão** — nada passa o
limiar sob nenhuma delas —, mas mudariam *quais genes ficam no topo*, e topo é o que alimenta
follow-up.

| | hoje | alternativa | custo de testar |
|---|---|---|---|
| **REVEL** | 0,5 | 0,6 | refazer a `pDM` e rodar · ~1 dia |
| **SpliceAI** | 0,2 | 0,5 | refazer a `pLOF` e rodar · ~1 dia |
| **Súbita idiopática** (`SO_396.5`) | fora | dentro, +203 casos | refazer a Fase 1 e rodar · ~1 dia |

**O que eu faria:** nenhuma delas agora. Rodar análise de sensibilidade para um resultado nulo
produz mais resultados nulos. Elas passam a valer **se** a coorte crescer — e aí viram as primeiras
coisas a variar.

A terceira é diferente das outras duas: não é escolha estatística, é **clínica**. Se perda súbita
idiopática pertence a um estudo genético de perda relacionada à idade é pergunta para o Doug, não
para mim.

---

## C · Lacunas de cobertura, com tamanho medido

Aqui o valor é saber **o tamanho**, mesmo sem resolver.

### Sem cromossomo X — custa 6% dos genes conhecidos

Medido: dos **100 genes ClinGen Definitive/Strong** de surdez, **6 estão no X**:

```
AIFM1   COL4A5   POU3F4   PRPS1   SMPX   TIMM8A
```

Então nossos testes cobrem 94 de 100. É ponto cego real num estudo de surdez — existe surdez ligada
ao X — mas é 6%, não um terço.

**Resolver exige:** máscaras do X (o VEP da release cobre), e o SAIGE tratar hemizigosidade em homens.
Não é trivial, e é a lacuna mais defensável de fechar se alguém perguntar o que falta.

### Cobertura do AlphaMissense — não medido

Para genes sem escore do AlphaMissense, só o REVEL decide se uma variante é `pDM`. Quantos genes
ficam nessa situação **nunca foi medido**. É barato de medir e ninguém mediu.

---

## D · Limitações que não se resolvem aqui

Não são dívida nossa. São propriedades da coorte.

**AFR com 490 casos, e `tau₂` = 0.** Não há componente poligênico detectável nesse tamanho. Rodou
porque não rodar seria viés de reporte. Não se resolve com análise — se resolve com mais gente.

**Os 15 `UNKNOWN` e os grupos pequenos.** EAS, SAS, AMR e os não classificados existem só na coorte
`combined`, sem N para braço próprio. Os PCs do combined é que descontam a mistura.

**Os 438 GB de genótipos fora da pasta.** Único caminho que aponta para fora. Copiar não tornaria
nada mais seguro — é conversão mecânica da release, não codifica decisão. Registrado no
[`PROVENANCE.md`](../PROVENANCE.md).

---

## E · Pendências técnicas que medi e não importam

Dois itens que eu levantei como preocupação e que **a medição desarmou**. Mantidos porque a medição
é o que os desarma — sem ela voltariam a parecer problema.

### O GRM da coorte combinada é magro — e a calibração está boa

38.833 marcadores contra 114.694 do AFR, porque `--hwe 1e-6` aplicado **entre** ancestralidades
remove 176 mil variantes onde o mesmo filtro dentro do EUR remove 19 mil. Desvio de Hardy-Weinberg
numa amostra estruturada é esperado; o filtro faz o trabalho errado ali.

**Mas o λ desta análise está bom:**

```
combined   pDM 0,944   pLOF 0,964   pLOF;pDM 0,952
EUR        pDM 0,999   pLOF 0,958   pLOF;pDM 0,984
AFR        pDM 0,843   pLOF 0,664   pLOF;pDM 0,972
```

Perto de 1, levemente conservador — que é o que um GRM magro produziria se estivesse atrapalhando, e
não está o bastante para importar. O AFR/pLOF em 0,664 é o único fora da faixa, e é a coorte com 490
casos: tamanho, não GRM.

**Conclusão:** vale registrar, não vale refazer. Reconstruir exigiria decidir como tratar HWE em
amostra multi-ancestral, que é pergunta para o time e não desbloqueia nada.

### O perfil estranho da PCA do AFR

O PC1 domina muito mais do que o esperado. Outliers e normalização foram descartados como causa, e
**não aparece no EUR**. Registrado em [`phase_3/results/PCA_NOTES.md`](../phase_3/results/PCA_NOTES.md).

Sem explicação — mas os PCs funcionam: o λ do AFR está em faixa aceitável, e a escolha de 3
componentes é sustentada por dois caminhos independentes. Curiosidade, não defeito.

---

## O que de fato muda o estudo

Nada desta lista. As três coisas que mudariam estão fora dela:

| | o que muda | estado |
|---|---|---|
| **Audiogramas** | o fenótipo deixa de ser código de diagnóstico e vira medida quantitativa | falta a ponte REDCap ↔ PMBB |
| **Coorte maior** | All of Us, UK Biobank — é o que o resultado nulo pede | não iniciado |
| **Replicação externa** | o terceiro sentido de "replicação" que o Doug nomeou | não iniciado |

As duas primeiras atacam a causa raiz — **3.164 casos é pouco, e código de diagnóstico é substituto
pobre de audiograma**. Todo o resto desta página é ajuste de margem.
