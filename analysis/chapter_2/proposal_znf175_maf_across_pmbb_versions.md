# Proposta — Comparação do MAF das variantes raras de ZNF175 entre versões do PMBB (v1→v4)

**Autor:** Andre Rico
**Data:** 2026-06-04
**Status:** Proposta (pré-execução) — para discussão com Molly / Nikki
**Contexto do projeto:** Capítulo 2. Resolução do enigma de perda de sinal ZNF175–tinnitus (11K → 45K).

---

## 1. Motivação — o enigma 11K → 45K

A associação **ZNF175 → tinnitus** foi significativa em ~11K exomas (Park et al., *Nat Med* 2021, pLOF) e **sumiu em ~45K** (replicação do Daniel — negativa). A pergunta central, conforme a reunião de 2026-06-03 (Doug/Molly):

> **A perda de sinal é artefato técnico ou null real?**

A diretriz da Molly é explícita: **descartar a explicação técnica antes de invocar biologia**. Andre já reproduziu o resultado dos 45K (negativo). Falta o outro lado do espelho.

Esta proposta é o **primeiro passo concreto e barato** dessa investigação (Action Item #1 da reunião): antes de re-rodar qualquer teste de associação, verificar se o **substrato de variantes do ZNF175 é estável** entre as versões do PMBB.

> **Princípio:** se o conjunto de variantes raras do ZNF175 (ou seus MAFs) não é estável entre os freezes, re-rodar a associação é inútil. A estabilidade do substrato é **pré-requisito**, não consequência.

---

## 2. A proposta em uma frase

Levantar **todas as variantes raras do ZNF175** e comparar sua **presença e MAF** entre os freezes do PMBB (v1=11K → v2 → v3 → v4 atual), de modo a **isolar o componente técnico** (pipeline/calling/build) das mudanças não-técnicas (crescimento da coorte e composição de ancestralidade).

---

## 3. O que isto diagnostica

Uma comparação variante-a-variante entre freezes responde perguntas que a re-execução do burden **não** responde:

1. **O substrato é estável?** Com o crescimento + re-processamento do PMBB (provável mudança de build GRCh37→GRCh38, joint-calling diferente), o conjunto de variantes pLOF do ZNF175 pode mudar — variantes do v1 somem no v4 (filtradas, re-genotipadas como ref, dropadas no QC) ou aparecem novas.
2. **As "raras" continuaram raras?** Uma variante a MAF 0,0008 no v1 pode virar 0,0012 no v4 e **cruzar o limiar de 0,1%** → sai do burden. Mecanismo concreto de mudança do conjunto.
3. **Os ~8 casos que dirigiam o sinal** seguem presentes e raros? (O sinal era *small-case* → território de **winner's curse**.)
4. **Multialélico / stop-loss:** o repo já tem `data/PMBB_Exome/addBack_multiallelic_stoploss/`, sinalizando que normalização de multialélicos/stop-loss do ZNF175 **já foi um problema conhecido**. Uma comparação de MAF entre freezes expõe exatamente esse tipo de artefato. **Prior alto de achar algo.**

---

## 4. ⚠️ O ajuste crítico de desenho — não confundir 3 coisas

O MAF de uma variante mudar entre v1 e v4 pode ter **três causas**, e só **uma** é o alvo:

| Fonte da mudança no MAF | É o alvo? |
|---|---|
| **(1) Pipeline / calling / build** (re-chamada, normalização, liftover) | ✅ **SIM** — o "artefato técnico" |
| **(2) Crescimento da coorte** (mais gente → estimativa refina; rara regride à freq. real) | ❌ mecânica do winner's curse, não artefato |
| **(3) Composição de ancestralidade** (mix EUR/AFR diferente por freeze) | ❌ confundidor puro (MAF é ancestria-dependente) |

### Como isolar o componente técnico (1):

- **Comparar MAF DENTRO de estrato de ancestralidade** (MAF_EUR e MAF_AFR separados), **nunca** o global → elimina (3).
- **Se os freezes forem aninhados** (v1 ⊂ v4 em pessoas): fazer **concordância de genótipo por-amostra** nos indivíduos compartilhados. Como a pessoa é a mesma, qualquer call diferente entre v1 e v4 é **100% técnico** → isola (1), removendo (2) e (3). **É o teste mais forte.**
- **Rastrear variantes por HGVS / rsID normalizado**, não por coordenada bruta (build/normalização mudam coordenadas; é onde o multialélico/stop-loss aparece).

---

## 5. A pegadinha: "incremental" ≠ "genótipo estável"

Hipótese de trabalho: as versões do PMBB são **incrementais** (recrutamento cumulativo, v_{n+1} ⊇ v_n). **Isto deve ser verificado, não assumido** (passo zero).

Ponto-chave (a favor da proposta): mesmo que aninhado em pessoas, o **genótipo da mesma pessoa pode mudar** entre versões por motivos técnicos:

1. **Joint calling** — adicionar amostras e re-rodar a chamada conjunta **recomputa** os genótipos dos indivíduos antigos. A mesma pessoa pode ter call diferente no v1 e no v4 sem ter mudado nada.
2. **Build / normalização** (GRCh37→GRCh38, multialélico/stop-loss).

→ Por isso a incrementalidade é um **presente**: torna a concordância por-amostra um teste **cirúrgico** de artefato técnico.

**Ressalva:** recrutamento incremental pode **não** virar aninhamento estrito na análise — amostras do v1 podem ser **dropadas no QC** de versões posteriores. Por isso o passo zero é intersectar as listas de IDs.

---

## 6. Passo zero (rápido e decisivo)

1. **Confirmar a relação entre freezes:** intersectar os IDs de amostra de v1/v2/v3/v4 → são aninhados? Qual o N compartilhado? *(Decide: concordância por-amostra vs. só MAF agregado.)*
2. **Amarrar o mapeamento** versão ↔ tamanho ↔ build (v1=11K/GRCh37? v4=?/GRCh38?).
3. **Localizar os assets existentes:** `data/PMBB_Exome/ZNF175/` e `data/PMBB_Exome/addBack_multiallelic_stoploss/` — o que já está levantado e em que build.

---

## 7. Entregável

Tabela-mestra das variantes do ZNF175:

| Variante (HGVS) | Build | v1 | v2 | v3 | v4 | … |
|---|---|---|---|---|---|---|
| presente? (Y/N) | | | | | | |
| MAF_EUR | | | | | | |
| MAF_AFR | | | | | | |
| N_portadores | | | | | | |
| pLOF? (sim/não) | | | | | | |

Mais, se os freezes forem aninhados:
- **Matriz de concordância de genótipo por-amostra** (mesmos IDs, v1 vs v4) para as variantes-chave.
- **Destaque dos ~8 casos** que dirigiam o sinal original.

---

## 8. Árvore de resultados (qualquer desfecho é informativo)

- **Variantes / MAF instáveis entre freezes** → encontramos (parte da) **explicação técnica** que a Molly pediu. 🎯 (provável foco: multialélico/stop-loss, drift de MAF cruzando 0,1%, drops de QC.)
- **Rock-stable em todos os freezes** → substrato **descartado** como culpado → a perda pende para **winner's curse / null real**, e restringimos o problema ao pipeline de *associação* (ex.: detalhe da ata de que as replicações viraram **missense**, não pLOF).

---

## 9. Encaixe na reunião (2026-06-03)

- É o **caminho mais barato dentro do Action Item #1** (reproduzir o pipeline 11K do Park p/ apples-to-apples), e **precede** a re-execução do burden.
- **Não depende** da audiometria (bloqueada pelo matching audiograma↔ID, Action #7) nem de localizar o Joe Park primeiro.
- Usa assets que já existem no repo (`ZNF175/`, `addBack_multiallelic_stoploss/`).
- Bom **quick-win diagnóstico** para levar à Molly/Nikki.

---

## 10. Escopo — o que esta proposta NÃO é

- **Não** re-testa a associação ZNF175–tinnitus (é forense de variante / QC; vem antes).
- **Não** usa audiometria (fase quantitativa, em v4, depende do matching).
- **Não** é o estudo novo exome-wide de HL em v4 (raia da Nikki/Elena, métodos SKAT).
- **Não** assume nesting — verifica.

---

## 11. Questões em aberto / a confirmar

- [ ] Quais freezes (v1–v4) estão de fato acessíveis no LPC? Tamanhos e builds?
- [ ] As versões são aninhadas em IDs de amostra? (passo zero)
- [ ] Quais variantes do ZNF175 já estão levantadas em `data/PMBB_Exome/ZNF175/`? Em que build?
- [ ] Definição do conjunto: só pLOF, ou pLOF + missense REVEL? (Para o enigma, o relevante é o **pLOF** que dirigiu o discovery do Park.)
- [ ] MAF de referência: interno (coorte) e/ou gnomAD (qual versão/ancestralidade)?

---

## Referências do projeto
- `docs/papers/paper_summary_park2021.md` — resumo do Park 2021 (fonte do ZNF175→tinnitus).
- `docs/papers/andre_notes_park2021.md` — notas de método (pLOF, MAF, burden, modelo estatístico).
- `docs/meetings/2026-06-03_research-vision-and-elena-onboarding.md` — ata (enigma 11K→45K, Action Items).
- `data/PMBB_Exome/ZNF175/`, `data/PMBB_Exome/addBack_multiallelic_stoploss/` — assets do Daniel.
- `analysis/daniel/` — runbook de replicação do Hui (pipeline irmão).
