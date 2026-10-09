# Fase 2 — quais variantes contam · análise

Braço: **ExWAS de variante comum, fenótipo bilateral neurossensorial**.
Rodado em 2026-10-09, job LSF `49576769`, 66 de 66 tarefas.

> **Nota sobre este documento.** Os conceitos desta fase — por que MAF depende de quem
> está na coorte, por que frequência de ALT não é MAF, por que o censo existe — estão
> explicados por extenso em
> `../exwas_elena_replication/andre_notes/phase_02_analise_pt.md`. Não repito aqui de
> propósito: explicação duplicada em dois arquivos é explicação que desincroniza quando
> um dos dois é corrigido. Aqui ficam os **números deste braço** e o que difere.

---

## O que muda em relação ao braço amplo

Uma coisa só: **quem está na coorte**.

| | amplo | bilateral |
|---|---:|---:|
| casos | 6.752 | 3.164 |
| controles | 50.746 | 50.746 |
| total | 57.498 | 53.910 |

Os controles são **as mesmas pessoas**. Os casos bilaterais são um subconjunto estrito
dos casos amplos. A restrição remove 3.588 casos e mais nada.

Como MAF é propriedade do grupo de pessoas, o censo tem que rodar de novo — mas como
as coortes diferem em 6% do total, o resultado é quase o mesmo.

---

## O censo

De 13.223.743 variantes no arquivo:

| coorte | amostras | censo plink2 | o que o SAIGE testou | barra secundária |
|---|---:|---:|---:|---:|
| combined | 53.910 | 70.374 | 70.350 | 7,107e-07 |
| EUR | 40.143 | 74.585 | 74.566 | 6,705e-07 |
| AFR | 10.539 | 118.083 | 118.049 | 4,236e-07 |

Comparando com o braço amplo (70.366 / 74.536 / 117.953): diferenças de dezenas de
variantes em dezenas de milhares. Esperado — trocar 3.588 casos por nada move pouquíssimas
variantes através do limiar de 1%.

A diferença entre a coluna "censo plink2" e "o que o SAIGE testou" é o mesmo mecanismo
documentado no braço amplo: o plink2 divide pelos genótipos **observados**, o SAIGE pela
coorte **inteira**. O censo é teto, não igualdade.

**Zero variantes com frequência de ALT acima de 95%** em qualquer coorte, igual ao outro
braço.

---

## As duas barras

- **Primária `5e-8`** — convenção, e a mesma do outro braço para os dois serem comparáveis
- **Secundária `0,05/n`** — Bonferroni sobre os testes realmente feitos

As duas fixadas antes do step 2 rodar.

Repare que a barra **não muda** por termos metade dos casos. Barra é sobre quantos testes
foram feitos, não sobre quanto poder temos. O que muda com menos casos é a chance de
chegar lá — e isso está na fase 5.

---

## Arquivos

| | |
|---|---|
| `phase_2/data/{combined,EUR,AFR}.keep` | listas de amostras, conferidas contra o `.fam` |
| `phase_2/results/census/<coorte>.chr<N>.afreq.alt1.bins` | os 66 histogramas |
| `phase_2/results/census_summary.json` | o agregado |
