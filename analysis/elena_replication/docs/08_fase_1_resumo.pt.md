# Fase 1 — quem entra no estudo · resumo

**Cópia de trabalho pessoal.** Detalhamento completo: [`phase_1/results/FINDINGS.md`](../phase_1/results/FINDINGS.md) ·
Registro com toda a procedência: [`02_phase_1_reference.md`](02_phase_1_reference.md) ·
Explicação em linguagem leiga: [`03_step_1_explained.pt.md`](03_step_1_explained.pt.md)

---

## Por que esta fase existe

Antes de perguntar *qual gene causa perda auditiva*, duas coisas precisam estar resolvidas: **quem
tem perda auditiva** e **quem está no estudo**. Tudo que vem depois herda essa resposta, e erro aqui
não tem conserto adiante.

## O que fizemos

Reconstruímos casos e controles do zero a partir do release, sem olhar o código dela, e seguimos a
coorte até o arquivo que o SAIGE de fato consumiu.

## O que encontramos

**O fenótipo está certo.** 70.925 de 70.925, pessoa a pessoa. Regra de 2, exclusão do meio ambíguo e
o tratamento do tinnitus na tabela `observation` — tudo reproduz.

**Mas a coorte analisada não é a que o fenótipo sustenta:**

```
57.507   analisáveis
  -431   nunca entraram (40 casos)  ← filtro contra o .fam do IMPUTADO
  +556   entraram indevidamente     ← fenótipo sem a tabela observation
57.632   o que rodou
```

Três defeitos independentes:

| | |
|---|---|
| **filtro errado** | cortou pelo `.fam` do array; a análise usa exoma do início ao fim |
| **arquivo sobrescrito** | o fenótipo foi corrigido 31/jul 20:57; o SAIGE rodou 3/ago nas covariáveis de 20:21 |
| **exclusão concentrada** | **15,6% dos participantes do Leste Asiático** contra 0,73% da coorte |

Todos empurram para o nulo. O risco é **achado perdido**, não achado falso.

## Como seguimos

**Usamos a nossa coorte**, sem filtro de dado imputado e sem adições externas:

| | reprodução | corrigido |
|---|---|---|
| combined | 57.632 · 6.712 casos | **57.507 · 6.752 casos** |
| EUR | 43.016 | **42.786** |
| AFR | 11.387 | **11.334** |

Arquivos em `phase_1/results/cohorts/`. O braço de reprodução existe para provar que os números dela
reproduzem; o corrigido é o que seguimos.

Corrigir **restaura 165 participantes do Leste Asiático** — aquele estrato cresce 17,6% enquanto
todos os outros mudam meio por cento.

⚠️ O arquivo de covariáveis que o SAIGE consumiu é o **único registro do que foi analisado**, porque o
fenótipo de origem foi sobrescrito. Se for apagado, o braço de reprodução acaba.

## Em aberto

- por que o filtro usou o `.fam` do imputado — Nikki e Elena
- por que o ramo de perda auditiva não foi refeito, se o de tinnitus foi
- 4 controles legítimos que somem sem explicação
- tinnitus como fenótipo — **fora de escopo**, verificamos só perda auditiva
