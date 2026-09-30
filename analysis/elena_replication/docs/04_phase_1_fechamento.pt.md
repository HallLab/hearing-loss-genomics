# O que fechamos na Fase 1

**Autor:** Andre Rico · **Data:** 2026-09-30
**Cópia de trabalho pessoal, em português.** Não é página de Confluence. O registro formal está em
[`phase_1/results/FINDINGS.md`](../phase_1/results/FINDINGS.md) e a referência completa em
[`02_phase_1_reference.md`](02_phase_1_reference.md).

---

## A pergunta

Antes de perguntar *qual gene causa perda auditiva*, é preciso saber **quem tem perda auditiva** e
**quem está no estudo**. A Fase 1 verifica essas duas coisas. Tudo que vem depois herda a resposta —
e um erro aqui não tem conserto lá na frente.

---

## O que está certo

**A regra que decide quem é caso reproduz exatamente.**

Reconstruí casos e controles do zero, a partir dos dados brutos do hospital, sem olhar o código da
Nikki. Deu igual: **70.925 de 70.925 pessoas, uma a uma.** Não é só o total bater — são as mesmas
pessoas.

Isso cobre a regra de 2 (diagnóstico em duas datas), a exclusão de quem tem doença de ouvido que não
é perda auditiva, e o tratamento do zumbido, que na versão 4 do PMBB mudou de tabela e é a armadilha
mais fácil de cair.

**A estratificação por ancestralidade também confere.** Os arquivos EUR e AFR contêm só quem o
release classifica como EUR e AFR, e o `combined` é exatamente a soma com os grupos menores.

Vale registrar: em um dos checks **eu estava errado e o pipeline estava certo**. Minha primeira
reconstrução discordou em 558 pessoas porque eu tinha esquecido a tabela nova do zumbido. A
conferência corre nos dois sentidos.

---

## O que está errado

### 1. Foram cortadas 431 pessoas que podiam ser analisadas

O estudo roda sobre **exoma**. Mas a lista de participantes foi filtrada por quem tinha dado de
**array** — um exame diferente, que a análise não usa em lugar nenhum.

Resultado: 431 pessoas com exoma completo ficaram de fora, **40 delas com perda auditiva**.

Não é divergência de método. É inconsistência interna: o pipeline barra com base num dado que
nenhuma etapa dele consome.

### 2. Entraram 556 pessoas que a própria regra exclui

Pessoas com doença de ouvido que **não** é perda auditiva — nem caso, nem controle. A regra manda
separá-las. Elas entraram como controles.

Importante: **nenhuma delas tem qualquer diagnóstico de perda auditiva.** Não são casos disfarçados.
A regra as separa por precaução, não por saber que estão doentes.

### 3. O arquivo analisado não existe mais

Este é o mais incômodo. A cronologia:

```
31/jul 20:21   covariáveis de perda auditiva construídas   (57.632 pessoas)
31/jul 20:57   todo o fenótipo é REGENERADO — a correção    (57.080)
01/ago 14:27   covariáveis de ZUMBIDO reconstruídas da versão corrigida  ✓
03/ago 14:38   modelos ajustados, usando as covariáveis VELHAS de perda auditiva
03/ago 16:22   o teste roda
```

A correção veio **primeiro**. A análise rodou **três dias depois**, em cima de um arquivo construído
antes dela.

E dá para ver que alguém percebeu a regeneração: as covariáveis de zumbido foram reconstruídas no dia
seguinte, a partir da versão corrigida. O ramo de perda auditiva simplesmente nunca foi refeito.

Isso não é descuido. É o modo de falha mais comum que existe num pipeline com etapas manuais: dois
ramos paralelos, um atualizado e o outro não, sem nada declarando que um depende do outro.

### 4. Quem ficou de fora não é uma amostra aleatória

Era a pergunta que ninguém tinha testado, e a resposta importa:

```
taxa de exclusão na coorte inteira     0,73%

Leste Asiático   208 / 1.333  =  15,60%   ← 21× a taxa geral
Europeu          233 / 51.867 =   0,45%
Africano          49 / 14.927 =   0,33%
```

**Uma em cada seis pessoas de ancestralidade do Leste Asiático foi excluída.** Também há desvio em
sexo (mulheres 1,01% contra homens 0,43%) e forte em lote de sequenciamento (92% dos excluídos vêm do
lote 1).

O padrão de lote e data sugere causa técnica — sequenciados por exoma cedo, array nunca concluído
para uma parte. Não é nada sobre as pessoas.

Mas o efeito é de representatividade, não de poder estatístico: um grupo que já era 1,6% da coorte
perdeu um sexto dos seus. Isso não melhora aumentando o N.

---

## O que fechamos: dois cenários, em arquivo

A Fase 1 termina entregando **duas coortes** que as fases seguintes consomem sem re-derivar:

| | combined | EUR | AFR |
|---|---|---|---|
| **reprodução** (o que rodou) | 57.632 · 6.712 casos | 43.016 · 5.160 | 11.387 · 1.279 |
| **corrigido** (o que deveria) | 57.507 · 6.752 casos | 42.786 · 5.183 | 11.334 · 1.285 |

**Por que dois e não só o corrigido.** Se rodássemos só a versão corrigida, qualquer diferença nos
resultados seria ambígua: veio do defeito do pipeline, ou da coorte que *nós* mudamos? Perderíamos a
capacidade de dizer "reproduzimos os números dela" — que é o produto principal deste trabalho.

**A diferença entre os dois é o resultado.** Transforma "existe um defeito" em "o defeito move este
gene nesta magnitude" — ou em "o defeito é real e não muda nada", que é uma recomendação diferente
para o grupo.

Um detalhe que precisa ficar guardado: o braço de reprodução foi **copiado**, não reconstruído,
porque o arquivo de onde ele veio foi sobrescrito. O arquivo de covariáveis sobrevivente é o **único
registro do que foi analisado**. Se ele for apagado, esse cenário se perde.

O corrigido também tem um número bom de contar: ele **restaura 165 participantes do Leste Asiático**,
fazendo aquele grupo crescer 17,6% enquanto todos os outros mudam meio por cento.

---

## O que fica em aberto

| | |
|---|---|
| Por que o filtro usou o `.fam` do array | pergunta para Nikki e Elena |
| Por que o ramo de perda auditiva não foi reconstruído, se o de zumbido foi | idem |
| Os 4 controles legítimos que somem sem explicação | resíduo pequeno, não perseguido |
| Zumbido | declarado fora de escopo — verificamos só perda auditiva |
| Se o grupo restaura os 431 | decisão deles, não nossa |

---

## Em duas frases

O jeito de decidir quem tem perda auditiva está **certo** — reproduzido pessoa a pessoa. O jeito de
decidir quem entra no estudo está **errado** em três pontos independentes, todos pequenos, todos
empurrando o resultado para o lado conservador, e **nenhum detectado por nada** no pipeline.

O problema não é 0,59% dos casos ter sumido. É que sumiu e o pipeline reportou sucesso.
