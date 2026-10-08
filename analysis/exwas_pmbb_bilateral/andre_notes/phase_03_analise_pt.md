# Fase 3 — o que é descontado · análise

**Autor:** Andre Rico · **Data:** 2026-10-08
**Cópia de trabalho pessoal, em português.** Não é página de Confluence.
Evidência herdada: [`phase_3/results/FINDINGS.md`](../phase_3/results/FINDINGS.md) (documento da replicação, copiado) ·
Premissas: [`PREMISES.md`](../PREMISES.md) · Fase anterior:
[`phase_02_analise_pt.md`](phase_02_analise_pt.md)

---

## O que esta fase faz

A Fase 1 separou **pessoas**, a Fase 2 separou **variantes**. Esta responde a terceira pergunta:
**o que mais, além do gene, explica quem tem perda auditiva?**

Porque se você não descontar isso, encontra o gene errado.

> **A PCA desta fase foi copiada**, não refeita — ela não depende do fenótipo, e explico por quê
> mais abaixo. As **covariáveis** foram construídas aqui, porque dependem.

---

## O problema que esta fase resolve

Imagine que você acha um gene associado à perda auditiva. Antes de comemorar:

**Idade.** Perda auditiva é fortemente ligada à idade. Se o seu grupo de casos for mais velho que o
de controles — e vai ser, porque gente mais velha tem mais diagnóstico de tudo — qualquer variante
mais comum em idosos vai parecer associada à surdez. Não é o gene; é a idade.

**Ancestralidade.** Esse é o clássico. Grupos populacionais diferem em **milhões** de variantes, por
razões históricas — migração, isolamento, deriva genética. Nada disso tem a ver com doença. Mas se
um grupo tiver, por qualquer motivo, mais diagnósticos de surdez, **toda variante mais comum nesse
grupo** vai aparecer associada.

Isso se chama **confusão por estrutura populacional**, e é o modo clássico de uma análise genética
produzir achado falso.

**Lote de sequenciamento.** O PMBB sequenciou as amostras em três lotes, em momentos diferentes.
Lotes diferentes podem ter viés técnico — um reagente, uma versão de software. Se casos e controles
não estiverem igualmente distribuídos entre os lotes, o viés técnico vira sinal falso.

**Parentesco.** Se dois irmãos estão na amostra, compartilham genes *e* tendem a compartilhar
diagnóstico, por motivos que não são o gene testado. Isso é tratado no passo seguinte, dentro do
SAIGE — mas é o mesmo problema.

A solução para os três primeiros é a mesma: **medir e pôr no modelo como covariável**, para o teste
perguntar "esse gene explica algo *além* disso?".

---

## Conceito — o que é um componente principal

Idade e sexo são fáceis de medir. Ancestralidade não — não existe uma coluna "ancestralidade" com um
número.

O que existe são as variantes genéticas da pessoa. Com 70.925 pessoas e meio milhão de variantes
comuns, é uma tabela de 70.925 × 500.000. Impossível de usar como covariável.

Mas a informação ali é **muito redundante**. Variantes próximas no genoma são herdadas juntas, e
populações inteiras carregam padrões comuns. A **análise de componentes principais (PCA)** comprime
essa redundância: procura as **direções de maior variação** na tabela e resume cada pessoa em poucos
números.

Uma analogia: para dizer onde alguém mora no Brasil você não precisa da lista de todas as ruas por
onde já passou. Latitude e longitude bastam. Os componentes principais são as coordenadas
geográficas do genoma.

**E ninguém diz à PCA o que é ancestralidade.** Ela só procura as maiores fontes de variação — e
acontece que, num conjunto de genomas humanos, a maior fonte de variação *é* ancestralidade. Ela
descobre sozinha.

---

## Quantos componentes usar, e como se decide

Esta é a única decisão real da fase.

**De menos** deixa confusão residual — sobra estrutura populacional não descontada, e ela vira
achado falso. **De mais** também custa: cada covariável gasta um grau de liberdade, então componentes
que não explicam nada reduzem o poder de detectar o que você procura. E pior — um componente do
"platô" pode correlacionar com o fenótipo por acaso, **introduzindo** viés em vez de remover.

O critério é o **autovalor**: quanto cada componente explica. Eles caem rápido no começo e depois
achatam. O ponto onde achata é o "joelho", e é onde se corta.

### Os números desta análise

```
EUR   autovalores: 194,2  42,4  29,1  13,4  12,4  11,9 ...
      quedas:       78,2%  31,2%  53,9%   8,0%   4,0%
                                          └── joelho: a partir daqui achatou

AFR   autovalores:  94,0   7,3   6,0   5,9   5,5   5,3 ...
      quedas:       92,2%  18,6%   1,4%   6,0%   3,7%
                                   └── joelho
```

**Adotado: 5 para o combined, 4 para o EUR, 3 para o AFR.**

O joelho é o primeiro PC cuja queda fica abaixo de 10%. Nenhum caso é de fronteira — mover o limiar
para 8% ou 12% não muda resposta nenhuma.

O **combined** usa 5 e vem de fonte diferente: os PCs da PCA global da release. Faz sentido, porque
ali a estrutura a descontar é **entre** ancestralidades, que é o que uma PCA global mede. Dentro do
EUR ou do AFR, um PC global está ocupado separando continentes e diz quase nada sobre estrutura
interna — por isso esses dois usam PCA feita dentro da própria ancestralidade.

> O histórico de como se chegou a 5/4/3, incluindo o que o estudo original usava, está em
> [`elena_replication/docs/07_pcs_explicados.pt.md`](../../elena_replication/docs/07_pcs_explicados.pt.md).

---

## Por que a PCA está copiada e não refeita

Porque ela **não depende de quem é caso**.

A PCA rodou sobre **todo mundo com exoma** em cada ancestralidade — 51.867 no EUR, 14.927 no AFR —
não sobre a coorte de análise. Ancestralidade é propriedade da pessoa, não do diagnóstico dela.
Mudar o fenótipo não muda a ancestralidade de ninguém.

Então refazer aqui produziria **exatamente os mesmos arquivos**. Foram copiados e conferidos byte a
byte ([`PROVENANCE.md`](../PROVENANCE.md)).

**Esses mesmos arquivos têm um segundo uso.** O conjunto de variantes podado por LD que a PCA usa
serve também como **GRM** — a matriz que o SAIGE usa para medir parentesco no próximo passo:

```
ALL    38.833 marcadores   70.925 pessoas    GRM da coorte combinada
EUR    53.920 marcadores   51.867 pessoas
AFR   114.694 marcadores   14.927 pessoas
```

---

## O que foi construído aqui: as covariáveis

Um arquivo por coorte, com as colunas que o SAIGE consome:

```
IID   PHENO   AGE   AGE2   SEX   Batch   PC1 ... PCn
```

**Por que `AGE2` (idade ao quadrado).** O efeito da idade sobre a audição não é uma reta. Entre 20 e
30 anos quase nada acontece; entre 70 e 80, muito. Só com `AGE` o modelo assume que cada ano pesa
igual. Acrescentar `AGE²` deixa a curva encurvar.

**`SEX` e `Batch` são declarados categóricos.** Isso importa mais do que parece para o `Batch`: os
valores são 1, 2 e 3 — rótulos de lote, não quantidades. Sem declarar, o SAIGE lê como número e
assume que o lote 2 fica exatamente no meio entre o 1 e o 3. No AFR isso seria interpolar entre
grupos de 7.853, 638 e 2.048 pessoas. Lote é rótulo.

### O resultado

| coorte | pedidos | escritos | casos | controles | PCs |
|---|---:|---:|---:|---:|---:|
| combined | 53.919 | **53.910** | 3.164 | 50.746 | 5 |
| EUR | 40.150 | **40.143** | 2.547 | 37.596 | 4 |
| AFR | 10.539 | **10.539** | 490 | 10.049 | 3 |

**Nove pessoas saem no combined, sete no EUR — todas por idade não registrada.** E essa exclusão é
**correta**, o que vale contrastar com a da Fase 1: idade é covariável que o modelo **consome**, e
uma linha sem ela não pode ser ajustada. Diferente de excluir alguém por faltar um dado que o modelo
nunca toca.

---

## Em aberto

- **O GRM da coorte combinada é magro.** 38.833 marcadores, contra 53.920 no EUR e 114.694 no AFR.
  A causa está no filtro: `--hwe 1e-6` aplicado **entre** ancestralidades remove 176 mil variantes,
  onde o mesmo filtro dentro do EUR remove 19 mil. Desvio de Hardy-Weinberg numa amostra estruturada
  é esperado e não indica erro de genotipagem — o filtro está fazendo o trabalho errado ali. Decidir
  como tratar HWE em amostra multi-ancestral é pergunta para o time.
- **Os 15 `UNKNOWN` e os outros grupos pequenos.** EAS, SAS, AMR e os não classificados existem só na
  coorte `combined`. Não há N para braço próprio, e os PCs do combined são o que desconta a mistura.
- **A PCA do AFR tem um perfil estranho**, com o PC1 dominando muito mais do que o esperado. Outliers
  e normalização foram descartados como causa; não aparece no EUR. Registrado em
  [`phase_3/results/PCA_NOTES.md`](../phase_3/results/PCA_NOTES.md), sem explicação até agora.
