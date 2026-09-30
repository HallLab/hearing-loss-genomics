# Fase 2 — quais variantes contam

**Autor:** Andre Rico · **Data:** 2026-09-30 · **Status:** completa
**Cópia de trabalho pessoal, em português.** Não é página de Confluence. O registro formal está em
[`phase_2/results/FINDINGS.md`](../phase_2/results/FINDINGS.md).

---

## A pergunta

A Fase 1 decidiu **quem** entra no estudo. A Fase 2 decide **quais variantes genéticas contam**.

Parece detalhe técnico. É a decisão mais consequente do pipeline inteiro, e vale entender por quê.

---

## Por que existem "máscaras"

Cada variante rara aparece em pouquíssima gente — às vezes uma pessoa só. Não há poder estatístico
nenhum para testar uma variante dessas isoladamente.

A saída é **agrupar por gene**: em vez de perguntar *"quem tem esta variante?"*, pergunta-se *"quem
tem **qualquer** variante danosa neste gene?"*. Aí os números somam e o teste fica possível.

Mas isso obriga a decidir o que é "danosa". Essa decisão é a **máscara** — uma lista, por gene, de
quais variantes entram no teste.

O pipeline monta quatro:

| máscara | o que agrupa |
|---|---|
| `pLOF` | variantes que **quebram** a proteína |
| `pDM` | trocas de aminoácido previstas como danosas |
| `pLOF_pDM` | as duas juntas |
| `ALL` | todas as variantes raras, sem filtro |

**Se a máscara estiver errada, todo resultado do gene está errado** — e não há como perceber olhando
o resultado, porque ele sai com aparência normal.

---

## O que significa "quebrar a proteína"

Um gene é lido em palavras de três letras, que viram a sequência de aminoácidos da proteína. Há três
jeitos clássicos de arruinar essa leitura:

- **stop prematuro** — a leitura para no meio, e sai meia proteína
- **frameshift** — uma letra a mais ou a menos desloca todas as palavras seguintes, e o resto sai como
  lixo
- **erro de splicing** — o gene vem com pedaços (íntrons) que precisam ser recortados antes da
  leitura. Se o ponto de corte quebra, o íntron fica e a proteína sai errada

Esses três são **perda de função** — pLOF, de *predicted loss of function*.

### E aqui mora o problema

Os pontos exatos de corte se chamam **doador** e **aceitador**. Variante neles = corte quebrado, com
alta probabilidade.

Mas existem variantes **perto** do ponto de corte, sem estar nele. O VEP — o programa que anota o que
cada variante faz — tem nomes específicos para elas:

| anotação | o que é | impacto VEP |
|---|---|---|
| `splice_donor_variant` | **no** ponto de corte | **HIGH** |
| `splice_acceptor_variant` | **no** ponto de corte | **HIGH** |
| `splice_polypyrimidine_tract_variant` | numa região que ajuda a marcar o corte | **LOW** |
| `splice_donor_region_variant` | perto do doador, fora dele | **LOW** |
| `splice_donor_5th_base_variant` | quinta base depois do doador | **LOW** |

As três de baixo **podem ou não** atrapalhar o splicing. Na maioria dos casos não atrapalham.

É exatamente para decidir esses casos que existe o **SpliceAI** — um preditor que dá uma nota de 0 a 1
para "esta variante atrapalha o corte?". É por isso que a definição do próprio plano de análise diz:

> pLOF = frameshift, stop-gained, start-lost, stop-lost, **ou variantes de splice (SpliceAI ≥ 0,2)**

Ou seja: as de impacto baixo entram **só se** o SpliceAI concordar.

---

## O que encontramos

### O número

No cromossomo 8, de 197.002 registros marcados como pLOF pelo pipeline:

```
 50.252  têm perda de função de verdade         25,5%
146.750  não têm — só uma anotação de splice
             de impacto BAIXO                   74,5%
      0  não têm nenhuma das duas
```

Zero exceções — nada entra na máscara por outro caminho.

**Atenção a um detalhe que muda o número.** Os 197.002 acima são *registros de anotação*, e cada
variante tem um registro por transcrito. A máscara contém **variantes**, não registros. Contando
variantes, o chr8 dá **64,5%**.

Os dois números estão certos para o que medem: o de registro serve para diagnosticar a regra, o de
variante descreve a máscara. **É o de variante que se cita.**

### E vale para o genoma inteiro

```
1.089.876   variantes na máscara pLOF
  395.074   com perda de função real      36,2%
  694.802   sem                           63,8%
```

Por cromossomo a faixa é **60,9% a 65,7%** — todos os 22 dentro de cinco pontos. O chr8 está em
64,5%, quase exatamente a média. O defeito é **uniforme**, não concentrado em lugar nenhum.

**Quase dois terços da máscara pLOF do estudo — 694.802 variantes — não é perda de função.**

### O portão do SpliceAI nunca foi usado

A definição exige SpliceAI ≥ 0,2 para as de splice. Na prática:

```
129.118  ENTRARAM com SpliceAI abaixo de 0,2
  9.094  FICARAM DE FORA com SpliceAI acima de 0,2
```

O SpliceAI está **descorrelacionado** da decisão. Ele foi calculado, foi escrito no arquivo, e nunca
foi consultado.

### A causa, no código

```python
lof_terms = [
    "frameshift_variant", "stop_gained", "start_lost", "stop_lost",
    "splice_acceptor_variant", "splice_donor_variant",
    "splice_donor_5th_base_variant",         # ← impacto BAIXO
    "splice_donor_region_variant",           # ← impacto BAIXO
    "splice_polypyrimidine_tract_variant"    # ← impacto BAIXO
]
is_pLOF = any(term in str(Consequence) for term in lof_terms)
```

Os três termos de impacto baixo estão **explicitamente na lista**, um por linha. Não foi acidente de
programação — foi a lista que alguém escreveu.

E `is_pLOF` sai **só** dessa lista. Não há nenhuma menção a SpliceAI nessa linha.

Então são dois erros independentes:

1. três anotações que não são perda de função entraram na definição
2. o portão que deveria filtrá-las nunca foi ligado

Só o `splice_polypyrimidine_tract_variant` responde por 122.036 registros — mais que o **dobro** de
todo o conjunto de perda de função legítima.

---

## Como sabemos que isso é defeito e não discordância

Essa distinção importa, porque duas pessoas podem anotar as mesmas variantes e discordar
legitimamente — thresholds diferentes, versões diferentes de programa.

Fechei essa porta de três formas:

**1. A comparação começou contra o release.** O PMBB publica os próprios *group files* com as mesmas
categorias. Foi assim que o problema apareceu.

**2. Mas a confirmação é interna.** O arquivo de classificação do próprio pipeline tem a consequência
que **ele** atribuiu, o SpliceAI que **ele** calculou, e a decisão `is_pLOF` — tudo na mesma linha.
Não precisei de referência externa: o pipeline discorda de si mesmo.

**3. E a regra do código reproduz o arquivo em 100%.** Reimplementei a expressão e comparei:
**6.093.288 de 6.093.288 linhas**, zero divergências. A lógica que li é a que rodou.

### A comparação com o `pDM` reforça

A máscara `pDM` também diverge do release — 16.244 variantes que o release não considera danosas.

Mas ali é **desacordo de limiar**: quão danosa é uma troca de aminoácido depende de qual preditor e
qual corte se usa, e a própria reunião de 2026-07-01 deixou isso aberto (REVEL 0,5 ou 0,6).

**Discordar de quão danosa é uma missense é juízo. Chamar variante intrônica de perda de função não
é.** É esse contraste que torna o achado do pLOF sustentável.

---

## O que isso significa

**A máscara `pLOF` do estudo é majoritariamente variante intrônica perto de sítio de splice.** Quem
leu um resultado "pLOF" desse pipeline leu, em três quartos, outra coisa.

Direção do efeito: encher a máscara de variantes inofensivas **dilui** o sinal. Se um gene tem 5
variantes que realmente quebram a proteína e 15 que não fazem nada, o teste mistura as 20 e o efeito
some no ruído.

Ou seja — como na Fase 1, o erro empurra para o **nulo**. Não fabrica achado; esconde.

Mas ao contrário da Fase 1, **este não é pequeno.** Lá eram 0,59% dos casos. Aqui são quase dois
terços de uma máscara inteira, em todos os 22 cromossomos.

---

## O que está pendente

| | |
|---|---|
| ~~Se o chr8 é representativo~~ | **resolvido — é.** 22 cromossomos entre 60,9% e 65,7% |
| Quanto isso muda os resultados | é pergunta da Fase 4, não desta |
| As máscaras `pLOF_pDM` e `ALL` | herdam o `pLOF`, então herdam o defeito — não quantificado |
| Refazer a máscara corrigida | decisão de escopo, ainda não tomada |

---

## Em duas frases

A definição de pLOF no plano de análise está certa. A implementação incluiu três anotações de impacto
baixo que não são perda de função, e nunca ligou o filtro de SpliceAI que existia justamente para
decidir esses casos.

Isso faz quase dois terços da máscara pLOF — 694.802 variantes, em todos os 22 cromossomos — ser
algo que a própria definição do estudo não admitiria.
