# Fase 1 — quem entra no estudo · análise

**Autor:** Andre Rico · **Data:** 2026-10-07
**Cópia de trabalho pessoal, em português.** Não é página de Confluence.
Versão canônica: [`phase_1/results/FINDINGS.md`](../phase_1/results/FINDINGS.md) ·
Premissas: [`PREMISES.md`](../PREMISES.md)

---

## O que esta fase faz

Separa os 70.925 participantes com exoma em três grupos: **casos**, **controles** e **excluídos**.

É a única fase desta análise que difere da replicação. Todo o resto — máscaras, PCA, SAIGE — é o
mesmo pipeline corrigido, sem alteração.

---

## De onde vêm os dados

Tudo da release institucional, nenhum arquivo intermediário de ninguém.

### Fonte 1 — diagnósticos por phecode

```
/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/
    PMBB-Release-2026-4.0_phenotype_conditions_phecode_x.txt
```

Puxamos toda linha cujo `condition_source_value` casa com `SO_39` seguido de um dígito — ou seja, a
**família auricular inteira**. São **655.946 linhas** (pessoa, data, código).

### Fonte 2 — a tabela `observation`

```
/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/
    PMBB-Release-2026-4.0_phenotype_observation.txt
```

Puxamos os códigos de tinnitus `388.3x` e `H93.1x`. São **25.094 linhas**.

**Por que essa segunda fonte existe:** a release v4 do PMBB *mudou de lugar* os códigos padrão de
tinnitus. Eles saíram dos arquivos de phecode e foram para a tabela `observation`. Quem monta a
família auricular olhando só os phecodes deixa **556 pessoas** com tinnitus invisíveis — e elas
acabam no grupo controle, que é justamente onde não podem estar.

### Fonte 3 — quem tem exoma, e de qual ancestralidade

```
/static/PMBB/PMBB-Release-2026-4.0/Exome/PCA/combined/
    PMBB-Release-2026-4.0_genetic_exome.commonsnps.samples_ancestries.tsv
```

As 70.925 pessoas com dado de exoma, e a coluna `Class` com a ancestralidade de cada uma.

### Repare que filtramos em duas camadas de código diferentes

Não é inconsistência — é imposição das tabelas:

| tabela | filtro | camada |
|---|---|---|
| `conditions_phecode_x` | `SO_39[0-9]` | **PhecodeX** |
| `observation` | `388.3x`, `H93.1x` | **ICD** |

A tabela `observation` **não tem coluna de phecode**. Ela guarda o código ICD cru. Então para
alcançar aqueles 25.094 eventos de tinnitus não há escolha: ou se filtra por ICD ali, ou se perdem
as 556 pessoas.

A regra geral vale guardar: **a camada de código é determinada pela tabela, não pela sua
preferência.** Detalhes em [`00_pmbb_codes.md`](00_pmbb_codes.md).

### Fonte auxiliar — a tabela mestra de códigos

```
analysis/elena/rarevariantExWAS/PMBBv4_phecodex/
    PMBBv4_hearing_tinnitus_case_code_master.csv
```

Montada pela Nikki. Não é usada no cálculo — serve para **traduzir** phecode em código ICD legível,
e é de onde saiu a tabela de códigos mais abaixo. É a mesma tabela que ela ficou de mandar ao Doug
na reunião de 02/10.

---

## Conceito 1 — o que é um phecode

Um diagnóstico no prontuário vem como **código ICD**: `H90.3`, `389.18`. São milhares de códigos,
muito específicos, e dois códigos diferentes podem significar praticamente a mesma coisa.

Um **phecode** agrupa códigos ICD que descrevem a mesma condição, para a análise não ficar
fragmentada. O `SO_396` é o phecode de *perda auditiva*.

E ele tem filhos, que é onde a coisa fica interessante.

## Conceito 2 — os filhos do `SO_396` são dois eixos, não subtipos

Essa foi a descoberta que destravou a fase. Os filhos **não** são categorias mutuamente exclusivas —
são **dois eixos independentes**, e um código ICD carrega um de cada:

| eixo | códigos |
|---|---|
| **tipo** | `.1` condutiva · `.2` neurossensorial · `.3` mista · `.5` súbita idiopática |
| **lateralidade** | `.8` bilateral · `.9` unilateral |

Então `H90.3 "Sensorineural hearing loss, bilateral"` carrega três coisas ao mesmo tempo:
`SO_396` (é perda auditiva), `SO_396.2` (é neurossensorial) e `SO_396.8` (é bilateral).

Isso é o que torna a exigência do Doug — "bilateral neurossensorial" — expressável diretamente:
basta pedir `.2` **e** `.8`.

## Conceito 3 — a regra de 2

Um diagnóstico isolado no prontuário pode ser erro de digitação, suspeita descartada, ou código posto
para justificar um exame. Dois diagnósticos **em datas diferentes** são muito mais confiáveis.

Então: caso precisa do diagnóstico em **duas datas distintas**. Quem tem só uma fica de fora — nem
caso, nem controle.

---

## Os códigos que usamos

### Os seis que tornam alguém CASO

São os códigos ICD que carregam `.2` e `.8` ao mesmo tempo:

| vocabulário | código | descrição | pessoas com ≥1 evento |
|---|---|---|---:|
| ICD10CM | `H90.3` | Sensorineural hearing loss, bilateral | 3.964 |
| ICD9CM | `389.18` | Sensorineural hearing loss, bilateral | 1.044 |
| ICD9CM | `389.16` | Sensorineural hearing loss, asymmetrical | 439 |
| ICD10CM | `H91.13` | **Presbycusis, bilateral** | 155 |
| ICD9CM | `389.12` | Neural hearing loss, bilateral | 85 |
| ICD9CM | `389.11` | Sensory hearing loss, bilateral | 26 |

Repare no `H91.13 Presbycusis` — **presbiacusia é a perda auditiva relacionada à idade**, que é
exatamente o fenótipo de interesse do estudo. Ele qualifica, e é bom que qualifique.

### Códigos que NÃO tornam caso, e que entravam antes

| código | descrição | o que falta |
|---|---|---|
| `H91.90` | Unspecified hearing loss, unspecified ear | tipo e lado |
| `H90.5` | Unspecified sensorineural hearing loss | o lado |
| `H91.93` | Unspecified hearing loss, bilateral | o tipo |
| `H90.0` | Conductive hearing loss, bilateral | é condutiva, não neurossensorial |
| `H91.92` | Unspecified hearing loss, left ear | é unilateral |

Cada um perde um dos dois eixos. Juntos, são a maior parte dos 3.588 casos que saíram.

**Vale separar os 2.695.** Essas pessoas **têm** bilateral neurossensorial registrado — falham só na
regra de 2, por terem o diagnóstico em uma data única. É motivo diferente de não ter o diagnóstico, e
é o maior grupo isolado que a restrição cria. Se algum dia a regra de 2 for afrouxada, é aí que se
olha primeiro. (A primeira versão deste documento somava os 2.695 com os 4.900 num número só, o que
escondia exatamente isso.)

### E o que desqualifica alguém como CONTROLE

**Qualquer** código da família auricular — `SO_390` (otite), `SO_391`, `SO_392`, `SO_394`,
`SO_397` (tinnitus), `SO_398` — mais o tinnitus da tabela `observation`.

Controle tem que ser alguém **sem nenhum sinal de problema de ouvido**.

---

## O resultado

```
70.925  participantes com exoma
 ───────
  3.164  CASOS       bilateral neurossensorial, em ≥2 datas
 50.755  CONTROLES   nenhuma evidência auricular
  2.695  excluídos   TÊM bilateral neurossensorial, mas em uma data só
  4.900  excluídos   têm perda auditiva, nunca bilateral neurossensorial
  9.411  excluídos   outra evidência auricular (tinnitus, otite, etc.)
```

Por ancestralidade:

| | casos | controles |
|---|---:|---:|
| EUR | 2.547 | 37.603 |
| **AFR** | **490** | 10.049 |
| EAS | 43 | 1.002 |
| AMR | 37 | 769 |
| SAS | 32 | 800 |
| `UNKNOWN0/1/2` | 15 | 532 |
| **total** | **3.164** | **50.755** |

As três linhas `UNKNOWN` são pessoas que a classificação de ancestralidade da release não conseguiu
atribuir a nenhum dos grupos de referência do 1000 Genomes. São 15 casos — pequeno, mas a tabela
precisa fechar com o total, senão quem somar as cinco primeiras acha 3.149 e fica procurando os 15
que faltam.

---

## O que isso custou, e por que vale

Os casos caíram de **6.752 para 3.164** — 53% a menos. É muita coisa, e perder metade dos casos
reduz o poder de detectar qualquer efeito.

A razão de fazer mesmo assim, nas palavras do Doug na reunião de 02/10:

> *"unilateral hearing loss we exclude, because it's less likely genetic, more likely environment
> related"*

Perda de um ouvido só costuma ter causa ambiental — trauma, infecção, exposição a ruído num lado.
Jogar isso no grupo de casos de um estudo **genético** adiciona ruído com cara de sinal. Menos casos,
mas casos que são de fato o que se está procurando.

E resolveu uma confusão da reunião: o Doug contava ~4.000 no browser do PMBB usando o código de
bilateral neurossensorial; a tabela mostra `H90.3` em **3.964 pessoas**. Ele e o pipeline contavam
coisas diferentes, e nenhum dos dois estava errado sobre a própria pergunta.

---

## Um erro meu no caminho, que vale registrar

"Bilateral neurossensorial em duas datas" tem **três leituras possíveis**, e elas não concordam:

| | leitura | casos |
|---|---|---:|
| A | ≥2 datas com `.2` **ou** `.8`, entre quem tem os dois | 4.029 |
| B | ≥2 datas de `.2` **e** ≥2 de `.8`, contadas em separado | 3.259 |
| **C** | **≥2 datas em que a MESMA data tem `.2` e `.8`** | **3.164** |

Minha primeira implementação usou a **A**.

### Uma pessoa real que mostra o problema

`PMBB2551063181373`, nove consultas em 16 anos:

| data | códigos | o que significa |
|---|---|---|
| 2008-07-07 | `.1` `.8` | condutiva, bilateral |
| 2010-03-13 | `.1` `.8` | condutiva, bilateral |
| 2010-12-17 | `.2` `.3` `.9` | neurossensorial + mista, **unilateral** |
| 2015-11-15 | `.1` `.8` | condutiva, bilateral |
| 2019-04-26 | `.1` `.8` | condutiva, bilateral |
| 2019-12-06 | `.1` `.9` | condutiva, unilateral |
| 2020-01-24 | `.3` `.9` | mista, unilateral |
| 2023-09-04 | `.8` | bilateral, tipo não registrado |
| 2024-11-02 | `.1` `.8` | condutiva, bilateral |

**Essa pessoa nunca foi diagnosticada com perda bilateral neurossensorial.** Nenhuma consulta tem
`.2` e `.8` juntos. O que ela tem é condutiva bilateral seis vezes, e uma única consulta com
neurossensorial *unilateral*.

```
leitura A — datas que carregam .2 OU .8    →  8 datas  →  VIRA CASO
leitura C — datas em que a MESMA data tem os dois  →  0  →  não é caso
```

O erro é tratar `.2` e `.8` como evidências do mesmo diagnóstico, quando são eixos independentes. O
`.8` significa apenas *"bilateral"* — e perda **condutiva** bilateral também o carrega. Então cada
consulta de otite bilateral entra na contagem de uma definição escrita para excluir perda condutiva.
A pessoa junta `.8` por seis consultas de condutiva, pega um `.2` solto de outra doença, e a soma de
dois eixos colhidos em **momentos diferentes** monta um diagnóstico que nunca existiu.

A **B** tem o mesmo defeito de forma mais branda: ainda deixa uma consulta neurossensorial-unilateral
somar com uma condutiva-bilateral.

A **C** é o que a frase significa: em pelo menos duas datas distintas, **o diagnóstico registrado
naquela data** foi bilateral neurossensorial.

**A regra:** dois eixos só significam um diagnóstico quando aparecem juntos **na mesma linha do
prontuário**. Separados no tempo, são duas doenças diferentes somadas por engano. São **865 pessoas**
que a leitura frouxa incluiria.

Peguei isso comparando as três entre si. **Nada mais adiante no pipeline teria notado** — as três
produzem uma coorte plausível, rodam até o fim e dão p-valores.

---

## Em aberto

- **AFR com 490 casos.** A replicação, com 1.285, já marcava 26 dos 30 primeiros do AFR como
  frágeis. Com 490 piora. Vai rodar — não rodar é um viés de reporte próprio —, mas tem que ser lido
  como nulo sem poder, não como busca.
- **Audiogramas substituiriam esta fase inteira.** Código de diagnóstico é um substituto pobre para
  um audiograma, que é quantitativo e resolve tipo e lateralidade sem inferência nenhuma. Falta a
  ponte REDCap ↔ PMBB. Quando existir, **só esta fase muda.**
- **A súbita idiopática (`SO_396.5`) ficou de fora.** 203 pessoas da definição ampla carregam esse
  código. Se perda súbita idiopática pertence a um estudo genético de perda relacionada à idade é
  pergunta clínica, não minha.
