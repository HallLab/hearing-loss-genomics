# Fase 2 — quais variantes contam · resumo

**Cópia de trabalho pessoal.** Detalhamento completo: [`phase_2/results/FINDINGS.md`](../phase_2/results/FINDINGS.md) ·
Página de revisão com referências de linha: [`06_plof_mask_defect.md`](06_plof_mask_defect.md) ·
Explicação em linguagem leiga: [`05_phase_2_fechamento.pt.md`](05_phase_2_fechamento.pt.md)

---

## Por que esta fase existe

Variante rara sozinha não tem poder estatístico, então o teste agrupa por gene. Isso obriga a decidir
o que é "danosa" — e essa decisão é a **máscara**. Variante inofensiva dentro do grupo **dilui o
sinal**: um gene com 5 variantes que quebram a proteína e 15 inertes é testado sobre as 20.

A máscara não é economia de cálculo. É o que concentra o sinal.

## O que fizemos

Comparamos as máscaras contra os *group files* que o próprio PMBB publica, depois confirmamos tudo
pelo arquivo de classificação interno do pipeline.

## O que encontramos — quatro defeitos, quatro causas

| defeito | afeta | tamanho | quem achou |
|---|---|---|---|
| três termos de splice de **impacto BAIXO** na lista de pLOF | `pLOF`, `pLOF_pDM` | — | nós |
| **portão do SpliceAI ≥ 0,2 calculado e nunca aplicado** | `pLOF`, `pLOF_pDM` | 63,8% das variantes sem perda de função | nós |
| **parsing do REVEL** — VEP escreve lista por transcrito, `to_numeric` virou NaN | **`pDM`**, `pLOF_pDM` | 5.391 onde deviam ser 20.618 (**73,9% perdido**) | **Nikki** |
| genes **não-codificantes** nas máscaras | `pLOF`, `pLOF_pDM` | 1.101 de 19.038 genes | **Nikki** |

O terceiro derruba o que tínhamos escrito. A Fase 2 afirmava que o `pDM` estava íntegro e que sua
divergência era desacordo de limiar. **Estava errado** — e não vimos porque comparamos de fora, onde
limiar diferente e parsing quebrado produzem o mesmo sintoma.

## Como seguimos

As quatro máscaras **reconstruídas da fonte** (`phase_2/results/masks_v2/`), não filtradas das dela —
o conserto do REVEL *acrescenta* variantes, e filtrar não acrescenta:

```
pLOF       1.002.120 →   462.144      pDM crescer é o REVEL recuperado
pDM          720.983 →   890.332  ↑
pLOF_pDM   1.717.883 → 1.340.936
ALL       21.415.507 → 15.951.289
```

Regra corrigida:

```
pLOF =   lista explícita: frameshift / stop_gained / start_lost / stop_lost
         / splice_acceptor / splice_donor / transcript_ablation
  OU     qualquer outra anotação de splice, mas só com SpliceAI ≥ 0,2
  E      BIOTYPE == protein_coding
```

Casamento de **termo exato**, não pelo campo `IMPACT` do VEP — que é por linha, e traria os termos de
impacto baixo de volta.

⚠️ **Uma decisão além dos consertos:** o filtro de biotype foi aplicado ao `ALL` também, que não tinha
defeito. Remove 6.532 genes. Então o `ALL` carrega **duas** diferenças entre braços, não uma.

## Em aberto

- de onde veio a lista `lof_terms` — três termos nomeados um a um não parecem esquecimento
- REVEL 0,5 vs 0,6 — pendente desde 2026-07-01, e agora importa mais
- o que isso faz com os resultados — **Fase 4**
