# Fase 3 — o que é descontado · resumo

**Cópia de trabalho pessoal.** Detalhamento completo: [`phase_3/results/FINDINGS.md`](../phase_3/results/FINDINGS.md) ·
PCs explicados do zero: [`07_pcs_explicados.pt.md`](07_pcs_explicados.pt.md) ·
Notebook com os gráficos: [`phase_3/scripts/01_pc_selection.pt.ipynb`](../phase_3/scripts/01_pc_selection.pt.ipynb)

---

## Por que esta fase existe

Grupos populacionais diferem em milhões de variantes por razões históricas, nada a ver com doença. Se
um gene for mais comum num grupo que por acaso tem mais perda auditiva, aparece associação falsa.

Os **componentes principais** medem ancestralidade para descontá-la. A pergunta é **quantos usar** —
de menos deixa confusão residual, de mais gasta graus de liberdade e pode introduzir viés.

## O que fizemos

Auditamos os autovalores das três PCAs do estudo e refizemos duas delas.

## O que encontramos

| coorte | rodou com | scree sustenta | fonte |
|---|---:|---:|---|
| combined | 5 | **5** ✓ | PCA do release |
| EUR | 9 | **4** | PCA dentro da ancestralidade |
| AFR | 10 | **3** | idem |

**O combined estava certo.** EUR e AFR usaram mais que o dobro.

Joelho = primeiro PC cuja queda fica abaixo de 10%. Nenhum caso é de fronteira — mover o limiar para
8% ou 12% não muda resposta nenhuma.

**E ela planejou o método certo.** Os scree plots são dela, de 14/jul, e o plano dizia *"vou usar um
gráfico de variância explicada — provavelmente 4 ou 5"*. O método foi planejado, o dado foi gerado, e
os números usados não vêm dele.

**Um bloqueio que obrigou a refazer a PCA.** A dela rodou **depois** do corte do `.fam`, então 203
pessoas do EUR e 34 do AFR restauradas pela Fase 1 não tinham PC. Refizemos sobre todos com exoma,
com os parâmetros dela inalterados — o AFR reproduz o conjunto LD-pruned dela com **2 variantes de
diferença em 114.696**.

## Como seguimos

```
reprodução :  5 / 9 / 10   fixo, é o que rodou
corrigido  :  5 / 4 /  3   confirmado por dois caminhos independentes
REVEL      :  0,5 primário + 0,6 sensibilidade   (o que a reunião decidiu)
```

Covariáveis em `phase_3/results/covariates/`, colunas iguais às dela: `IID PHENO AGE AGE2 SEX Batch PC1..PCn`

| | reprodução | corrigido |
|---|---|---|
| combined | 57.632 · 5 PCs | 57.498 · 5 PCs |
| EUR | 43.016 · 9 PCs | **42.779 · 4 PCs** |
| AFR | 11.387 · 10 PCs | **11.334 · 3 PCs** |

O corrigido perde **9 pessoas sem idade registrada** — e esse é o contraste que vale guardar: idade é
covariável que o modelo **usa**, então excluir quem não tem é correto. Os 517 da Fase 1 caíram por
falta de PC imputado, que o modelo **nunca consome**. Mesmo mecanismo, justificativa oposta.

## Em aberto

- os scree justificam 9 e 10? — pergunta para a Nikki, embora nosso dado diga que não
- **anomalia do AFR:** nosso PC1/PC2 dá 12,8, o dela 1,8. Outliers e normalização descartados; o
  `eigenvec` bruto dela não sobreviveu, então não dá para saber sobre o que foi calculado. **Não
  aparece no EUR**, onde os PC1 batem (194,24 contra 197,37)
- a PCA dela foi bem construída? — esta fase checa se os PCs foram **usados** certo, não se a PCA é sólida
