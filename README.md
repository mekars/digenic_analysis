## Introduction
A case-control digenic analysis method including codes to:
- Convert input plink files into .raw format and prepare intermediate files to extract the binary combinations of input variants
- Generate binary combinations of input variants per sample
- Extract variant pairs observed in the proband (child) but not in their parents
- Extract gene pairs from variant pairs and select those to be used in burden testing
- Perform a covariate-adjusted burden test using Firth's logistic regression

## Dependencies
- `python 3.7.3`
- `R version 4.2.0`
- `plink 1.9`

## Input
- PLINK Input Files: **input.bed, input.bim, input.fam**
- Variant IDs and Gene Names File: **variantIDs_genes.txt**, sorted according to chromosome and position
```sh
1:970725:G:A    PLEKHN1
1:970735:G:A    PLEKHN1
1:1048963:G:A   AGRN
1:1087151:G:A   C1orf159
1:1087504:C:G   C1orf159
1:1087541:C:T   C1orf159
1:1087552:C:T   C1orf159
1:1087577:C:T   C1orf159
1:1181796:G:A   TTLL10
1:1184020:G:A   TTLL10
.
.
.
```
- Phenotype Data File: **phenotypes_with_pcs.txt** with columns containing sample ID, PC1, PC2, PC3 and case (1) - control (0) status, respectively

```sh
Sample1 -0.43889992     -1.362575357    0.536021629     1
Sample2 -0.420863071    -1.256101564    0.55213641      0
Sample3 -0.495990495    -1.377131679    0.626447258     1
.
.
.
```

## Output
- **final_burden_testing_result.txt**

```sh
GenePair   P  OR  L95  U95  CaseCount  ControlCount
DNAH1_DNAH10  0.418048892880962  1.42330896590859  0.609112484936727  3.51712165598446  6  5
DNAH1_FAT1  0.304256891852805  0.633636676912625  0.252970078481912  1.51006718098772  3  20
HSPG2_DNAH10  0.367147809116862  1.55106111980962  0.603235609632432  4.35220944565729  9  4
FAT2_DNHD1  0.305292398260725  1.70214961514179  0.621912692044574  5.20803707504133  12  10
.
.
.
```

