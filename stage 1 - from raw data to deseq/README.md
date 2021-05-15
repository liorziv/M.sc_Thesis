# From Esat to Normalized Files

Esat res - contains east res from 4 experiments

Esat assemble number of reads per window (gene)

- Pic 10 + 14(1/2)
- Kcl 10

## Code:

1.main_normalization_step#1_Gauss.R - load and normalize the esat data

2.load_and_normalize_esat_Gauss.R - helper functions for 1

3.data_processing#2_Gauss.R - create some structs from the normalized data Mean values,fc(fold change) values etc..