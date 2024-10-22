# NLP_viral_escape
Code and instructions for the analysis in the manuscript "A systematic evaluation of the language-of-viral-escape model"

# Running the language models
All of the code that was used to generate grammaticality and semantic change values using the bespoke Hie BiLSTM model and [ESM2](https://github.com/facebookresearch/esm) can be found in the folder `language_models`. Most of the code was modified from [Hie et al (2021)](https://doi.org/10.1126/science.abd7331) and we will highlight scripts written and modified for this work. The results from this code can be found in `language_models/results/`

Install conda environment by downloading the `environment.yml` file.
```
conda env create -f environment.yml
```

## BiLSTM (Hie et al)

Note that while the script has "ESM2" in the title, the `model_name` parameter defines which model will be used. 
```
# run on the coronavirus spike single mutants
screen -d -m -L -Logfile cov_starr_screen.log python language_models/bin/cov_fasta_ESM2.py \
     data/cov/cov2_spike_wt.fasta \
     data/cov/starr_library.fasta \
     --model_name bilstm \
     --checkpoint language_models/models/cov.hdf5 \
     --dim 640 | tail -n+31 & \
     > cov_starr_fasta.log &

# run on the coronavirus Omicron-defining mutations
screen -d -m -L -Logfile cov_omicron_muts_screen.log python language_models/bin/cov_fasta_ESM2.py \
     data/cov/cov2_spike_wt.fasta \
     data/cov/wuhan_mutants.fasta \
     --model_name bilstm \
     --checkpoint language_models/models/cov.hdf5 \
     --dim 640 | tail -n+31 & \
     > cov_omicron_muts.log &
```



## ESM2 (Lin et al)

```
# run on the coronavirus Starr mutations
screen -d -m -L -Logfile cov_Starr_ESM2.log python bin/cov_fasta_ESM2.py \
     data/cov/cov2_spike_wt.fasta \
     data/cov/starr_library.fasta \
     --model_name esm2_t30_150M_UR50D \
     --checkpoint language_models/models/cov.hdf5 \
     --dim 640 --max_len 1274 | tail -n+31 & \
     > cov_Starr_ESM2_screen.log &

# run on the coronavirus Omicron-defining mutations
screen -d -m -L -Logfile cov_omicron_muts_ESM2.log python bin/cov_fasta_ESM2.py \
     data/cov/cov2_spike_wt.fasta \
     data/cov/wuhan_mutants.fasta \
     --model_name esm2_t30_150M_UR50D \
     --checkpoint language_models/models/cov.hdf5 \
     --dim 640 --max_len 1274 | tail -n+31 & \
     > cov_omicron_muts_ESM2_screen.log &

# run on the flu H1 HA mutations
screen -d -m -L -Logfile flu_ESM2.log python bin/cov_fasta_ESM2.py \
     data/flu/WSN1993_H1_HA.fa \
     data/flu/flu_H1_mutant_library.fasta \
     --model_name esm2_t30_150M_UR50D \
     --checkpoint language_models/models/flu.hdf5 \
     --dim 640 --max_len 1274 | tail -n+31 & \
     > flu_ESM2_screen.log &

# run on the hiv BG505 mutations
screen -d -m -L -Logfile hiv_ESM2.log python bin/cov_fasta_ESM2.py \
     data/hiv/hiv_env_translated.fasta \
     data/hiv/hiv_BG505_mutant_library.fasta \
     --model_name esm2_t30_150M_UR50D \
     --checkpoint language_models/models/hiv.hdf5 \
     --dim 640 --max_len 1274 | tail -n+31 & \
     > hiv_ESM2_screen.log &
```

### Masked grammaticality
We wrote a script to extract the grammaticality scores for each mutated sequence titled `masked_ESM2.py` and the script can be run in the environment using the command:
```
python language_models/ESM_probability_calcs.py <reference.fasta>
```
This script can be used on any protein reference sequence to compute the masked probabilities at every site along the sequence. The reference sequences used in this study are available at [data/cov/cov2_spike_wt.fasta](https://github.com/allmanbrent/NLP_viral_escape/blob/main/data/cov/cov2_spike_wt.fasta), [data/flu/WSN1933_H1_HA.fa](https://github.com/allmanbrent/NLP_viral_escape/blob/main/data/flu/WSN1933_H1_HA.fa), [data/hiv/BG505_env_translated.fasta](https://github.com/allmanbrent/NLP_viral_escape/blob/main/data/hiv/BG505_env_translated.fasta).

# Running the structural models

## MutComputeX
The original Github repository for MutComputeX ([d'Oelsnitz et al. 2024](https://www.nature.com/articles/s41467-024-46356-y)) can be found at [https://github.com/danny305/MutComputeX](https://github.com/danny305/MutComputeX). The inference pipeline for MutComputeX can be found here: [https://github.com/danny305/MutComputeX/blob/master/scripts/generate_predictions.py](https://github.com/danny305/MutComputeX/blob/master/scripts/generate_predictions.py).

## MutRank
MutRank ([Gong et al. 2024](https://openreview.net/forum?id=XblaAN1jq6&referrer=%5Bthe%20profile%20of%20Tianlong%20Chen%5D(%2Fprofile%3Fid%3D~Tianlong_Chen1))) is currently under peer review. The inference pipeline can be found here: [https://github.com/danny305/StabilityOracle/blob/master/scripts/run_stability_oracle.py](https://github.com/danny305/StabilityOracle/blob/master/scripts/run_stability_oracle.py)

## Stability Oracle
The original Github repository for Stability Oracle ([Diaz et al. 2024](https://www.nature.com/articles/s41467-024-49780-2)) can be found at [https://github.com/danny305/StabilityOracle](https://github.com/danny305/StabilityOracle).
