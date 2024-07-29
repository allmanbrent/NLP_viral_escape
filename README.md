# NLP_viral_escape
Code and instructions for the analysis in the manuscript "A systematic evaluation of the language-of-viral-escape model"

# Running the language models
All of the code that was used to generate grammaticality and semantic change values using the bespoke Hie BiLSTM model and [ESM2](https://github.com/facebookresearch/esm) can be found in the folder `language_models`. Most of the code was modified from [Hie et al (2021)](https://doi.org/10.1126/science.abd7331) and we will highlight scripts written and modified for this work. The results from this code can be found in `language_models/results/`

Install both the pip and conda environment by downloading the `environment.yml` and `requirements.txt` files.
```
# pip environment
pip install
# conda environment
```

## BiLSTM (Hie et al)

```
# run on the coronavirus spike single mutants
screen -d -m -L -Logfile cov_starr_screen.log python language_models/bin/cov_fasta_ESM2.py \
     data/cov/cov2_spike_wt.fasta \
     data/cov/starr_library.fasta \
     --model_name bilstm \
     --checkpoint models/cov.hdf5 \
     --dim 640 | tail -n+31 & \
     > cov_starr_fasta.log &

# run on the coronavirus Omicron-defining mutations
screen -d -m -L -Logfile cov_omicron_muts_screen.log python language_models/bin/cov_fasta_ESM2.py \
     data/cov/cov2_spike_wt.fasta \
     data/cov/wuhan_mutants.fasta \
     --model_name bilstm \
     --checkpoint models/cov.hdf5 \
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
     --checkpoint models/cov.hdf5 \
     --dim 640 --max_len 1274 | tail -n+31 & \
     > cov_Starr_ESM2_screen.log &

# run on the coronavirus Omicron-defining mutations
screen -d -m -L -Logfile cov_omicron_muts_ESM2.log python bin/cov_fasta_ESM2.py \
     data/cov/cov2_spike_wt.fasta \
     data/cov/wuhan_mutants.fasta \
     --model_name esm2_t30_150M_UR50D \
     --checkpoint models/cov.hdf5 \
     --dim 640 --max_len 1274 | tail -n+31 & \
     > cov_omicron_muts_ESM2_screen.log &

# run on the flu H1 HA mutations
screen -d -m -L -Logfile flu_ESM2.log python bin/cov_fasta_ESM2.py \
     data/flu/WSN1993_H1_HA.fa \
     data/flu/flu_H1_mutant_library.fasta \
     --model_name esm2_t30_150M_UR50D \
     --checkpoint models/flu.hdf5 \
     --dim 640 --max_len 1274 | tail -n+31 & \
     > flu_ESM2_screen.log &

# run on the hiv BG505 mutations
screen -d -m -L -Logfile hiv_ESM2.log python bin/cov_fasta_ESM2.py \
     data/hiv/hiv_env_translated.fasta \
     data/hiv/hiv_BG505_mutant_library.fasta \
     --model_name esm2_t30_150M_UR50D \
     --checkpoint models/hiv.hdf5 \
     --dim 640 --max_len 1274 | tail -n+31 & \
     > hiv_ESM2_screen.log &
```

### Masked grammaticality
We wrote a script to extract the grammaticality scores for each mutated sequence titled `masked_ESM2.py` and the script can be run in the environment using the command:
```

```
# Running the structural models

## MutComputeX

## MutRank

## Stability Oracle
