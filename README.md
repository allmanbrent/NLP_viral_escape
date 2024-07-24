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

## ESM2 (Lin et al)

### Masked grammaticality
We wrote a script to extract the grammaticality scores for each mutated sequence titled `masked_ESM2.py` and the script can be run in the environment using the command:
```

```
# Running the structural models

## MutComputeX

## MutRank

## Stability Oracle
