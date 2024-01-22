# Mitochcheck-MAP-analysis

Utilizing the mean average precision (MAP) metric to assess  reproducibility and perturbation effect on single-cell profiles in the MitoCheck dataset.

## Installing Repo

To use this repo, clone the repository into your directory:
```bash
git clone
```

## Downloading the data

Within the `./data` folder, you can find a script named `download.py`.
Running this script facilitates the download of the MitoCheck data repository, accessible [here](https://zenodo.org/records/7967386).
The downloaded data, comprising both the training and control datasets, is stored in the `./data/raw` directory.
The training dataset contains labeled cells, each providing information about its phenotypic state (e.g., interphase, prophase, etc.).
Additionally, a log is generated within the `./data` directory to provide insights into the download process.

## Applying mAP to single-cell mitchech data

### Data Modifications

All the analysis is condcuted in the `./notebooks` directory.
The first notebooks that is used to analyize and generate the mAP scores is the `./notebooks/mitocheck-map-analysis.ipynb` notebook.

Firt we load in both the the negative controls and the training (phenotypically lalbed) dataset.
There were some minor modificiations to the dataset.
All phenotypes were stored in the `Mitocheck_Phenotypic_Class` in the trianing data.
Since there are not labeles in the negative control, we added the `Mitocheck_Phenotypic_Class` column that contains the phenotypic labeled and labeled the negative control cells as `neg_control`.

### Executing mAP

To perform the mAP analysis, the following parameters were used in the `copairs.map.run_pipeline()` function:

```python
pos_sameby = [
    "Mitocheck_Phenotypic_Class",
]
pos_diffby = ["Cell_UUID"]

neg_sameby = []
neg_diffby = ["Mitocheck_Phenotypic_Class"]

null_size = 1000
batch_size = 1000

# number of resampling
n_resamples = 10
```

These parameters allows `copairs` to distinguish cells belonging to different phenotypes.
The `Mitocheck_Phenotype_Class` is used to distinguish which phenotypes are presenting within the dataset.
The `Cell_UUID` enables the selection of cells and the formation of unique pairs both within and between the same phenotypic groups.

Further more, to mitigate label bias, a 1:1 ratio comparison was employed.
This ensures a balanced comparison between positive and negative instances.

However, we use the parameter `n_resamples`, which resamples the number of control cells that one wants to use when comparing a phenotype.

The processing pipeline begins by selecting cells and concatenating the labeled dataset with the control.
The dataset is then split into metadata and features.
Notably, the feature separation step is unique as it involves choosing from three different types of feature spaces: CellProfiler, DeepProfiler, and Both features, resulting in three distinct groups.

Each feature group undergoes two shuffling methods:

1. Phenotype shuffled: Phenotypic labels are shuffled within the concatenated feature space.
2. Feature space shuffled: Actual values within the feature spaces are shuffled around.

Subsequently, all of the data (shuffled and not shuffled) is fed into the `copairs.map.run_pipeline()` function, and the single-cell average precision scores (AP) are stored within the `./data/processed/sc_ap_scores` directory.

Once all the single-cell AP scores are computed, we aggregate the scores for cells using the `CellUUID` calculating the mean of the AP scores.
It's important to note that this mean does not represent the overall Average Precision (mAP) score.
Due to the application of subsampling techniques, each labeled cell appears multiple times in the training dataset, depending on the number of subsampling iterations.
In this analysis, we conducted 10 subsamples, resulting in labeled cells appearing 10 times.
The aggregated results are stored in the `./data/processed/sc_ap_scores/merged_sc_agg_ap_scores.csv` file.

Next, we utilize the `copairs.map.aggregate()` function to further aggregate the scores based on phenotypes.
The `sameby` parameter is set to `Mitocheck_Phenotype_Class`, and a threshold of `0.5` is applied.
This aggregation process is performed on the `merged_sc_agg_ap_scores.csv` file, grouping single-cell average precision scores based on phenotypes and assessing whether they surpass the p-value threshold.
The outcome of this step is the generation of mAP score files located in `./data/processed/aggregate_mAPs/sc_mAP_scores.csv`.

Finally, using these two files, we create visual representations in the `./figures` directory.
Bar plots showcase the mAP scores, while box plots illustrate the distribution of single-cell average precision scores.
These figures provide a comprehensive overview of the analysis results.
