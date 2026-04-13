# CMAscan

CMAscan is a repository for curated chaperone-mediated autophagy (CMA) motif data and the notebooks used to analyse the dataset and derive scoring resources.

## Repository contents

- `dataset/cma_motif_dataset.csv` and `dataset/cma_motif_dataset.xlsx`
  - curated CMA motif dataset with literature annotations
- `dataset/dataset_seqs.fasta` and `dataset/dataset_hq_seqs.fasta`
  - FASTA files used for motif-to-sequence consistency checks
- `PSSM/cPSSM.pkl`
  - canonical PSSM built from validated canonical CMA motifs
- `PSSM/ePSSM.pkl`
  - expanded PSSM built from validated canonical and atypical CMA motifs
- `ExternalSoftware/iupred3.tar.gz`
  - bundled IUPred3 package used by the workflow
- `Notebooks/CMAscan_DB_analysis.ipynb`
  - main analysis notebook prepared for GitHub and Google Colab

## Dataset terminology

### Type
- `C`: canonical motif
- `A`: atypical motif

### Dataset class
- `P`: positive canonical
- `P*`: putative canonical
- `PA`: positive atypical
- `PA*`: putative atypical
- `E`: excluded

### Experimental support
- `MDM`: mutagenesis support flag (`Yes`, `No`, `In silico`)
- `NI`: not investigated
- `CoLoc`: co-localization
- `IsoLyso`: isolated lysosome assay
- `Lyso U`: lysosomal uptake
- `L2A Puncta`: LAMP2A-dependent puncta formation

## Notebook scope

`Notebooks/CMAscan_DB_analysis.ipynb` performs:

- dataset loading and reporting
- FASTA consistency checks against the curated motif table
- selection of high-quality motif subsets
- positional amino acid analysis
- background amino acid composition analysis
- cPSSM, ePSSM, and charge-based PSSM construction
- motif permutation scoring
- LOOCV benchmarking and TPR plotting

The notebook is written so it can run locally or in Google Colab. When run in Colab from GitHub, it fetches required files directly from this repository and downloads the reviewed human reference proteome from UniProt when needed.

## Google Colab

After the repository is pushed, the notebook can be opened in Colab with:

`https://colab.research.google.com/github/Marcin-Lubocki/CMAscan/blob/main/Notebooks/CMAscan_DB_analysis.ipynb`

Recommended usage in Colab:

1. Open the notebook from GitHub.
2. Choose `Runtime` -> `Run all`.
3. Wait for package installation and notebook execution to finish.

## External tools

CMAscan uses external tools that should be cited separately if relevant to a publication.

### MusiteDeep
- [https://www.musite.net/](https://www.musite.net/)
- Wang D. et al. 2020. MusiteDeep: a deep-learning based webserver for protein post-translational modification site prediction and visualization. *Nucleic Acids Research* 48(W1):W140-W146.
- Wang D. et al. 2019. Capsule network for protein post-translational modification site prediction. *Bioinformatics* 35(14):2386-2394.
- Wang D. et al. 2017. MusiteDeep: a deep-learning framework for general and kinase-specific phosphorylation site prediction. *Bioinformatics* 33(24):3909-3916.

### IUPred3
- [https://iupred3.elte.hu/](https://iupred3.elte.hu/)
- Erdős G., Pajkos M., Dosztányi Z. 2021. IUPred3: prediction of protein disorder enhanced with unambiguous experimental annotation and visualization of evolutionary conservation. *Nucleic Acids Research* 49(W1):W297-W303.
- Mészáros B., Erdős G., Dosztányi Z. 2018. IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. *Nucleic Acids Research* 46(W1):W329-W337.
- Erdős G., Dosztányi Z. 2020. Analyzing Protein Disorder with IUPred2A. *Current Protocols in Bioinformatics* 70(1):e99.

## Contact

If you find an issue in the dataset or notebook, contact: `marcin.lubocki@phdstud.ug.edu.pl`
