# CMAscan

**Detection and scoring of CMA-targeting motifs in protein sequences.**

CMAscan is a browser-based pipeline for screening protein sequences for chaperone-mediated autophagy (CMA) targeting motifs. The tool runs entirely in Google Colab — no installation required.

## Run CMAscan

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Marcin-Lubocki/CMAscan/blob/main/Notebooks/CMAscan_colab.ipynb)

> *Click the badge above to open the notebook in Google Colab. Runs in your browser; no installation needed.*

## What CMAscan does

For each protein sequence in your input, CMAscan reports:

- All candidate CMA-targeting pentamer motifs (canonical and PTM-dependent variants)
- Two PSSM scores per motif:
  - **cPSSM** — score against canonical CMA motifs
  - **ePSSM** — score against canonical + atypical CMA motifs (extended PSSM)
- **IUPred3** disorder score for each motif's local sequence context
- For PTM-dependent motifs: simulated post-translational modification (acetyl K→Q; phospho S/T/Y→E or D) and re-scoring of the modified motif

The output is a CSV file plus a set of summary visualizations (motif type frequencies, score distributions, disorder–score relationships, an interactive 3D map of motif hits, and per-protein density heatmap).

## Quick start

1. Click the **Open in Colab** badge above.
2. Click **Runtime → Run all**.
3. Check **"Use demo input"** in Step 1 to try CMAscan with a built-in
   example (no upload needed) — recommended for first-time users.
4. To analyze your own sequences, leave "Use demo input" unchecked and
   upload your FASTA file when prompted.
5. To enable disorder scoring, check **"Use IUPred3"** and upload
   `iupred3.tar.gz` when prompted (download from
   https://iupred3.elte.hu/ — academic license, free).
6. Results appear at the end of the notebook and download automatically.

For details, troubleshooting, and parameter descriptions, see the notebook itself — every step is documented inline.

## Input format

A standard FASTA file with one or more protein sequences. Sequences must contain only the 20 standard amino acids (plus optional `X` for unknown residues). CMAscan strips whitespace and normalizes case automatically; sequences with non-standard residues are reported with a clear error message.

## Output

A single CSV file (`motif_hits_scored.csv`) with one row per motif candidate. Columns:

| Column | Description |
|---|---|
| `protein_name` | Identifier from FASTA header |
| `mer` | 5-residue motif candidate |
| `type` | `canonical`, `phospho`, `acetyl`, or `No motifs found` |
| `localization` | 1-based start position in the protein sequence |
| `cPSSM` | Score from canonical PSSM |
| `ePSSM` | Score from extended PSSM (canonical + atypical motifs) |
| `IUPred3` | Mean disorder score over the 5 motif residues (0 = ordered, 1 = disordered) |
| `mer_PTM_mimetic` | Motif sequence after simulated PTM |
| `cPSSM_PTM` | cPSSM score for the PTM mimetic |
| `ePSSM_PTM` | ePSSM score for the PTM mimetic |
| `IUPred3_PTM` | Disorder score recomputed on the PTM mimetic sequence |

## Repository contents

```
CMAscan/
├── Notebooks/
│   ├── CMAscan_colab.ipynb                              Main pipeline notebook
│   ├── CMAscan_DB_analysis.ipynb                        Motif database analysis
│   └── CMAscan_DB_analysis_sequence_and_structural_     Validation analysis
│       features_validation.ipynb
├── PSSM/
│   ├── cPSSM.pkl                                        Canonical PSSM matrix
│   └── ePSSM.pkl                                        Extended PSSM matrix
├── Results/                                             Pipeline outputs (gitignored)
├── dataset/
│   ├── cma_motif_dataset.csv                            Curated CMA motif dataset
│   ├── cma_motif_dataset.xlsx
│   ├── demo_input.fasta                                 Example input (25 proteins)
│   ├── motif_structures_sasa_v0.6.xlsx
│   └── motif_structures_v0.6.csv
├── LICENSE                                              MIT License
├── README.md
└── requirements.txt
```

## Dependencies

CMAscan installs its dependencies automatically when you run the notebook. For reference, these are:

- BioPython
- pandas, NumPy, SciPy
- matplotlib, seaborn, plotly
- tqdm
- IUPred3 (optional external tool)

IUPred3 is NOT bundled in this repository (license restriction). See Quick Start for download instructions.

## External tools and citations

CMAscan uses [IUPred3](https://iupred3.elte.hu/) as an optional external tool, downloaded separately by the user under their own academic license from https://iupred3.elte.hu/. **If you publish results obtained with CMAscan, please cite the IUPred3 authors:**

- Erdős G, Pajkos M, Dosztányi Z. IUPred3: prediction of protein disorder enhanced with unambiguous experimental annotation and visualization of evolutionary conservation. *Nucleic Acids Research* 2021;49(W1):W297–W303.
- Mészáros B, Erdős G, Dosztányi Z. IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. *Nucleic Acids Research* 2018;46(W1):W329–W337.
- Erdős G, Dosztányi Z. Analyzing Protein Disorder with IUPred2A. *Current Protocols in Bioinformatics* 2020;70(1):e99.

## Data and software availability

The CMA motif dataset and all analysis scripts used in the accompanying publication are openly available in this repository. CMAscan is implemented as a Google Colab notebook and runs directly in the browser; no installation is required.

## License

MIT License. See [LICENSE](LICENSE) for the full text.

## Contact

For questions, bug reports, or feature requests, please open an issue on this repository.
