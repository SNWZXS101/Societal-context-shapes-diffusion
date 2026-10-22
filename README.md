# Societal context shapes diffusion: code and figure inputs

This repository groups the manuscript code and required input files by analysis stage and figure. Figure folders such as `Fig2a`, `Fig3c`, and `FigS7` contain the source code and data selected for that figure. Rendered images and PDFs are excluded. The three PowerPoint source files for the conceptual `Fig1` are included.

| Directory | Contents |
|---|---|
| `PINN-SMFs-of-210waves/` | Existing repository content for the wave-level PINN analysis. |
| `210waves/` | Dynamic susceptible fraction, Extended Data Fig. 1, Fig. S1, and Table S1. |
| `feature extraction/` | Main feature extraction pipeline and the selected feature MAT file. |
| `Figure1/Fig1/` | Conceptual figure source files, including PowerPoint, MATLAB, and MAT files. |
| `FPCA and cluster/` | Final FPCA and clustering workflow, Fig. 2a–g, and supplementary diagnostics. |
| `Hierarchical Bayesian modelling/` | Covariate integration, main HBM, burden model, Fig. 3a–f, and Extended Data Figs. 2–8. |
| `Counterfactuals and target-guided control/` | Counterfactual workflow and Extended Data Fig. 9. |
| `Provinces Ranking/` | Province scoring workflow, Fig. 4a–e, and Extended Data Figs. 10–11. |
| `Uncertainty and robustness/` | PINN training sensitivity and Figs. S2–S14. |

The local `210waves/210波数据集` and `210waves/逐波PINN` folders were intentionally omitted from this upload. The existing `PINN-SMFs-of-210waves/` directory and outbreak-distribution workbook were retained. Ten Chinese directory names in the uploaded copy were translated to English; file contents were not edited. Some scripts still refer to their original local paths, so running them from this repository may require path configuration.