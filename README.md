**System Requirements**<br>
All scripts were executed using R version 4.4.1 (2024-06-14) and Python 3.12.4. Details about package versions are provided in the paper.

**Installation Guide**<br>
There is no specific installation guide other than installing common R or Python packages. Each package should be installed according to the instructions provided by its respective developers. The data required to run each script is available on GEO, with the accession codes provided in the paper. If all packages are installed correctly, the scripts should run without any issues in RStudio or Terminal using the Rscript command. Each script is well-documented and clearly explains the purpose of each code block. The following table shows R and Python packages along with their version for this study.<br><br>

| R                        | R_version               | Python           | Py_version  |
|:--------------------------|:------------------------|:------------------|:------------|
| cowplot                 | 1.1.3                  | adjustText       | 1.3.0      |
| dbscan                  | 1.2.0                  | anndata          | 0.11.1     |
| doParallel              | 1.0.17                 | IPython          | 8.17.2     |
| dplyr                   | 1.1.4                  | kneed            | 0.8.5      |
| ggplot2                 | 3.5.1                  | matplotlib       | 3.9.2      |
| ggpubr                  | 0.6.0                  | matplotlib-inline | 0.1.6      |
| ggrepel                 | 0.9.6                  | numpy            | 1.26.2     |
| ggridges                | 0.5.6                  | pandas           | 2.2.3      |
| grid                    | 4.4.1                  | scanpy           | 1.10.4     |
| gridExtra               | 2.3                    | scipy            | 1.11.4     |
| hacksig                 | 0.1.2                  | seaborn          | 0.13.2     |
| infercnv                | 1.20.0                 | scikit-learn     | 1.5.2      |
| Matrix                  | 1.7.1                  | statsmodels      | 0.14.4     |
| MatrixGenerics          | 1.16.0                 |                  |            |
| miloR                   | 2.0.0                  |                  |            |
| monocle3                | 1.3.7                  |                  |            |
| parallel                | 4.4.1                  |                  |            |
| patchwork               | 1.3.0                  |                  |            |
| pheatmap                | 1.0.12                 |                  |            |
| rstatix                 | 0.7.2                  |                  |            |
| scater                  | 1.32.1                 |                  |            |
| Seurat                  | 5.2.0,4.1.4,4.4.0      |                  |            |
| SeuratDisk              | 0.0.0.9021             |                  |            |
| SeuratObject            | 5.0.2                  |                  |            |
| SeuratWrappers          | 0.3.5                  |                  |            |
| SingleCellExperiment    | 1.26.0                 |                  |            |
| statmod                 | 1.5.0                  |                  |            |
| survival                | 3.8.3                  |                  |            |
| survminer               | 0.5.0                  |                  |            |
| tidyr                   | 1.3.1                  |                  |            |
| viridis                 | 0.6.5                  |                  |            |
| xlsx                    | 0.6.5                  |                  |            |

<br>

**Instructions for Use**<br>

To run each figure script:
1. clone this GitHub repository  (https://github.com/IzarLab/PDAC_tropism.git) to a local directory (e.g., in Documents).
2. download **compartments.zip** and **Misc.zip** from Gene Expression Omnibus cited in the paper.
3. unzip them and move **Misc** and **compartments** to the directory (like **Documents/PDAC_tropism/**) containing **scripts**. This directory should now contain three folders: **scripts**, **Misc** and **compartments**.
4. open RStudio (or any other R-enabled IDE) and set **working directory**, `setwd()`, to **Documents/PDAC_tropism/**. Each **Fig_x** subfolder contains R scripts corresponding to figure panel described in the paper.

**Notes:**
1. For generating Fig. 2c, download **PN14_S134_L002_tumor.cna.seg** and **PN14_S134_L002_tumor.seg** available on GEO lpWGS repository (accession number given in the paper), then run `ichorCNA` with the parameters specified in **scripts/Fig_2/c.pdf**.
2. To generate the circos plot (Fig. 4j) there are two ways.
   - **_Easy way_**:  
     - open a Terminal session and change the directory to **scripts/Fig_4/**. Type in `jupyter notebook` and open **j.ipynb** in a web browser (e.g., Chrome).  
     - install the Python packages listed in **j.ipynb**, **j.py**, and also install _circos_ (https://circos.ca/).  
     - run cells in **j.ipynb**. This should create **output** subdirectory inside **Fig_4** with necessary files for generating circos plot, and possibly the plot itself. If the script throws _FileNotFoundError_, then open a new Terminal session, navigate to **output**, and type `circos` that should generate **circos.svg** in the **output** subfolder.  
   - **_Hard way_**: follow the instructions on the ContactTracing GitHub page (https://github.com/LaughneyLab/ContactTracing_tutorial), and replace their jupyter notebook (tutorial.ipynb) and python script (python/contactTracing_library.py) with **j.ipynb** and **j.py** respectively.
3. For running Echidna, please refer to their GitHub repository: https://github.com/azizilab/echidna

The following video demo shows how to reproduce Figure 1 as an example.<br><br>

[![Video Thumbnail](https://img.youtube.com/vi/zvmdHKROiBA/0.jpg)](https://www.youtube.com/watch?v=zvmdHKROiBA)
