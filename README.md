# METI

## METI: Deep profiling of tumor ecosystems by integrating cell morphology and spatial transcriptomics

METI (Morphology-Enhanced Spatial Transcriptome Analysis Integrator) is an novel analytic framework that systematically analyzes cancer cells and cells of the TME by incorporating spatial gene expression, tissue histology, and prior knowledge of cancer and TME cells. METI starts with the identification of key cellular components and their states within the TME, including various immune cells and their transcriptional states, tumor stromal components such as cancer-associated fibroblasts (CAFs), and the epithelial compartment. Meanwhile, METI offers complementary information on cell morphology for various cell type from the H&E images. The combined results from gene expression and histology features provide a comprehensive understanding of the spatial cellular composition and organization within the tissue. 
<br>

![METI workflow](doc/workflow.png)

<br>

## Usage

With [**METI**](https://github.com/Flashiness/METI) package, you can:
- Map normal and premalignant cells;
- Identify cancer cell domain and heterogeneity;
- Map and phenotype T cell;
- Analyze other immune cells including B cells, plasma cells, neutrophils;
- Analyze stromal cells.
<br>
For tutorial, please refer to: https://github.com/Flashiness/METI/tree/main/tutorial.md
<br>

## System Requirements
Python support packages: torch, pandas, numpy, scipy, scanpy > 1.5, anndata, sklearn, cv2, TESLAforST.

## Versions the software has been tested on
- System: Anaconda
- Python: 3.7.9
- Python packages: pandas = 1.1.3, numpy = 1.20.2, python-igraph=0.8.3, torch=1.6.0,louvain=0.7.0, scipy = 1.5.2, scanpy = 1.6.0, anndata = 0.7.4,  sklearn = 0.23.3, cv2=4.5.1, TESLAforST=1.2.4.
<br>

## Contributing

Source code: [Github](https://github.com/Flashiness/METI)  

We are continuing adding new features. Bug reports or feature requests are welcome. 

## References

Please consider citing the following reference:

- https://doi.org/10.1101/2023.10.06.561287
<br>




