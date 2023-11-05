# OmiHier: simultaneous learning of class label and hierarchy for omics data multi-class classification
*Jiemin Xie, Keyi Li, Zhanyu Liang, Bozhen Ren, Xuemei Liu, Yunhui Xiong, Li C. Xia**

## Summary
We proposed Omics Label Hierarchy Learning (OmiHier), a data-driven hierarchical structure learning algorithm. OmiHier adopts a bottom-up iterative framework, interlacing classification error minimization with successive label clustering, thus enables automatic and simultaneous learning of both class hierarchy and sample labels. We evaluated OmiHier on a number of simulated and real-world multi-omics datasets, including complex disease, microbiome, single-cell and spatial transcriptomics data. The benchmark demonstrated OmiHier’s high performance in classification accuracy, as well in inferring the true biological hierarchy.

The repository contains all the data (https://github.com/labxscut/OmiHier/main/Data) and code (https://github.com/labxscut/OmiHier/main/code.md) used in this study, as well as some important results (please refer to the figures(https://github.com/labxscut/OmiHier/main/Figures)).

# Data 

## Original data 

We downloaded six real-world multi-omics datasets for benchmark.

* ① Breast and gastric cancer data were downloaded from The Cancer Genome Atlas (TCGA) and the Molecular Taxonomy of Breast Cancer International Consortium (METABRIC), including mutation, copy number aberration and methylation data, and through both the cBioPortal (https://www.cbioportal.org/) and the TCGA data portal.

* ① Colon cancer microbiome data was downloaded from the ENA database (PRJEB7774).

* ① Gastric cancer cell line NCI-N87 scRNA-seq data was downloaded from Gene Expression Omnibus (GSE142750) and National Institute of Health’s SRA (PRJNA498809).

* ① Lymphoid cell scRNA-seq data was downloaded from human Ensemble Cell Atlas (hECA) system.

* ① Liver cancer spatial transcriptome data was downloaded from the HCCDB database (Integrative Molecular Database of Hepatocellular Carcinoma).

These datasets all come with sample- or cell-level multi-omics profiles and known biological labels, such as disease subtype, cell lineage and disease status.




# Contact & Support:

* Li C. Xia: email: [lcxia@scut.edu.cn](mailto:lcxia@scut.edu.cn)
* Jiemin Xie: email: [202120130808@mail.scut.edu.cn](mailto:202120130808@mail.scut.edu.cn)
