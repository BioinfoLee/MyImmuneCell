# Human Circulating Neutrophil Dynamics Across the Lifespan
Large-scale single-cell profiling of human circulating neutrophils across the lifespan addresses a major gap in the understanding of immune homeostasis and aging. Here, My Immune Cell Landscape (MICL) is presented as a multimodal atlas of peripheral blood from 954 healthy Chinese individuals (2–100 years of age), integrating single-cell RNA, surface epitope, TCR, and BCR profiles from the same cells, along with bulk proteomics, whole-genome sequencing, and clinical laboratory measurements. The atlas resolves 99 biologically fine-grained cell subsets, comprising ~28.35 million cells from the granulocyte layer and ~23.54 million cells from the mononuclear cell layer.
This repository stores the steps used to integrates granulocytes into the circulating immune landscape based on multi-omic profiling. These notebooks can be run sequentially by Python or R kernel, as both languages are utilized in this analysis.


Stage1 Reference atlas：Prior to conducting the large-scale cohort study, we performed scRNA-seq, paired scTCR-seq and scBCR-seq of PBMCs, whole blood, and matched hybrid samples (comprising equal numbers of PBMCs and whole blood cells) from eight donors. This analysis was performed to characterize the cell-type composition of different samples and to establish a reference framework for subsequent cohort analyses.

Stage2 MyImmuneCell：Based on results of pilot experiments, hybrid samples (comprising equal cell numbers of PBMCs and whole-blood cells) were used for sequencing. For each sample, two technical replicates were set up to address potential technical variations and enable detection of larger number of cells per donor, followed by integrated cellular indexing of transcriptomes and epitopes (CITE-seq), scTCR-seq, and scBCR-seq. The code covers data preprocessing, cell annotation, and downstream analysis.

Legal Information

License
The license for this package is available on Github in the file LICENSE.txt in this repository.
