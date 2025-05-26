# Cytotoxic NK cells impede response to checkpoint immunotherapy in melanoma with an immune-excluded phenotype

This repository contains scripts used to analyse scRNAseq and spatial transcriptomic data in the doi: xxx publication.

## List of scripts / steps of the analysis

1. Selection of immune cells (GC_ALL_IMMUNE.R)
2. Reclustering of C8T and NK cells (CD8_and_NK_cells_subset.R)
3. Reclustering of Macrophages and DCs (Macro_DC_subset.R)
4. Reclustering of only Macrophages (Macrophages_subtypes.R)
5. Testing differences of immune cells abundances between responders and non-responders (Detailed_annotated_percentages.R)
6. Validation of NK cells in the Sade-Feldman et al. data (NK_validation_Sade_Feldman.R)
7. Validation of NK cells in the Riaz et al. data (NK_validation_in_Riaz_data.R)
9. scRNAseq analysis of mouse data (NRAS_YUMM_scRNAseq.R)
10. CITEseq analysis of mouse data (CITEseq.R)
11. Xenium data cell annotation (Xenium_data_annotation.ipynb)
12. Xenium data analysis after cell annotation (Xenium_data_after_annotation.ipynb)


*1 Preprocessing of all samples including DoubletFinder (GC_ALL.R) can be found here: https://github.com/MarineLab/Pozniak-et-al
*2 Inference of copy number variations (CNV; HB_GC_ALL.R) can be found here: https://github.com/MarineLab/Pozniak-et-al

## Location of the raw data

The human data was deposited to the EGA portal with the study number: EGAS00001006488
The mouse data was deposited to the GEO portal with the study number: GSE207592

### Location of the rds files of the entire tumor micoenvironment and the immune cells

https://doi.org/10.48804/GSAXBN
