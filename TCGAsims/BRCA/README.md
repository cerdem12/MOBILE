# MOBILE on TCGA data

We show the applicability of MOBILE using tumor transcriptomic (RNAseq) and proteomic (RPPA) data from TCGA database. 

Preprocess TCGA downloaded data files using 'prepBRCAdata' script and start MOBILE pipeline with correctly formatted data matrices as input. 

In Breast Cancer (BRCA) analysis of 878 paired primary tumor cases with 27797 transcript and 457 protein levels, we determined three subtypes by looking at the transcript levels:

- HER2-amplified (HER2-amp)
- Triple-negative (TNBC)
- Estrogen and progesterone receptor positive (ER+/PR+)

Following MOBILE Lasso Module procedure and using all samples, FULL breast cancer integrated association network (FULL-TCGA-IAN) is generated. 
Then, by excluding subtype sample columns one at a time, we ran LOGO module to obtain three subtype-specific IANs.
