# MOBILE on TCGA data

We show the applicability of MOBILE using tumor transcriptomic (RNAseq) and proteomic (RPPA) data from TCGA database. 

Preprocess TCGA downloaded data files using 'prepdata' scripts and start MOBILE pipeline with correctly formatted data matrices as input. 

Examplary files for different cancer types are located under corresponding folders.

Following MOBILE Lasso Module procedure using all samples, integrated association networks (FULL-IANs) are generated. Then, by excluding subtype sample columns one at a time, run LOGO module to obtain subtype-specific IANs.
