# MOBILE on PRAD data

Preprocess TCGA downloaded data files using 'prepPRADdata' script and start MOBILE pipeline with correctly formatted data matrices as input. 

In PRAD analysis of XXX paired primary tumor cases with XXX transcript and XXX protein levels, we determined four subtypes by looking at the transcript levels:

- ERG-positive
- ETS-positive
- SPINK1-positive
- Triple-negative (TNPC)

Following MOBILE Lasso Module procedure and using all samples, FULL prostate cancer integrated association network (FULL-IAN) is generated. 
Then, by excluding subtype sample columns one at a time, we ran LOGO module to obtain four subtype-specific IANs.
