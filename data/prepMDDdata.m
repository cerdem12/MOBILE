% Update the below filenames, can be .csv, .txt, etc.
% Given are examples from MDD datasets

%% Read-in RPPA lvl4 data
data_RPPA = readtable('MDD_RPPA_Level4.csv');
txt_RPPA = readtable('MDD_RPPA_antibodyAnnotations.csv');

% Reorder the MDD data into coherent sample-column-order (PBS-->TGFB)
colmorder_RPPA = [6,28,30,8,10,13,15,23,25,18,20,2,4,33,35]; 

RPPA_lvl4 = struct;
RPPA_lvl4.data = data_RPPA{:,colmorder_RPPA+1}; % Ordered data columns
RPPA_lvl4.samples = data_RPPA.Properties.VariableNames(colmorder_RPPA+1); % Ordered data column names (aka conditions)
RPPA_lvl4.antibody = data_RPPA{:,1};
RPPA_lvl4.antibodyfull = txt_RPPA{:,1};
RPPA_lvl4.hgnc_id = txt_RPPA{:,2};
RPPA_lvl4.p_sites = txt_RPPA{:,3};
RPPA_lvl4.effect = txt_RPPA{:,4};
RPPA_lvl4.pathway = txt_RPPA{:,5};

save('RPPA_lvl4.mat','RPPA_lvl4''-v7.3')


%% Read-in RNAseq lvl4 data
data_RNA = readtable('MDD_RNAseq_Level4.csv');
txt_RNA = readtable('MDD_RNAseq_geneAnnotations.csv');

% Reorder the MDD data into coherent sample-column-order (PBS-->TGFB)
colmorder_RNA = [1,2,3,14,15,10,11,12,13,6,7,4,5,8,9];

RNAseq_lvl4 = struct;
RNAseq_lvl4.data = data_RNA{:,colmorder_RNA+1}; % Ordered data columns
RNAseq_lvl4.samples = data_RNA.Properties.VariableNames(colmorder_RNA+1); % Ordered data columns
RNAseq_lvl4.ensembl_id = txt_RNA{:,1};
RNAseq_lvl4.hgnc_id = txt_RNA{:,2};

% Find nnz for each transcript and remove if not measured at least in 3 samples
numsamples = size(RNAseq_lvl4.data,2);
numZEROS_RNA = numsamples-(sum(RNAseq_lvl4.data~=0,2)); 
nnz3_RNA = find(numZEROS_RNA<=(numsamples-3)); % number of transcripts measured at least in 3 samples
NOT_nnz3_RNA = (setdiff([1:size(RNAseq_lvl4.data,1)],nnz3_RNA))';

% Create the new struct with low-measured transcripts removed
RNAseq_lvl42 = RNAseq_lvl4; 
RNAseq_lvl42.data(NOT_nnz3_RNA,:) = [];
RNAseq_lvl42.ensembl_id(NOT_nnz3_RNA,:) = [];
RNAseq_lvl42.hgnc_id(NOT_nnz3_RNA,:) = [];

save('RNAseq_lvl4.mat','RNAseq_lvl42''-v7.3')


%% Read-in ATACseq lvl4 data
clc
data_ATAC = readtable('MDD_ATACseq_Level4.csv');
opts = detectImportOptions("MDD_ATACseq_peakMetadata.csv",'NumHeaderLines',0,'VariableNamingRule','preserve');
txt_ATAC = readtable('MDD_ATACseq_peakMetadata.csv',opts);

% Reorder the MDD data into coherent sample-column-order (PBS-->TGFB)
colmorder_ATAC = [1,2,3,14,15,10,11,12,13,6,7,4,5,8,9];

ATACseq_lvl4 = struct;
ATACseq_lvl4.data = data_ATAC{:,colmorder_ATAC+1}; % Ordered data columns
ATACseq_lvl4.samples = data_ATAC.Properties.VariableNames(colmorder_ATAC+1); % Ordered data columns
ATACseq_lvl4.peak = data_ATAC{:,1};
ATACseq_lvl4.peakstart = num2cell(txt_ATAC{:,3});
ATACseq_lvl4.peakend = num2cell(txt_ATAC{:,4});
ATACseq_lvl4.annot = txt_ATAC{:,7};
ATACseq_lvl4.dist2TSS = num2cell(txt_ATAC{:,13});
ATACseq_lvl4.ensembl_id = txt_ATAC{:,15};
ATACseq_lvl4.entrez_id = num2cell(txt_ATAC{:,14});
ATACseq_lvl4.hgnc_id = txt_ATAC{:,16};

% Find and remove any peak regions with no ENS IDs
noENS_ATAC = find(strcmp('NA',ATACseq_lvl4.ensembl_id));
ATACseq_lvl42 = ATACseq_lvl4;
ATACseq_lvl42.data(noENS_ATAC,:) = [];
ATACseq_lvl42.peak(noENS_ATAC,:) = [];
ATACseq_lvl42.peakstart(noENS_ATAC,:) = [];
ATACseq_lvl42.peakend(noENS_ATAC,:) = [];
ATACseq_lvl42.annot(noENS_ATAC,:) = [];
ATACseq_lvl42.dist2TSS(noENS_ATAC,:) = [];
ATACseq_lvl42.ensembl_id(noENS_ATAC,:) = [];
ATACseq_lvl42.entrez_id(noENS_ATAC,:) = [];
ATACseq_lvl42.hgnc_id(noENS_ATAC,:) = [];

save('ATACseq_lvl4.mat','ATACseq_lvl42''-v7.3')

