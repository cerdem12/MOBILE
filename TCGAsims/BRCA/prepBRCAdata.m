%%%% Below script prepares TCGA data into a MOBILE friendly format.
%%%% Users need to have a "Sample Sheet" file, downloaded from the TCGA 
%%%% repository about the data files downloaded for the analysis.
%%%% 
%%%% One also needs a reference file to get gene identifiers (HGNC) for
%%%% the RPPA antibodies. Can be downloaded from:
%%%% https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
%%

%% Get all file/folder names of TCGA data
clc
namingsFile = 'TCGA_sample_sheet_file_name.txt';
files = tdfread(namingsFile);

%% Find case-matched RNAseq and RPPA data folder names
clc
disp('Start')
rnaFiles = {};
protFiles = {};
caseFiles = {};

Cases = cellstr(files.Case_ID);
DataTypes = cellstr(files.Data_Category);
FileNames = cellstr(files.File_Name);
SampleTypes = cellstr(files.Sample_Type);

qq = 1;
while qq <= size(Cases,1)
    tempC = Cases{qq,1};
    idxSame0 = find(strcmp(Cases,tempC));
    idxPriTum = [];
    for ss = 1:size(idxSame0,1)
        % Find-out if the case has Primary tumor data, will exclude Metastaic and solid normal tissue datasets here
        if (strcmp(SampleTypes{idxSame0(ss),1},"Primary Tumor")) 
            idxPriTum(ss,1) = idxSame0(ss);
        end
    end
    idxSame = intersect(idxSame0,idxPriTum);
    if size(idxSame,1)>1
        idProtPr = 0; % Find-out if the case has proteomics data
        for bb = 1:size(idxSame,1)
            if (convertCharsToStrings(DataTypes{idxSame(bb),1})=="Proteome Profiling") 
                idProtPr = bb;
            end
        end
        if idProtPr==1
            rnaFiles{end+1,1} = FileNames{idxSame(2),1};
            protFiles{end+1,1} = FileNames{idxSame(idProtPr),1};
            caseFiles{end+1,1} = Cases{idxSame(1),1};
        elseif idProtPr>1
            rnaFiles{end+1,1} = FileNames{idxSame(1),1}; % Even if there could be more than one RNAseq data, we store the first only
            protFiles{end+1,1} = FileNames{idxSame(idProtPr),1};
            caseFiles{end+1,1} = Cases{idxSame(1),1};
        end        
        qq = qq + size(idxSame0,1);
    else
        qq = qq + 1;
    end
end

save('TCGA_metadata.mat','caseFiles','Cases','DataTypes', ...
     'FileNames','files','namingsFile','protFiles','rnaFiles', ...
     'SampleTypes','-v7.3')

disp('End')

%% Read-in RNAseq fkpm and RPPA arrays and save separate files
clc
RNAdata_genesAll = {};
for qq = 1:size(caseFiles,1)
    rr = readtable(rnaFiles{qq,1},"FileType","text");
    geneIDs = rr.gene_id(5:end,1);
    genes = rr.gene_name(5:end,1);
    fpkm = rr.fpkm_unstranded(5:end,1);
    RNAdata_genesAll = union(RNAdata_genesAll,genes,'stable');
    save(['RNAdata_ens_' num2str(qq) '.mat'],'geneIDs','-v7.3')
    save(['RNAdata_genes_' num2str(qq) '.mat'],'genes','-v7.3')
    save(['RNAdata_fpkm_' num2str(qq) '.mat'],'fpkm','-v7.3')
    disp([1 qq])
end
save('RNAdata_genesAll.mat','RNAdata_genesAll','geneIDs','-v7.3')

%% RPPA
clc
RPPAdata_anstAll = {};
for qq = 1:size(caseFiles,1)
    pp = tdfread(protFiles{qq,1});
    agid = cellstr(pp.AGID);
    ants = cellstr(pp.peptide_target);
    expLevels = cellstr(pp.protein_expression);
    catNOs = cellstr(pp.catalog_number);
    RPPAdata_anstAll = union(RPPAdata_anstAll,ants,'stable');
    save(['RPPAdata_agid_' num2str(qq) '.mat'],'agid','-v7.3')
    save(['RPPAdata_ants_' num2str(qq) '.mat'],'ants','-v7.3')
    save(['RPPAdata_exps_' num2str(qq) '.mat'],'expLevels','-v7.3')
    save(['RPPAdata_cats_' num2str(qq) '.mat'],'catNOs','-v7.3')
    disp([2 qq])
end
save('RPPAdata_antsAll.mat','RPPAdata_anstAll','-v7.3')

%% Read-in RNAseq fkpm and RPPA data and merge them into matrices
clc
load('TCGA_metadata.mat')
load('RNAdata_genesAll.mat')
load('RNAdata_genes_1.mat')
genes0 = genes;
RNAseqData = zeros(size(genes0,1),size(caseFiles,1));
for qq = 1:size(caseFiles,1)
    load(['RNAdata_genes_' num2str(qq) '.mat'])
    load(['RNAdata_fpkm_' num2str(qq) '.mat'])
    if ~strcmp(genes0,genes)
        disp([9999 qq]);
    end
    RNAseqData(:,qq) = fpkm;
    disp([1 qq])
end
save('RNAseqData_fpkm.mat','RNAseqData','-v7.3')

%% Construct and save RPPA data matrix
clc
load('TCGA_metadata.mat')
load('RPPAdata_antsAll.mat')
load('RPPAdata_ants_1.mat')
ants0 = ants;
RPPAData = zeros(size(ants0,1),size(caseFiles,1));
for qq = 1:size(caseFiles,1)
    load(['RPPAdata_ants_' num2str(qq) '.mat'])
    load(['RPPAdata_exps_' num2str(qq) '.mat'])
    if (size(ants,1)~=size(ants0,1))
        [~,~,ib] = intersect(ants,ants0);
        RPPAData(ib,qq) = str2double(expLevels);
    else
        if (~strcmp(ants0,ants))
            disp([9999 qq]);
        else
            RPPAData(:,qq) = str2double(expLevels);
        end
    end
    disp([2 qq])
end
save('RPPAData_explvl.mat','RPPAData','-v7.3')

%% Find and remove any transcript with no HGNC or other annotaions
clc
ele2rmv1 = {'_','.','LINC','MIR','Telomerase-vert','hsa-mir-1253', ...
            'hsa-mir-423','snoZ196','7SK'};
ele2rmv2 = {'^U\d{1,7}$'};
ele2rmv3 = {'RN7SK','RNU','RPL','RPSAP'};
 
rows2rmv1 = {};
counter1 = 1;
for qq = ele2rmv1
    rows2rmv1{1,counter1} = find(~cellfun(@isempty,(strfind(genes0,qq))));
    counter1 = counter1 + 1;
end

rows2rmv2 = {};
counter1 = 1;
for qq = ele2rmv2
    rows2rmv2{1,counter1} = find(~cellfun(@isempty,regexp(genes0,qq)));
    counter1 = counter1 + 1;
end

rows2rmv3 = {};
counter1 = 1;
for qq = ele2rmv3
    rows2rmv3{1,counter1} = find(startsWith(genes0,qq));
    counter1 = counter1 + 1;
end
disp('Done')

% Concataneta all to be removed rows
row2rmAll12 = [];
row2rmAll123 = [];
for qq = rows2rmv1
    row2rmAll12 = [row2rmAll12;qq{1,1}];
    row2rmAll123 = [row2rmAll123;qq{1,1}];
end
for qq = rows2rmv2
    row2rmAll12 = [row2rmAll12;qq{1,1}];
    row2rmAll123 = [row2rmAll123;qq{1,1}];
end
for qq = rows2rmv3
    row2rmAll123 = [row2rmAll123;qq{1,1}];
end
row2rmAll12 = (sortrows(row2rmAll12));
row2rmAll123 = (sortrows(row2rmAll123));
disp('Done2')

%% Save the new RNAseq data after removed rows
RNAseqData1 = RNAseqData;
RNAseqData1(row2rmAll12,:) = [];
genes1 = genes0;
genes1(row2rmAll12,:) = [];
geneIDs1 = geneIDs;
geneIDs1(row2rmAll12,:) = [];
RNAseqData10 = RNAseqData1;

num0sallowed = 90; % ~10% of all columns (N=878)
numZEROS_RNA = size(RNAseqData1,2)-(sum(RNAseqData1~=0,2)); 
nnzs_RNA = find(numZEROS_RNA<(size(RNAseqData1,2)-num0sallowed+1)); 
NOT_nnzs_RNA = (setdiff([1:size(RNAseqData1,1)],nnzs_RNA))';

RNAseqData1(NOT_nnzs_RNA,:) = [];
genes1(NOT_nnzs_RNA,:) = [];
geneIDs1(NOT_nnzs_RNA,:) = [];

% Save the raw input data for MOBILE
TCGA_RNAseq_v1 = struct;
TCGA_RNAseq_v1.data = log2(RNAseqData1+1); % save as log2(fpkm+1)
TCGA_RNAseq_v1.hgnc = genes1;
TCGA_RNAseq_v1.ens = geneIDs1;
save('TCGA_RNAseq_v1','TCGA_RNAseq_v1','-v7.3')

%% Clean-up RPPA data - remove NaN rows and add HGNC and ENSEMBL IDs
antGenesFile = 'TCGA_antibodies_descriptions_gencode_v36.txt';
antInfo = tdfread(antGenesFile);
antGenes = cellstr(antInfo.gene_name);
antAGIDs = cellstr(antInfo.AGID);
antENSGs = cellstr(antInfo.gene_id);
hgncIDs = {};
ensgIDs = {};
for qq = 1:size(ants0,1)
    aa = find(strcmp(agid,antAGIDs{qq,1}));
    hgncIDs{qq,1} = antGenes{aa,1};
    ensgIDs{qq,1} = antENSGs{aa,1};
end

RPPAData1 = RPPAData;
ants1 = ants0;
idx2rmv = find(sum(isnan(RPPAData1),2));
RPPAData1(idx2rmv,:) = [];
ants1(idx2rmv,:) = [];
hgncIDs(idx2rmv,:) = [];
ensgIDs(idx2rmv,:) = [];

TCGA_RPPA_v1 = struct;
TCGA_RPPA_v1.data = RPPAData1;
TCGA_RPPA_v1.hgnc = hgncIDs;
TCGA_RPPA_v1.ens = ensgIDs;
TCGA_RPPA_v1.antibody = ants1;
save('TCGA_RPPA_v1','TCGA_RPPA_v1','-v7.3')

%% The rest are BRCA data specific, needs to be modified for others
clc
XXraw = TCGA_RNAseq_v1.data; % Right-hand-side (RHS) matrix: mRNAs
YYraw = TCGA_RPPA_v1.data; % Left-hand-side (LHS) matrix: Proteins
XXvar = flipud(sortrows([std(XXraw,1,2),(1:size(XXraw,1))'])); % calculate variance (std) of raw data
YYvar = flipud(sortrows([std(YYraw,1,2),(1:size(YYraw,1))']));
XXindices2keep = floor(1*size(XXraw,1)); % Find top 10% of highly variant mRNAs 
YYindices2keep = floor(1*size(YYraw,1)); % Find top 20% of highly variant proteins
XX = XXraw(sortrows(XXvar(1:XXindices2keep,2)),:); % Retain only top 10% of mRNAs 
YY = YYraw(sortrows(YYvar(1:YYindices2keep,2)),:); % Retain only top 20% of proteins/phosphoproteins
RNAseqIDs = sortrows(XXvar(1:XXindices2keep,2)); % Find indeces of retained mRNAs in raw data
RPPAIDs = sortrows(YYvar(1:YYindices2keep,2)); % Find indeces of retained proteins in raw data
RNAseqgenes = ;
RPPAgenes = ;

%% Find samples for different subtypes, based on RNA expression levels
genes1 = TCGA_RNAseq_v1.hgnc(:,1);

gene = 'ERBB2';
idx_rna = find(~cellfun(@isempty,(strfind(genes1(RNAseqIDs,1),gene))));
AA_H2lowIDs = find(XX(idx_rna,:)<7); % 7 is the break point in the histogram

gene = 'PGR';
idx_rna = find(~cellfun(@isempty,(strfind(genes1(RNAseqIDs,1),gene))));
AA_PRlowIDs = find(XX(idx_rna,:)<1); 

gene = 'ESR1';
idx_rna = find(~cellfun(@isempty,(strfind(genes1(RNAseqIDs,1),gene))));
AA_ERlowIDs = find(XX(idx_rna,:)<1); 

AA_TNBCids1 = intersect(AA_H2lowIDs,AA_PRlowIDs);
AA_TNBCids2 = intersect(AA_TNBCids1,AA_ERlowIDs);

%% Plot transcriptomic data
genes2 = genes1(RNAseqIDs,1);
figure;
subplot(1,3,1); hold on;
gene = 'ERBB2';
idx_rna = find(~cellfun(@isempty,(strfind(genes2,gene))));
histogram(XX(idx_rna(1),AA_TNBCids2))
ylabel('frequency')
xlabel('log2(fpkm+1)')
title(genes2(idx_rna(1)))
set(gca,'FontWeight','bold','FontSize',14,'FontName','Arial')

subplot(1,3,2); hold on;
gene = 'PGR';
idx_rna = find(~cellfun(@isempty,(strfind(genes2,gene))));
histogram(XX(idx_rna(1),AA_TNBCids2))
ylabel('frequency')
xlabel('log2(fpkm+1)')
title(genes2(idx_rna(1)))
set(gca,'FontWeight','bold','FontSize',14,'FontName','Arial')

subplot(1,3,3); hold on;
gene = 'ESR1';
idx_rna = find(~cellfun(@isempty,(strfind(genes2,gene))));
histogram(XX(idx_rna(1),AA_TNBCids2))
ylabel('frequency')
xlabel('log2(fpkm+1)')
title(genes2(idx_rna(1)))
set(gca,'FontWeight','bold','FontSize',14,'FontName','Arial')

%% Plot proteomics data
genes3 = TCGA_RPPA_v1.hgnc(RPPAIDs,1);

figure;
subplot(1,3,1); hold on;
gene = 'ERBB2';
idx_rna = find(~cellfun(@isempty,(strfind(genes3,gene))));
histogram(YY(idx_rna(1),:))
ylabel('frequency')
xlabel('arb. units')
title(genes3(idx_rna(1)))
set(gca,'FontWeight','bold','FontSize',14,'FontName','Arial')

subplot(1,3,2); hold on;
gene = 'PGR';
idx_rna = find(~cellfun(@isempty,(strfind(genes3,gene))));
histogram(YY(idx_rna(1),AA_TNBCids2))
ylabel('frequency')
xlabel('arb. units')
title(genes3(idx_rna(1)))
set(gca,'FontWeight','bold','FontSize',14,'FontName','Arial')

subplot(1,3,3); hold on;
gene = 'ESR1';
idx_rna = find(~cellfun(@isempty,(strfind(genes3,gene))));
histogram(YY(idx_rna(1),:))
ylabel('frequency')
xlabel('arb. units')
title(genes3(idx_rna(1)))
set(gca,'FontWeight','bold','FontSize',14,'FontName','Arial')











