close all
clear all
clc

addpath('/pfs/stor10/users/home/y/yujiao/MOBILE')
addpath('/pfs/stor10/users/home/y/yujiao/MOBILE/glmnet','/pfs/stor10/users/home/y/yujiao/MOBILE/funcs') % Make sure locations are correct
%% Data/matrix pre-processing
 % RPPA and RNAseq data
load('TCGA_RPPA_v1.mat')
load('TCGA_RNAseq_v1.mat')
rng(6); % Set random number generating (rng) seed
warning('off','all') % Turn-off warning messages, mostly from glmnet package

XXraw = TCGA_RNAseq_v1.data; % Right-hand-side (RHS) matrix: mRNAs
YYraw = TCGA_RPPA_v1.data; % Left-hand-side (LHS) matrix: Proteins
XXvar = flipud(sortrows([std(XXraw,1,2),(1:length(XXraw))'])); % calculate variance (std) of raw data
YYvar = flipud(sortrows([std(YYraw,1,2),(1:length(YYraw))']));
XXindices2keep = floor(0.1*size(XXraw,1)); % Find top 10% of highly variant mRNAs 
YYindices2keep = floor(0.2*size(YYraw,1)); % Find top 20% of highly variant proteins
XX = XXraw(sortrows(XXvar(1:XXindices2keep,2)),:); % Retain only top 10% of mRNAs 
YY = YYraw(sortrows(YYvar(1:YYindices2keep,2)),:); % Retain only top 20% of proteins/phosphoproteins
RNAseqIDs = sortrows(XXvar(1:XXindices2keep,2)); % Find indeces of retained mRNAs in raw data
RPPAIDs = sortrows(YYvar(1:YYindices2keep,2)); % Find indeces of retained proteins in raw data

genes1 = TCGA_RNAseq_v1.hgnc(:,1);
gene = 'ERG';
idx_rna = find(matches(genes1(RNAseqIDs,1),gene));
AA_ERGlowIDs = find(XX(idx_rna,:)<1.5); % 3 is the break point in the histogram
AA_ERGlowIDs = AA_ERGlowIDs(:);
gene = 'ETV1';
idx_rna = find(matches(genes1(RNAseqIDs,1),gene));
AA_ETV1lowIDs = find(XX(idx_rna,:)<1.5); % 3 is the break point in the histogram
AA_ETV1lowIDs = AA_ETV1lowIDs(:);
gene = 'ETV4';
idx_rna = find(matches(genes1(RNAseqIDs,1),gene));
AA_ETV4lowIDs = find(XX(idx_rna,:)<1.5); % 3 is the break point in the histogram
AA_ETV4lowIDs = AA_ETV4lowIDs(:);
gene = 'ETV5';
idx_rna = find(matches(genes1(RNAseqIDs,1),gene));
AA_ETV5lowIDs = find(XX(idx_rna,:)<1.5); % 3 is the break point in the histogram
AA_ETV5lowIDs = AA_ETV5lowIDs(:);
AA_ETVlowIDs1 = union(AA_ETV1lowIDs,AA_ETV4lowIDs);
AA_ETVlowIDs = union(AA_ETVlowIDs1,AA_ETV5lowIDs);
AA_ETVlowIDs = AA_ETVlowIDs(:);
gene = 'SPINK1';
idx_rna = find(matches(genes1(RNAseqIDs,1),gene));
AA_SPINK1lowIDs = find(XX(idx_rna,:)<1.5); % 3 is the break point in the histogram
AA_SPINK1lowIDs = AA_SPINK1lowIDs(:);
TNPCIDs1 =  intersect(AA_ERGlowIDs,AA_ETVlowIDs);
TNPCIDs =  intersect(TNPCIDs1,AA_SPINK1lowIDs);
%%%%% Remove TNPC and run: LOGO happens here
XX(:,TNPCIDs) = []; 
YY(:,TNPCIDs) = [];
% Normalize the input matrices 
XX3 = prepnormmats(XX,6,1); 
YY3 = prepnormmats(YY,6,1);
%% Simulation clean-up and save Lasso Coefficient Matrices
clc
addpath('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD/TNPC')
listTLasFiles = dir('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD/TNPC/TLas1_r*'); 
numFiles = size(listTLasFiles,1);
TLasMats = cell(numFiles,2);
qqtodel = [];
for qq = 1:numFiles % Loop for saving coefficient matrices
    if listTLasFiles(qq,1).bytes > 100000
        name1 = listTLasFiles(qq,1).name;
        lname1 = load(name1,'TLas1');
        TLasMats{qq,1} = sparse(lname1.TLas1);
        clear lname1
        aname1 = find(name1=='_');
        NumName1 = name1((aname1(end)+2):end-4);
        TLasMats{qq,2} = str2double(NumName1);
        disp(qq)
    else
        qqtodel = [qqtodel;qq]; % If for some reason there was an error in saving files
    end
end
a0Mats = cell(numFiles,2);
for qq = 1:numFiles % Loop for saving y-intercept value arrays (usually all zero)
    if listTLasFiles(qq,1).bytes > 100000
        name1 = listTLasFiles(qq,1).name;
        lname1 = load(name1,'a0Vals1');
        a0Mats{qq,1} = sparse(lname1.a0Vals1);
        clear lname1 
        aname1 = find(name1=='_');
        NumName1 = name1((aname1(end)+2):end-4);
        a0Mats{qq,2} = str2double(NumName1);
        disp(qq)
    end
end
listTLasFiles(qqtodel,:) = [];
a0Mats(qqtodel,:) = [];
TLasMats(qqtodel,:) = [];
save('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD/TNPC/PRAD_TNPC.mat','a0Mats','TLasMats','listTLasFiles','qqtodel','-v7.3') ;

%% Finding Robust Lasso Coefficient Matrix
clc
% Remember to run/load the "Data pre-processing" section above if needed (e.g. working on a different day after simulations are done)
load('PRAD_TNPC.mat') % Load the Lasso coefficient matrices data if not in workspace already

[intranksR_PRAD_TNPC,intsii_PRAD_TNPC,numedges_PRAD_TNPC,TLasTarget_PRAD_TNPC] = ...
    runAssocRankerTCGA(0.5,TLasMats,TCGA_RNAseq_v1,RNAseqIDs,TCGA_RPPA_v1,RPPAIDs);

% Change the below name according to inputs: For example, for EGF-LOGO -> TLasBestRpR_EGFLOGO
save('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD/TNPC/TLasBestPRAD_TNPC.mat','intranksR_PRAD_TNPC','intsii_PRAD_TNPC','numedges_PRAD_TNPC','TLasTarget_PRAD_TNPC','-v7.3') 
disp('DONE and saved')
