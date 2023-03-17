%%%%% MOBILE -omics integrator (Code by Cemal Erdem, Ph.D.)
%%%%% Last update: June 2022
%%%%% Integration of proteomic and transcriptomic datasets
addpath('glmnetNew','funcs')
% Data/matrix pre-processing
% RPPA and RNAseq data
clc
load('TCGA_RPPA_v1.mat')
load('TCGA_RNAseq_v1.mat')
rng(6); % Set random number generating (rng) seed
warning('off','all') % Turn-off warning messages, mostly from glmnet package

XXraw = TCGA_RNAseq_v1.data; % Right-hand-side (RHS) matrix: mRNAs
YYraw = TCGA_RPPA_v1.data; % Left-hand-side (LHS) matrix: Proteins
XXvar = flipud(sortrows([std(XXraw,1,2),(1:size(XXraw,1))'])); % calculate variance (std) of raw data
YYvar = flipud(sortrows([std(YYraw,1,2),(1:size(YYraw,1))']));
XXindices2keep = floor(0.1*size(XXraw,1)); % Find top 10% of highly variant mRNAs 
YYindices2keep = floor(0.2*size(YYraw,1)); % Find top 20% of highly variant proteins
XX = XXraw(sortrows(XXvar(1:XXindices2keep,2)),:); % Retain only top 10% of mRNAs 
YY = YYraw(sortrows(YYvar(1:YYindices2keep,2)),:); % Retain only top 20% of proteins/phosphoproteins
RNAseqIDs = sortrows(XXvar(1:XXindices2keep,2)); % Find indeces of retained mRNAs in raw data
RPPAIDs = sortrows(YYvar(1:YYindices2keep,2)); % Find indeces of retained proteins in raw data

%%%%% Find HER2 amplified samples ( log2(fpkm+1)>7 )
% gene = 'ERBB2';
% idx_rna = find(~cellfun(@isempty,(strfind(genes1(RNAseqIDs,1),gene))));
% H2ampIDs = find(XX(idx_rna,:)>7); % 7 is the break point in the histogram
%%%%% Remove HER2-amplified OR HER2-normal samples and run: LOGO happens here
% XX(:,H2ampIDs) = []; % For instance, excluding columns 4 & 5 makes EGF-LOGO
% YY(:,H2ampIDs) = [];

% Normalize the input matrices 
XX3 = prepnormmats(XX,6,1); 
YY3 = prepnormmats(YY,6,1);

% Simulations
toplot = 0; % flag to plot real YY data vs calculated YY
rtimes = 2000; % Number of Lasso modules to run
CVnfolds = 4; % Cross-validation number for glmnet package
% mkdir FULL_TCGA % Create a new directory for outputs. Naming according to input data is advised
addpath('FULL_TCGA') % Add the new directory to the PATH
tic
parfor rr = 1:rtimes % Run on parallel
    rndgntr = randi(1e8); % Set the rng seed for each iteration 
    rng(rndgntr)
    % Run the glmnet package
    [TLas1,BLas1,a0Vals1,YYpg,YYpgdev] = glmnetRunner(XX3,YY3,toplot,CVnfolds);
    % Save the Lasso coefficients matrix and other info for each iteration
    parsavef(['FULL_TCGA/TLas1_r' num2str(rr)],TLas1,BLas1,a0Vals1,YYpg,YYpgdev); % See the directory naming above
    % Clear large parameters from the workspace for faster simulations
    [TLas1,BLas1,a0Vals1,YYpg,YYpgdev] = parclearf(TLas1,BLas1,a0Vals1,YYpg,YYpgdev);
end
toc