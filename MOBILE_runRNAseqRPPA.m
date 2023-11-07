%%%%% MOBILE -omics integrator (Code by Cemal Erdem, Ph.D.)
%%%%% Last update: June 2022
%%%%% Integration of proteomic and transcriptomic datasets
close all
clear all
clc

addpath('glmnet','data','funcs')

%% Data/matrix pre-processing

load('RNAseq_RPPA_lvl4_data') % RPPA and RNAseq data
rng(6); % Set random number generating (rng) seed
warning('off','all') % Turn-off warning messages, mostly from glmnet package

XXraw = RNAseq_lvl42.data; % Right-hand-side (RHS) matrix: mRNAs
YYraw = RPPA_lvl4.data; % Left-hand-side (LHS) matrix: Proteins
XXvar = flipud(sortrows([std(XXraw,1,2),(1:length(XXraw))'])); % calculate variance (std) of raw data
YYvar = flipud(sortrows([std(YYraw,1,2),(1:length(YYraw))']));
XXindices2keep = floor(0.1*size(XXraw,1)); % Find top 10% of highly variant mRNAs 
YYindices2keep = floor(0.2*size(YYraw,1)); % Find top 20% of highly variant proteins
XX = XXraw(sortrows(XXvar(1:XXindices2keep,2)),:); % Retain only top 10% of mRNAs 
YY = YYraw(sortrows(YYvar(1:YYindices2keep,2)),:); % Retain only top 20% of proteins/phosphoproteins
RNAseqIDs = sortrows(XXvar(1:XXindices2keep,2)); % Find indeces of retained mRNAs in raw data
RPPAIDs = sortrows(YYvar(1:YYindices2keep,2)); % Find indeces of retained proteins in raw data
% Remove 24hr and 48hr data here: LOGO happens here
% XX = XX(:,[1:3,6:15]); % For instance, excluding columns 4 & 5 makes EGF-LOGO
% YY = YY(:,[1:3,6:15]);
% Shuffle here OR not: Shuffle input data matrices as a control condition 
% XX = reshape(XX(randperm(numel(XX))),size(XX,1),size(XX,2));
% YY = reshape(YY(randperm(numel(YY))),size(YY,1),size(YY,2));
% Normalize the input matrices 
XX3 = prepnormmats(XX,6,1); 
YY3 = prepnormmats(YY,6,1);

%% Simulations
clc
toplot = 0; % flag to plot real YY data vs calculated YY
rtimes = 10000; % Number of Lasso modules to run
CVnfolds = 4; % Cross-validation number for glmnet package

mkdir FULL_RpR % Create a new directory for outputs. Naming according to input data is advised
addpath('FULL_RpR') % Add the new directory to the PATH

parfor rr = 1:rtimes % Run on parallel
    rndgntr = randi(1e8); % Set the rng seed for each iteration 
    rng(rndgntr)
    % Run the glmnet package
    [TLas1,BLas1,a0Vals1,YYpg,YYpgdev] = glmnetRunner(XX3,YY3,toplot,CVnfolds);
    % Save the Lasso coefficients matrix and other info for each iteration
    parsavef(['FULL_RpR/TLas_r' num2str(rr)],TLas1,BLas1,a0Vals1,YYpg,YYpgdev) % See the directory naming above
    % Clear large parameters from the workspace for faster simulations
    [TLas1,BLas1,a0Vals1,YYpg,YYpgdev] = parclearf(TLas1,BLas1,a0Vals1,YYpg,YYpgdev)
end

%% Simulation clean-up and save Lasso Coefficient Matrices
clc
listTLasFiles = dir('FULL_RpR/TLas_r*'); 
numFiles = size(listTLasFiles,1);
TLasMats = cell(numFiles,2);
qqtodel = [];
for qq = 1:numFiles % Loop for saving coefficient matrices
    if listTLasFiles(qq,1).bytes > 30000
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
    if listTLasFiles(qq,1).bytes > 30000
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
save('FULL_RpR/RpR_FULL.mat','a0Mats','TLasMats','listTLasFiles','qqtodel','-v7.3') ;

%% Finding Robust Lasso Coefficient Matrix
clc
% Remember to run/load the "Data pre-processing" section above if needed (e.g. working on a different day after simulations are done)
% load('RpR_FULL.mat') % Load the Lasso coefficient matrices data if not in workspace already

[intranksR_RpR_FULL,intsii_RpR_FULL,numedges_RpR_FULL,TLasTarget_RpR_FULL] = ...
    runAssocRankerRpR(0.5,TLasMats,RNAseq_lvl42,RNAseqIDs,RPPA_lvl4,RPPAIDs);

% Change the below name according to inputs: For example, for EGF-LOGO -> TLasBestRpR_EGFLOGO
save('FULL_RpR/TLasBestRpR_FULL.mat','intranksR_RpR_FULL','intsii_RpR_FULL','numedges_RpR_FULL','TLasTarget_RpR_FULL','-v7.3') 
disp('DONE and saved')


