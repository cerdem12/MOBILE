%%%%% MOBILE -omics integrator (Code by Cemal Erdem, Ph.D.)
%%%%% Last update: June 2022
%%%%% Integration of transcriptomic and epigenomic datasets 
close all
clear all
clc

addpath('glmnet','data','funcs')

%% Data/matrix pre-processing

load('ATACseq_RNAseq_lvl4_data')
rng(6);
warning('off','all')

XXraw = ATACseq_lvl42.data;
YYraw = RNAseq_lvl42.data;
XXvar = flipud(sortrows([std(XXraw,1,2),(1:length(XXraw))']));
YYvar = flipud(sortrows([std(YYraw,1,2),(1:length(YYraw))']));
XXindices2keep = floor(0.1*size(XXraw,1));
YYindices2keep = floor(0.1*size(YYraw,1));
XX = XXraw(sortrows(XXvar(1:XXindices2keep,2)),:);
YY = YYraw(sortrows(YYvar(1:YYindices2keep,2)),:);
ATACseqIDs = sortrows(XXvar(1:XXindices2keep,2));
RNAseqIDs = sortrows(YYvar(1:YYindices2keep,2));
% Remove 24hr and 48hr data here: LOGO happens here
% XX = XX(:,[1:3,6:15]); % For instance, excluding columns 4 & 5 makes EGF-LOGO
% YY = YY(:,[1:3,6:15]);
% Shuffle here OR not: Shuffle input data matrices as a control condition 
% XX = reshape(XX(randperm(numel(XX))),size(XX,1),size(XX,2));
% YY = reshape(YY(randperm(numel(YY))),size(YY,1),size(YY,2));
XX3 = prepnormmats(XX,6,1);
YY3 = prepnormmats(YY,6,1);

%% Simulations
clc
toplot = 0;
rtimes = 2; Number of Lasso modules to run. Use 2 to verify installation, set to 10000 for MOBILE
CVnfolds = 4;

mkdir FULL_RA % Create a new directory for outputs. Naming according to input data is advised
addpath('FULL_RA') % Add the new directory to the PATH

parfor rr = 1:rtimes
    rndgntr = randi(1e8);
    rng(rndgntr)
    % Run the glmnet package
    [TLas1,BLas1,a0Vals1,YYpg,YYpgdev] = glmnetRunner(XX3,YY3,toplot,CVnfolds);
    % Save the Lasso coefficients matrix and other info for each iteration
    parsavef(['FULL_RA/TLas_r' num2str(rr)],TLas1,BLas1,a0Vals1,YYpg,YYpgdev)
    % Clear large parameters from the workspace for faster simulations
    [TLas1,BLas1,a0Vals1,YYpg,YYpgdev] = parclearf(TLas1,BLas1,a0Vals1,YYpg,YYpgdev)
end

%% Simulation clean-up and save Lasso Coefficient Matrices
clc
listTLasFiles = dir('FULL_RA/TLas_r*'); 
numFiles = size(listTLasFiles,1);
TLasMats = cell(numFiles,2);
qqtodel = [];
for qq = 1:numFiles % Loop for saving coefficient matrices
    if listTLasFiles(qq,1).bytes > 1000000
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
    if listTLasFiles(qq,1).bytes > 1000000
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
save('FULL_RA/RA_FULL.mat','a0Mats','TLasMats','listTLasFiles','qqtodel','-v7.3') ;

%% Finding Robust Lasso Coefficient Matrix
clc
% Remember to run/load the "Data pre-processing" section above if needed (e.g. working on a different day after simulations are done)
% load('RA_FULL.mat') % Load the Lasso coefficient matrices data if not in workspace already

[intranksR_RA_FULL,intsii_RA_FULL,numedges_RA_FULL,TLasTarget_RA_FULL] = ...
    runAssocRankerRA(0.5,TLasMats,ATACseq_lvl42,ATACseqIDs,RNAseq_lvl42,RNAseqIDs);

% Change the below name according to inputs: For example, for EGF-LOGO -> TLasBestRA_EGFLOGO
save('FULL_RA/TLasBestRA_FULL.mat','intranksR_RA_FULL','intsii_RA_FULL','numedges_RA_FULL','TLasTarget_RA_FULL','-v7.3') 
disp('DONE and saved')
































