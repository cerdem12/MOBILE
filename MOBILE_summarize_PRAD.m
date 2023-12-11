%%%% This script creates a comprehensive list of all associations inferred
%%%% for FULL and all LOGO cases of TCGA_PRAD, including RNAseq-RPPA 

%% Here load all the data files FULL and LOGO (i.e. TLasBestRA_FULL.mat)
clc
rng(6);
warning('off','all')

% addpath('/pfs/stor10/users/home/y/yujiao/MOBILE')
% addpath('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD')
% addpath('/pfs/stor10/users/home/y/yujiao/MOBILE/glmnet','/pfs/stor10/users/home/y/yujiao/MOBILE/funcs')
% addpath('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD/TNPC')
% addpath('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD/SPINK1')
% addpath('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD/ERG')
% addpath('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD/ETS')
% addpath('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD/FULL')

load('TCGA_RPPA_v1.mat')
load('TCGA_RNAseq_v1.mat')
load('TLasBestPRAD_ERG.mat')
load('TLasBestPRAD_FULL.mat')
load('TLasBestPRAD_ETS.mat')
load('TLasBestPRAD_SPINK1.mat')
load('TLasBestPRAD_TNPC.mat')

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
% Normalize the input matrices 
XX3 = prepnormmats(XX,6,1); 
YY3 = prepnormmats(YY,6,1);


%%
% First,get the "FULL data" interactions matrix, will be a 35 column cell-array
clc 
TMat1 = TLasTarget_PRAD_FULL;
TMat2 = TLasTarget_PRAD_FULL;

[ZZ_FULL_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TMat1,TMat2,TCGA_RNAseq_v1,RNAseqIDs,TCGA_RPPA_v1,RPPAIDs,0);

ZZ_FULL_RpR(:,21) = num2cell(ones(size(ZZ_FULL_RpR,1),1));
ZZ_FULL_RpR(:,22:35) = num2cell(zeros(size(ZZ_FULL_RpR,1),14));


disp('DONE1')

%%
% Second,check each LeGO data matrices against the FULL case, and add the 
% interaction as a row if it is new, add info if it already exists.

%% ERG
[ZZ_test_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TLasTarget_PRAD_ERG,TLasTarget_PRAD_ERG,TCGA_RNAseq_v1,RNAseqIDs,TCGA_RPPA_v1,RPPAIDs,0);

testlen1 = size(ZZ_test_RpR,1);

for qq = 1:testlen1
    stest1 = strcmp(ZZ_test_RpR(qq,4),ZZ_FULL_RpR(:,4));
    stest2 = strcmp(ZZ_test_RpR(qq,7),ZZ_FULL_RpR(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RpR(ss2,22) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(ss2,23) = num2cell(1);
    else
        ZZ_FULL_RpR(end+1,1:20) = ZZ_test_RpR(qq,:);
        ZZ_FULL_RpR(end,22) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(end,23) = num2cell(1);
    end
end

disp('DONE2_1')

%% ETS
[ZZ_test_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TLasTarget_PRAD_ETS,TLasTarget_PRAD_ETS,TCGA_RNAseq_v1,RNAseqIDs,TCGA_RPPA_v1,RPPAIDs,0);

testlen1 = size(ZZ_test_RpR,1);

for qq = 1:testlen1
    stest1 = strcmp(ZZ_test_RpR(qq,4),ZZ_FULL_RpR(:,4));
    stest2 = strcmp(ZZ_test_RpR(qq,7),ZZ_FULL_RpR(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RpR(ss2,24) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(ss2,25) = num2cell(1);
    else
        ZZ_FULL_RpR(end+1,1:20) = ZZ_test_RpR(qq,:);
        ZZ_FULL_RpR(end,24) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(end,25) = num2cell(1);
    end
end

disp('DONE2_2')

%% SPINK1
[ZZ_test_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TLasTarget_PRAD_SPINK1,TLasTarget_PRAD_SPINK1,TCGA_RNAseq_v1,RNAseqIDs,TCGA_RPPA_v1,RPPAIDs,0);

testlen1 = size(ZZ_test_RpR,1);

for qq = 1:testlen1
    stest1 = strcmp(ZZ_test_RpR(qq,4),ZZ_FULL_RpR(:,4));
    stest2 = strcmp(ZZ_test_RpR(qq,7),ZZ_FULL_RpR(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RpR(ss2,26) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(ss2,27) = num2cell(1);
    else
        ZZ_FULL_RpR(end+1,1:20) = ZZ_test_RpR(qq,:);
        ZZ_FULL_RpR(end,26) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(end,27) = num2cell(1);
    end
end

disp('DONE2_3')

%% TNPC
[ZZ_test_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TLasTarget_PRAD_TNPC,TLasTarget_PRAD_TNPC,TCGA_RNAseq_v1,RNAseqIDs,TCGA_RPPA_v1,RPPAIDs,0);

testlen1 = size(ZZ_test_RpR,1);

for qq = 1:testlen1
    stest1 = strcmp(ZZ_test_RpR(qq,4),ZZ_FULL_RpR(:,4));
    stest2 = strcmp(ZZ_test_RpR(qq,7),ZZ_FULL_RpR(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RpR(ss2,28) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(ss2,29) = num2cell(1);
    else
        ZZ_FULL_RpR(end+1,1:20) = ZZ_test_RpR(qq,:);
        ZZ_FULL_RpR(end,28) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(end,29) = num2cell(1);
    end
end

disp('DONE2_4')
disp('DONE2_All')

%% Check the core of cores network - RpR
clc
[aa,bb] = find(cellfun('isempty', ZZ_FULL_RpR(:,21:end)));
ind = sub2ind(size(ZZ_FULL_RpR), aa, bb+20);
ZZ_FULL_RpR(ind) = num2cell(0);

coreIntsRpR = [(1:size(ZZ_FULL_RpR,1))', sum(cell2mat(ZZ_FULL_RpR(:,[21,23,25,27,29])),2)];
coreIntsRpR2 = (sortrows(coreIntsRpR,2));
coreIntsRpRList = ZZ_FULL_RpR(coreIntsRpR2(find(coreIntsRpR2(:,2)==5)),:);
disp('DONE3')

%% Rearrange final lists
ZZ_RpRList = {};
ZZ_RpRList(:,1:4) = ZZ_FULL_RpR(:,1:4);
ZZ_RpRList(:,5) = ZZ_FULL_RpR(:,16);
ZZ_RpRList(:,6:8) = ZZ_FULL_RpR(:,5:7);
ZZ_RpRList(:,9) = ZZ_FULL_RpR(:,17);
ZZ_RpRList(:,10) = ZZ_FULL_RpR(:,20);
ZZ_RpRList(:,11) = num2cell(coreIntsRpR(:,2));
ZZ_RpRList(:,12:16) = ZZ_FULL_RpR(:,21:2:29);
ZZ_RpRList(:,17:18) = [ZZ_FULL_RpR(:,8) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,8))))];
ZZ_RpRList(:,19:20) = [ZZ_FULL_RpR(:,22) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,22))))];
ZZ_RpRList(:,21:22) = [ZZ_FULL_RpR(:,24) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,24))))];
ZZ_RpRList(:,23:24) = [ZZ_FULL_RpR(:,26) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,26))))];
ZZ_RpRList(:,25:26) = [ZZ_FULL_RpR(:,28) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,28))))];


disp('DONE5')

%% Write the lists into files
clc
% Add Java POI Libs to matlab javapath(Maybe only necessary for Mac users)
javaaddpath(which('poi_library/poi-3.8-20120326.jar'));
javaaddpath(which('poi_library/poi-ooxml-3.8-20120326.jar'));
javaaddpath(which('poi_library/poi-ooxml-schemas-3.8-20120326.jar'));
javaaddpath(which('poi_library/xmlbeans-2.3.0.jar'));
javaaddpath(which('poi_library/dom4j-1.6.1.jar'));
javaaddpath(which('poi_library/stax-api-1.0.1.jar'));
%addpath('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD')
%directory = '/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD';
submodelfile = 'coreIntsRpRList.xlsx' ;
%fullFilePath = fullfile(directory, submodelfile);
xlwrite(submodelfile, coreIntsRpRList)

%directory = '/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD';
submodelfile = 'ZZ_RpRList.xlsx' ; 
%fullFilePath = fullfile(directory, submodelfile);
xlwrite(submodelfile, ZZ_RpRList)

disp('DONE6')

%% Save the lists
clc
%addpath('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD')
%save('/pfs/proj/nobackup/fs/projnb10/ecsbl/MOBILE/TCGA/PRAD/TLasBestRpRA_Summ.mat','ZZ_RpRList','-v7.3') 
save('/Users/root1/Downloads/PRAD_sum/TLasBestRpRA_Summ.mat','ZZ_RpRList','-v7.3') 

disp('DONE7')
