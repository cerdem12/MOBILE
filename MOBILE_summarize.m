%%%% This script creates a comprehensive list of all associations inferred
%%%% for FULL and all LOGO cases, including RNAseq-ATACseq & RNAseq-RPPA 


%% Here load all the data files FULL and LOGO (i.e. TLasBestRA_FULL.mat)
clc
rng(6);
warning('off','all')
load('RNAseq_RPPA_lvl4_data')
load('ATACseq_RNAseq_lvl4_data')
load('TLasBestRA_FULL.mat')
load('TLasBestRpR_FULL.mat')

load('The rest of files')

%%
% First,get the "FULL data" interactions matrix, will be a 35 column cell-array
clc 
TMat1 = TLasTarget_RpR_FULL;
TMat2 = TLasTarget_RpR_FULL;
TMat1ra = TLasTarget_RA_FULL;
TMat2ra = TLasTarget_RA_FULL;

[ZZ_FULL_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TMat1,TMat2,RNAseq_lvl42,RNAseqIDs,RPPA_lvl4,RPPAIDs,0);
[ZZ_FULL_RA,~,~,~] = ...
    TTcolorInts_RA_lvl4_v2(TMat1ra,TMat2ra,ATACseq_lvl42,ATACseqIDs,RNAseq_lvl42,RNAseqIDs,0);

ZZ_FULL_RpR(:,21) = num2cell(ones(size(ZZ_FULL_RpR,1),1));
ZZ_FULL_RpR(:,22:35) = num2cell(zeros(size(ZZ_FULL_RpR,1),14));

ZZ_FULL_RA(:,21) = num2cell(ones(size(ZZ_FULL_RA,1),1));
ZZ_FULL_RA(:,22:35) = num2cell(zeros(size(ZZ_FULL_RA,1),14));
disp('DONE1')

%%
% Second,check each LeGO data matrices against the FULL case, and add the 
% interaction as a row if it is new, add info if it already exists.

%% PBS
[ZZ_test_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TLasTarget_RpR_PBSLeGO,TLasTarget_RpR_PBSLeGO,RNAseq_lvl42,RNAseqIDs,RPPA_lvl4,RPPAIDs,0);
[ZZ_test_RA,~,~,~] = ...
    TTcolorInts_RA_lvl4_v2(TLasTarget_RA_PBSLeGO,TLasTarget_RA_PBSLeGO,ATACseq_lvl42,ATACseqIDs,RNAseq_lvl42,RNAseqIDs,0);
testlen1 = size(ZZ_test_RpR,1);
testlen2 = size(ZZ_test_RA,1);
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
for qq = 1:testlen2
    stest1 = strcmp(ZZ_test_RA(qq,4),ZZ_FULL_RA(:,4));
    stest2 = strcmp(ZZ_test_RA(qq,7),ZZ_FULL_RA(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RA(ss2,22) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(ss2,23) = num2cell(1);
    else
        ZZ_FULL_RA(end+1,1:20) = ZZ_test_RA(qq,:);
        ZZ_FULL_RA(end,22) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(end,23) = num2cell(1);
    end
end
disp('DONE2_1')

% EGF
[ZZ_test_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TLasTarget_RpR_EGFLeGO,TLasTarget_RpR_EGFLeGO,RNAseq_lvl42,RNAseqIDs,RPPA_lvl4,RPPAIDs,0);
[ZZ_test_RA,~,~,~] = ...
    TTcolorInts_RA_lvl4_v2(TLasTarget_RA_EGFLeGO,TLasTarget_RA_EGFLeGO,ATACseq_lvl42,ATACseqIDs,RNAseq_lvl42,RNAseqIDs,0);
testlen1 = size(ZZ_test_RpR,1);
testlen2 = size(ZZ_test_RA,1);
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
for qq = 1:testlen2
    stest1 = strcmp(ZZ_test_RA(qq,4),ZZ_FULL_RA(:,4));
    stest2 = strcmp(ZZ_test_RA(qq,7),ZZ_FULL_RA(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RA(ss2,24) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(ss2,25) = num2cell(1);
    else
        ZZ_FULL_RA(end+1,1:20) = ZZ_test_RA(qq,:);
        ZZ_FULL_RA(end,24) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(end,25) = num2cell(1);
    end
end
disp('DONE2_2')

% HGF
[ZZ_test_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TLasTarget_RpR_HGFLeGO,TLasTarget_RpR_HGFLeGO,RNAseq_lvl42,RNAseqIDs,RPPA_lvl4,RPPAIDs,0);
[ZZ_test_RA,~,~,~] = ...
    TTcolorInts_RA_lvl4_v2(TLasTarget_RA_HGFLeGO,TLasTarget_RA_HGFLeGO,ATACseq_lvl42,ATACseqIDs,RNAseq_lvl42,RNAseqIDs,0);
testlen1 = size(ZZ_test_RpR,1);
testlen2 = size(ZZ_test_RA,1);
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
for qq = 1:testlen2
    stest1 = strcmp(ZZ_test_RA(qq,4),ZZ_FULL_RA(:,4));
    stest2 = strcmp(ZZ_test_RA(qq,7),ZZ_FULL_RA(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RA(ss2,26) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(ss2,27) = num2cell(1);
    else
        ZZ_FULL_RA(end+1,1:20) = ZZ_test_RA(qq,:);
        ZZ_FULL_RA(end,26) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(end,27) = num2cell(1);
    end
end
disp('DONE2_3')

% OSM
[ZZ_test_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TLasTarget_RpR_OSMLeGO,TLasTarget_RpR_OSMLeGO,RNAseq_lvl42,RNAseqIDs,RPPA_lvl4,RPPAIDs,0);
[ZZ_test_RA,~,~,~] = ...
    TTcolorInts_RA_lvl4_v2(TLasTarget_RA_OSMLeGO,TLasTarget_RA_OSMLeGO,ATACseq_lvl42,ATACseqIDs,RNAseq_lvl42,RNAseqIDs,0);
testlen1 = size(ZZ_test_RpR,1);
testlen2 = size(ZZ_test_RA,1);
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
for qq = 1:testlen2
    stest1 = strcmp(ZZ_test_RA(qq,4),ZZ_FULL_RA(:,4));
    stest2 = strcmp(ZZ_test_RA(qq,7),ZZ_FULL_RA(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RA(ss2,28) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(ss2,29) = num2cell(1);
    else
        ZZ_FULL_RA(end+1,1:20) = ZZ_test_RA(qq,:);
        ZZ_FULL_RA(end,28) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(end,29) = num2cell(1);
    end
end
disp('DONE2_4')

% IFNG
[ZZ_test_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TLasTarget_RpR_IFNGLeGO,TLasTarget_RpR_IFNGLeGO,RNAseq_lvl42,RNAseqIDs,RPPA_lvl4,RPPAIDs,0);
[ZZ_test_RA,~,~,~] = ...
    TTcolorInts_RA_lvl4_v2(TLasTarget_RA_IFNGLeGO,TLasTarget_RA_IFNGLeGO,ATACseq_lvl42,ATACseqIDs,RNAseq_lvl42,RNAseqIDs,0);
testlen1 = size(ZZ_test_RpR,1);
testlen2 = size(ZZ_test_RA,1);
for qq = 1:testlen1
    stest1 = strcmp(ZZ_test_RpR(qq,4),ZZ_FULL_RpR(:,4));
    stest2 = strcmp(ZZ_test_RpR(qq,7),ZZ_FULL_RpR(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RpR(ss2,30) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(ss2,31) = num2cell(1);
    else
        ZZ_FULL_RpR(end+1,1:20) = ZZ_test_RpR(qq,:);
        ZZ_FULL_RpR(end,30) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(end,31) = num2cell(1);
    end
end
for qq = 1:testlen2
    stest1 = strcmp(ZZ_test_RA(qq,4),ZZ_FULL_RA(:,4));
    stest2 = strcmp(ZZ_test_RA(qq,7),ZZ_FULL_RA(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RA(ss2,30) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(ss2,31) = num2cell(1);
    else
        ZZ_FULL_RA(end+1,1:20) = ZZ_test_RA(qq,:);
        ZZ_FULL_RA(end,30) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(end,31) = num2cell(1);
    end
end
disp('DONE2_5')

% BMP2
[ZZ_test_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TLasTarget_RpR_BMP2LeGO,TLasTarget_RpR_BMP2LeGO,RNAseq_lvl42,RNAseqIDs,RPPA_lvl4,RPPAIDs,0);
[ZZ_test_RA,~,~,~] = ...
    TTcolorInts_RA_lvl4_v2(TLasTarget_RA_BMP2LeGO,TLasTarget_RA_BMP2LeGO,ATACseq_lvl42,ATACseqIDs,RNAseq_lvl42,RNAseqIDs,0);
testlen1 = size(ZZ_test_RpR,1);
testlen2 = size(ZZ_test_RA,1);
for qq = 1:testlen1
    stest1 = strcmp(ZZ_test_RpR(qq,4),ZZ_FULL_RpR(:,4));
    stest2 = strcmp(ZZ_test_RpR(qq,7),ZZ_FULL_RpR(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RpR(ss2,32) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(ss2,33) = num2cell(1);
    else
        ZZ_FULL_RpR(end+1,1:20) = ZZ_test_RpR(qq,:);
        ZZ_FULL_RpR(end,32) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(end,33) = num2cell(1);
    end
end
for qq = 1:testlen2
    stest1 = strcmp(ZZ_test_RA(qq,4),ZZ_FULL_RA(:,4));
    stest2 = strcmp(ZZ_test_RA(qq,7),ZZ_FULL_RA(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RA(ss2,32) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(ss2,33) = num2cell(1);
    else
        ZZ_FULL_RA(end+1,1:20) = ZZ_test_RA(qq,:);
        ZZ_FULL_RA(end,32) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(end,33) = num2cell(1);
    end
end
disp('DONE2_6')

% TGFB
[ZZ_test_RpR,~,~,~] = ...
    TTcolorInts_RpR_lvl4_v2(TLasTarget_RpR_TGFBLeGO,TLasTarget_RpR_TGFBLeGO,RNAseq_lvl42,RNAseqIDs,RPPA_lvl4,RPPAIDs,0);
[ZZ_test_RA,~,~,~] = ...
    TTcolorInts_RA_lvl4_v2(TLasTarget_RA_TGFBLeGO,TLasTarget_RA_TGFBLeGO,ATACseq_lvl42,ATACseqIDs,RNAseq_lvl42,RNAseqIDs,0);
testlen1 = size(ZZ_test_RpR,1);
testlen2 = size(ZZ_test_RA,1);
for qq = 1:testlen1
    stest1 = strcmp(ZZ_test_RpR(qq,4),ZZ_FULL_RpR(:,4));
    stest2 = strcmp(ZZ_test_RpR(qq,7),ZZ_FULL_RpR(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RpR(ss2,34) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(ss2,35) = num2cell(1);
    else
        ZZ_FULL_RpR(end+1,1:20) = ZZ_test_RpR(qq,:);
        ZZ_FULL_RpR(end,34) = ZZ_test_RpR(qq,8);
        ZZ_FULL_RpR(end,35) = num2cell(1);
    end
end
for qq = 1:testlen2
    stest1 = strcmp(ZZ_test_RA(qq,4),ZZ_FULL_RA(:,4));
    stest2 = strcmp(ZZ_test_RA(qq,7),ZZ_FULL_RA(:,7));
    ss2 = find(stest1+stest2==2);
    if ss2
        ZZ_FULL_RA(ss2,34) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(ss2,35) = num2cell(1);
    else
        ZZ_FULL_RA(end+1,1:20) = ZZ_test_RA(qq,:);
        ZZ_FULL_RA(end,34) = ZZ_test_RA(qq,8);
        ZZ_FULL_RA(end,35) = num2cell(1);
    end
end
disp('DONE2_7')

disp('DONE2_All')

%% Check the core of cores network - RpR
clc
[aa,bb] = find(cellfun('isempty', ZZ_FULL_RpR(:,21:end)));
ind = sub2ind(size(ZZ_FULL_RpR), aa, bb+20);
ZZ_FULL_RpR(ind) = num2cell(0);

coreIntsRpR = [(1:size(ZZ_FULL_RpR,1))', sum(cell2mat(ZZ_FULL_RpR(:,[21,23,25,27,29,31,33,35])),2)];
coreIntsRpR2 = (sortrows(coreIntsRpR,2));
coreIntsRpRList = ZZ_FULL_RpR(coreIntsRpR2(find(coreIntsRpR2(:,2)==8)),:);
disp('DONE3')

%% Check the core of cores network - RA
clc
[aa,bb] = find(cellfun('isempty', ZZ_FULL_RA(:,21:end)));
ind = sub2ind(size(ZZ_FULL_RA), aa, bb+20);
ZZ_FULL_RA(ind) = num2cell(0);

coreIntsRA = [(1:size(ZZ_FULL_RA,1))', sum(cell2mat(ZZ_FULL_RA(:,[21,23,25,27,29,31,33,35])),2)];
coreIntsRA2 = (sortrows(coreIntsRA,2));
coreIntsRAList = ZZ_FULL_RA(coreIntsRA2(find(coreIntsRA2(:,2)==8)),:);
disp('DONE4')

%% Rearrange final lists
ZZ_RpRList = {};
ZZ_RpRList(:,1:4) = ZZ_FULL_RpR(:,1:4);
ZZ_RpRList(:,5) = ZZ_FULL_RpR(:,16);
ZZ_RpRList(:,6:8) = ZZ_FULL_RpR(:,5:7);
ZZ_RpRList(:,9) = ZZ_FULL_RpR(:,17);
ZZ_RpRList(:,10) = ZZ_FULL_RpR(:,20);
ZZ_RpRList(:,11) = num2cell(coreIntsRpR(:,2));
ZZ_RpRList(:,12:19) = ZZ_FULL_RpR(:,21:2:35);
ZZ_RpRList(:,20:21) = [ZZ_FULL_RpR(:,8) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,8))))];
ZZ_RpRList(:,22:23) = [ZZ_FULL_RpR(:,22) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,22))))];
ZZ_RpRList(:,24:25) = [ZZ_FULL_RpR(:,24) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,24))))];
ZZ_RpRList(:,26:27) = [ZZ_FULL_RpR(:,26) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,26))))];
ZZ_RpRList(:,28:29) = [ZZ_FULL_RpR(:,28) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,28))))];
ZZ_RpRList(:,30:31) = [ZZ_FULL_RpR(:,30) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,30))))];
ZZ_RpRList(:,32:33) = [ZZ_FULL_RpR(:,32) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,32))))];
ZZ_RpRList(:,34:35) = [ZZ_FULL_RpR(:,34) num2cell(abs(cell2mat(ZZ_FULL_RpR(:,34))))];

ZZ_RAList = {};
ZZ_RAList(:,1:4) = ZZ_FULL_RA(:,1:4);
ZZ_RAList(:,5) = ZZ_FULL_RA(:,16);
ZZ_RAList(:,6:8) = ZZ_FULL_RA(:,5:7);
ZZ_RAList(:,9) = ZZ_FULL_RA(:,17);
ZZ_RAList(:,10) = ZZ_FULL_RA(:,20);
ZZ_RAList(:,11) = num2cell(coreIntsRA(:,2));
ZZ_RAList(:,12:19) = ZZ_FULL_RA(:,21:2:35);
ZZ_RAList(:,20:21) = [ZZ_FULL_RA(:,8) num2cell(abs(cell2mat(ZZ_FULL_RA(:,8))))];
ZZ_RAList(:,22:23) = [ZZ_FULL_RA(:,22) num2cell(abs(cell2mat(ZZ_FULL_RA(:,22))))];
ZZ_RAList(:,24:25) = [ZZ_FULL_RA(:,24) num2cell(abs(cell2mat(ZZ_FULL_RA(:,24))))];
ZZ_RAList(:,26:27) = [ZZ_FULL_RA(:,26) num2cell(abs(cell2mat(ZZ_FULL_RA(:,26))))];
ZZ_RAList(:,28:29) = [ZZ_FULL_RA(:,28) num2cell(abs(cell2mat(ZZ_FULL_RA(:,28))))];
ZZ_RAList(:,30:31) = [ZZ_FULL_RA(:,30) num2cell(abs(cell2mat(ZZ_FULL_RA(:,30))))];
ZZ_RAList(:,32:33) = [ZZ_FULL_RA(:,32) num2cell(abs(cell2mat(ZZ_FULL_RA(:,32))))];
ZZ_RAList(:,34:35) = [ZZ_FULL_RA(:,34) num2cell(abs(cell2mat(ZZ_FULL_RA(:,34))))];

disp('DONE5')

%% Write the lists into files
clc
submodelfile = 'coreIntsRpRList.xlsx' ; 
xlswrite(submodelfile,coreIntsRpRList)
submodelfile = 'coreIntsRAList.xlsx' ; 
xlswrite(submodelfile,coreIntsRAList)

submodelfile = 'ZZ_RpRList.xlsx' ; 
xlswrite(submodelfile,ZZ_RpRList)
submodelfile = 'ZZ_RAList.xlsx' ; 
xlswrite(submodelfile,ZZ_RAList)

disp('DONE6')

%% Save the lists
clc
save('TLasBestRpRA_Summ.mat','ZZ_RAList','ZZ_RpRList','-v7.3') 

disp('DONE7')
