%%%% This script helps to subset large associations list using specific
%%%% gene lists or based on single ligand dependence

%% Here load the data files
clc
rng(6);
warning('off','all')

addpath('glmnet','data','funcs')

load('RNAseq_RPPA_lvl4_data')
load('ATACseq_RNAseq_lvl4_data')
load('TLasBestRA_FULL.mat')
load('TLasBestRpR_FULL.mat')

load('The data and summary files')

%% color scheme for plots
colororder = [0.0000    0.4470    0.7410
              0.8500    0.3250    0.0980
              0.9290    0.6940    0.1250
              0.4940    0.1840    0.5560
              0.4660    0.6740    0.1880
              0.3010    0.7450    0.9330
              0.6350    0.0780    0.1840];
colororder = [colororder;colororder];

%% Pre-process input data for MOBILE input format
ATACseqraw = ATACseq_lvl42.data;
ATACvar = flipud(sortrows([std(ATACseqraw,1,2),(1:length(ATACseqraw))']));
ATACindices2keep = floor(0.1*size(ATACseqraw,1));
ATACseq0 = ATACseqraw(sortrows(ATACvar(1:ATACindices2keep,2)),:);
ATACseqIDs = sortrows(ATACvar(1:ATACindices2keep,2));
ATACgenes = ATACseq_lvl42.hgnc_id(ATACseqIDs,1);

RNAseqraw = RNAseq_lvl42.data;
RNAvar = flipud(sortrows([std(RNAseqraw,1,2),(1:length(RNAseqraw))']));
RNAindices2keep = floor(0.1*size(RNAseqraw,1));
RNAseq0 = RNAseqraw(sortrows(RNAvar(1:RNAindices2keep,2)),:);
RNAseqIDs = sortrows(RNAvar(1:RNAindices2keep,2));
RNAgenes = RNAseq_lvl42.hgnc_id(RNAseqIDs,1);

RPPAraw = RPPA_lvl4.data;
RPPAvar = flipud(sortrows([std(RPPAraw,1,2),(1:length(RPPAraw))']));
RPPAindices2keep = floor(0.2*size(RPPAraw,1));
RPPA0 = RPPAraw(sortrows(RPPAvar(1:RPPAindices2keep,2)),:);
RPPAIDs = sortrows(RPPAvar(1:RPPAindices2keep,2));
RPPAgenes = RPPA_lvl4.hgnc_id(RPPAIDs,1);

ATACseq3 = prepnormmats2(ATACseq0,6,1);
RNAseq3 = prepnormmats2(RNAseq0,6,1);
RPPA3 = prepnormmats2(RPPA0,6,1);

%% Check if there was reading error (genes as dates :D )
ENGids = {'ENSG00000108387';'ENSG00000184702';'ENSG00000125354'};
HGNCids = {'SEPT4';'SEPTIN5';'SEPT6'};
RNAgenesENS = RNAseq_lvl42.ensembl_id(RNAseqIDs,1);
for qq = 1:size(ENGids)
    tempg = ENGids(qq,1);
    aa = strcmp(RNAgenesENS,tempg);
    if find(aa)
        RNAgenes(find(aa)) = HGNCids(qq);
    end
end

%% Create custom input lists
% % Find the IFNG-dependent associations
Nn1 = ZZ_RpRList;  
Nn1(:,38) = num2cell(cell2mat(Nn1(:,12))+cell2mat(Nn1(:,17)));
Nn1(:,39) = num2cell(cell2mat(Nn1(:,21))+cell2mat(Nn1(:,31)));
Nn2 = Nn1(find(cell2mat(Nn1(:,38))==1),:); 
Nn3 = ZZ_RAList;  
Nn3(:,38) = num2cell(cell2mat(Nn3(:,12))+cell2mat(Nn3(:,17)));
Nn3(:,39) = num2cell(cell2mat(Nn3(:,21))+cell2mat(Nn3(:,31)));
Nn4 = Nn3(find(cell2mat(Nn3(:,38))==1),:); 
IFNG_LOGO = [Nn2;Nn4];

%% Find weighted coefficient lists for RNA-RPPA and RNA-ATAC matrices 
%%%% FULL=12 / PBS=13 / EGF=14 / HGF=15 / OSM=16 / IFNG=17 / BMP2=18 / TGFB1=19
%%%% FULL=21 / PBS=23 / EGF=25 / HGF=27 / OSM=29 / IFNG=31 / BMP2=33 / TGFB1=35
clc
colno4logo = 38; % full or logo case assoc. existence column no
colno4mags = 39; % assoc. magnitude column no
magcutoff = 0.01; % value cutoff for individual coeffs to pass 

uqList_F1 = {};
uqList_F2 = {};

Network2check = Nn2; % ZZ_RpRList; %      
num2normlz = 9321*59*3062; % number of possible coeffs in RNA-RPPA matrix
[uqList_F1,~] = findNetworkGeneWeights(Network2check,colno4logo,colno4mags,num2normlz,magcutoff);
Network2check = Nn4; % ZZ_RAList; % NnBTc3; % Nn4; %     
num2normlz = 9321*59*3062; % number of possible coeffs in RNA-ATAC matrix
[uqList_F2,~] = findNetworkGeneWeights(Network2check,colno4logo,colno4mags,num2normlz,magcutoff);

%% Combine the two (RNA-RPPA + RNA-ATAC) lists and sum the weighted coeffs
clc
uqList_F3 = uqList_F2;
for qq = 1:length(uqList_F1)
    tempG = uqList_F1(qq,1);
    aa = find(strcmp(uqList_F3(:,1),tempG));
    if aa
        uqList_F3(aa,2) = num2cell(cell2mat(uqList_F3(aa,2))+cell2mat(uqList_F1(qq,2)));
    else
        uqList_F3(end+1,:) = uqList_F1(qq,:);
    end
end

%% Write to tsv file
clc
flname = 'Test123.txt';
wrtlist = uqList_F3;

fid = fopen(flname,'w');
for k=1:length(wrtlist)
   fprintf(fid,'%s\t%s\n',wrtlist{k,1},wrtlist{k,2});
end
fclose(fid);

%% Remove genes without HGNC
clc
TempL0 = readtable(flname);
TempL1 = TempL0{:,1};
TempL2 = TempL0{:,2};
aaf = find(~cellfun(@isempty,strfind(TempL1,'ENSG0')));
TempL1(aaf,:) = [];
TempL2(aaf,:) = [];

fid = fopen([flname(1:end-4) '_filt.txt'],'w');
for k=1:length(TempL1)
   fprintf(fid,'%s\t%s\n',TempL1{k,1},TempL2(k,1));
end
fclose(fid);

%% Find interactions of custom genes in the Lasso/LOGO lists
Network2check = IFNG_LOGO; 
MagCutOff = 0.01;
colofMags = 39;
size(Network2check)
listofGenes = {'IFNG';'IFNGR1';'IFNGR2';'STAT1';'STAT3';'JAK1';'JAK2'; ...
               'PD-L1';'CD274';'IRF1';'IRF9';'PD-1';'CD279';'PDCD1'};
% listofGenes = {'STAT1';'CD274';'IRF1';'CD279'};
listofGenes(strcmp('',listofGenes)) = [];
listofColms2Check = [4;7;8];

Network2check = flipud(sortrows(Network2check,colofMags)); % sort based on max magnitude of in Full and LOGO
Network2check = Network2check(find(cell2mat(Network2check(:,colofMags))>=MagCutOff),:);
size(Network2check)

listofRows = [];
for qq = 1:size(listofGenes)
    for ww = 1:size(listofColms2Check)
        listofRows = [listofRows;find(strcmp(listofGenes{qq,1},Network2check(:,listofColms2Check(ww))))];
    end 
end
listofRows2 = unique(listofRows);
submodel = Network2check(listofRows2,:);

disp('DONE finding')

%% Eliminate non-curated gene associations
aa1 = find(~cellfun(@isempty,strfind(submodel(:,4),'ENSG0')));
aa2 = find(~cellfun(@isempty,strfind(submodel(:,7),'ENSG0')));
aa3 = unique([aa1;aa2]);
submodel(aa3,:) = [];

%%
clc
flname = 'IFNG_PDL1_mgt01.txt';
wrtlist = cell2table(submodel2);
writetable(wrtlist,flname,'Delimiter','\t')  

disp('DONE writing')







