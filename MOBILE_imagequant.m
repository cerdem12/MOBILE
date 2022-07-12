%% Load data
load('MOBILE_IFexptsData.mat') % MOBILE valdiation experiment data

LINCSexpData = readtable('MDD_IF_Level3.csv'); % LINCS IF data
reodridx = [14,15,16,17,18,19,20,73,74,75,76,77,78,79,80,81,82,83,84,85,21, ...
            22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42, ...
            43,44,45,46,60,61,62,63,64,65,66,67,68,69,70,71,72,47,48,49,50, ...
            51,52,53,54,55,56,57,58,59,1,2,3,4,5,6,7,8,9,10,11,12,13,86,87, ...
            88,89,90,91,92,93,94,95,96,97,98]; % reorder columns for ligand/time point pairs
LINCSexpData2 = table2cell(LINCSexpData);
LINCSexpData2mes = LINCSexpData2(:,1);
LINCSexpData2 = LINCSexpData2(:,2:end); 
LINCSexpData2ord = LINCSexpData2(:,reodridx); % Finalized LINCS data

%% MOBILE data pre-processing
run1nonctExc = [];
run2nonctExc = [7,10,13,12,15]; % Images excluded because of low quality (out of focus or partial)
run1DataNonCT = cell2mat(LassoExpData1(:,setdiff(1:size(LassoExpData1,2),run1nonctExc)));
run2DataNonCT = cell2mat(LassoExpData3(:,setdiff(1:size(LassoExpData3,2),run2nonctExc)));

%%% Ctrl condition: EGF @48hr
% Well cell count metric
WC_ctrlmean1 = mean((run1DataNonCT(15,1:6))); % EGF cond as CTRL
WC_ctrlmean2 = mean((run2DataNonCT(15,1:6)));
WellCounts1 = (run1DataNonCT(15,:))./WC_ctrlmean1;
WellCounts2 = (run2DataNonCT(15,:))./WC_ctrlmean2;

% Mean cells per cluster metric
MCpC_ctrlmean1 = mean((run1DataNonCT(14,1:6))); % EGF cond as CTRL
MCpC_ctrlmean2 = mean((run2DataNonCT(14,1:6)));
CellperClster1 = (run1DataNonCT(14,:))./MCpC_ctrlmean1;
CellperClster2 = (run2DataNonCT(14,:))./MCpC_ctrlmean2;

% Distance to second nearest neighbor metric
SND_ctrlmean1 = mean((run1DataNonCT(8,1:6))); % EGF cond as CTRL
SND_ctrlmean2 = mean((run2DataNonCT(8,1:6)));
SecNeighDist1 = (run1DataNonCT(8,:))./SND_ctrlmean1;
SecNeighDist2 = (run2DataNonCT(8,:))./SND_ctrlmean2;

%%% Ligand condition indexing - MOBILE
Eidx_nCT1 = 1:6; % EGF
EBidx_nCT1 = 7:11; % EGF+BMP2
ETidx_nCT1 = 12:14; % EGF+TGFB1
Eidx_nCT2 = 1:6; % EGF
EBidx_nCT2 = 7:9; % EGF+BMP2
ETidx_nCT2 = 10; % EGF+TGFB1

% Concatanate experimental replicates based on ligands
WellCountsCTnon_E = [WellCounts1(Eidx_nCT1),WellCounts2(Eidx_nCT2)]; 
WellCountsCTnon_EB = [WellCounts1(EBidx_nCT1),WellCounts2(EBidx_nCT2)];
WellCountsCTnon_ET = [WellCounts1(ETidx_nCT1),WellCounts2(ETidx_nCT2)];
CellperClsterCTnon_E = [CellperClster1(Eidx_nCT1),CellperClster2(Eidx_nCT2)];
CellperClsterCTnon_EB = [CellperClster1(EBidx_nCT1),CellperClster2(EBidx_nCT2)];
CellperClsterCTnon_ET = [CellperClster1(ETidx_nCT1),CellperClster2(ETidx_nCT2)];
SecNeighDistCTnon_E = [SecNeighDist1(Eidx_nCT1),SecNeighDist2(Eidx_nCT2)];
SecNeighDistCTnon_EB = [SecNeighDist1(EBidx_nCT1),SecNeighDist2(EBidx_nCT2)];
SecNeighDistCTnon_ET = [SecNeighDist1(ETidx_nCT1),SecNeighDist2(ETidx_nCT2)];

%% LINCS data pre-processing - normalize to EGF @48hr
ctrl_C1mean_WC = mean(str2double(LINCSexpData2ord(6,28:30)));
ctrl_C2mean_WC = mean(str2double(LINCSexpData2ord(6,31:33)));

ctrl_C1mean_MCpC = mean(str2double(LINCSexpData2ord(74,28:30)));
ctrl_C2mean_MCpC = mean(str2double(LINCSexpData2ord(74,31:33)));

ctrl_C1mean_SND = mean(str2double(LINCSexpData2ord(70,28:30)));
ctrl_C2mean_SND = mean(str2double(LINCSexpData2ord(70,31:33)));

%%% LINCS cells were collected in different times called Collections 1 & 2
%%% Define different collection data columns and use them for normalization
C1sampleIdx = [1,2,3,8,9,10,15,16,17,21,22,23,28,29,30,34,35,36,41,42,43, ...
               47,48,49,54,55,56,60,61,62,67,68,69,73,74,75,80,81,82,86,87,88,93,94,95];
C2sampleIdx = setdiff((1:98),C1sampleIdx); 

% Well cell count metric
WellCountsC1 = str2double(LINCSexpData2ord(6,:))./ctrl_C1mean_WC;
WellCountsC2 = str2double(LINCSexpData2ord(6,:))./ctrl_C2mean_WC;
WellCountsLINCS = zeros(1,98);
WellCountsLINCS(C1sampleIdx) = WellCountsC1(C1sampleIdx);
WellCountsLINCS(C2sampleIdx) = WellCountsC2(C2sampleIdx);

% Mean cells per cluster metric
CellperClsterC1 = str2double(LINCSexpData2ord(74,:))./ctrl_C1mean_MCpC;
CellperClsterC2 = str2double(LINCSexpData2ord(74,:))./ctrl_C2mean_MCpC;
CellperClsterLINCS = zeros(1,98);
CellperClsterLINCS(C1sampleIdx) = CellperClsterC1(C1sampleIdx);
CellperClsterLINCS(C2sampleIdx) = CellperClsterC2(C2sampleIdx);

% Distance to second nearest neighbor metric
SecNeighDistC1 = str2double(LINCSexpData2ord(70,:))./ctrl_C1mean_SND;
SecNeighDistC2 = str2double(LINCSexpData2ord(70,:))./ctrl_C2mean_SND;
SecNeighDistLINCS = zeros(1,98);
SecNeighDistLINCS(C1sampleIdx) = SecNeighDistC1(C1sampleIdx);
SecNeighDistLINCS(C2sampleIdx) = SecNeighDistC2(C2sampleIdx);

% Ligand condition indexes - LINCS
ctrl_C1 = [1,2,3];	
ctrl_C2  = [4,5,6,7];
PBS_C1 = [15,16,17];	
PBS_C2 = [18,19,20];				
EGF_C1 = [28,29,30];
EGF_C2 = [31,32,33];
IFNG_C1 = [67,68,69];
IFNG_C2 = [70,71,72];				
BMP2_C1 = [80,81,82];
BMP2_C2 = [83,84,85];
TGFB1_C1 = [93,94,95];
TGFB1_C2 = [96,97,98];

%% Data plotting
% Grouping parameter for box plots 
g1 = repmat({'BMP2-LINCS'},length((WellCountsLINCS([BMP2_C1,BMP2_C2]))),1);
g2 = repmat({'TGFB1-LINCS'},length((WellCountsLINCS([TGFB1_C1,TGFB1_C2]))),1);
g3 = repmat({'BMP2'},length(WellCountsCTnon_EB),1);
g4 = repmat({'TGFB1'},length(WellCountsCTnon_ET),1);
gg = [g1; g2; g3; g4];

% Plot well cell count data in the first panel
aa = [(WellCountsLINCS([BMP2_C1,BMP2_C2]))';
      (WellCountsLINCS([TGFB1_C1,TGFB1_C2]))';
      WellCountsCTnon_EB';
      WellCountsCTnon_ET'];
% Plot cells per cluster data in the second panel
bb = [(CellperClsterLINCS([BMP2_C1,BMP2_C2]))';
(CellperClsterLINCS([TGFB1_C1,TGFB1_C2]))';
CellperClsterCTnon_EB';
CellperClsterCTnon_ET'];
% Plot neighbor distance data in the third panel
cc = [(SecNeighDistLINCS([BMP2_C1,BMP2_C2]))';
(SecNeighDistLINCS([TGFB1_C1,TGFB1_C2]))';
SecNeighDistCTnon_EB';
SecNeighDistCTnon_ET'];

% custom color scheme
colors = [0.4375,0.1875,0.6250;
          1,0.7529,0;
          0.4375,0.1875,0.8250;
          1,0.7529,0.2];
ids = [4;3;2;1];
FntSze = 16;

% Definitions for patching
x1 = [0,0,0,0];
x2 = [2.5,2.5,2.5,2.5];
y1 = [0,1,2,5];
y2 = y1;

figure;
subplot(1,3,1); hold on;
% Area patch color to separate LINCS (coated plates) from MOBILE (non-coated exp. plates) data 
patch([x1 fliplr(x2)], [y1 fliplr(y2)], [0.1,0.1,0.1],'FaceAlpha',.15,'LineStyle','none')
boxplot(aa,gg,'ColorGroup',gg,'colors',colors,'symbol', '')
ylim([0 2])
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    tt = ids(j);
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(tt,:),'FaceAlpha',.5); % Individual bar patching
end
[C, ~, ic]= unique([gg],'stable');
h2 = findobj(gca,'tag','Median');
set(h2,'linewidth',3);
set(h2,'Color',[0 0 0])
% Individual data points
scatter(ic,aa,50,colors(ic,:),'filled','MarkerFaceAlpha',0.7','jitter','on','jitterAmount',0.15);
set(gca,'FontWeight','bold','FontSize',FntSze,'FontName','Arial')
set(gca,'XTickLabels',[])
set(gca,'YTickLabels',[])
plot(0:0.1:5,ones(size(0:0.1:5)),'k-.','LineWidth',1.3) % y=1 line (EGF condition)

% Plot second metric
subplot(1,3,2); hold on;
patch([x1 fliplr(x2)], [y1 fliplr(y2)], [0.1,0.1,0.1],'FaceAlpha',.15,'LineStyle','none')
boxplot(bb,gg,'ColorGroup',gg,'colors',colors,'symbol','')
ylim([0 4])
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    tt = ids(j);
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(tt,:),'FaceAlpha',.5);
end
h2 = findobj(gca,'tag','Median');
set(h2,'linewidth',3);
set(h2,'Color',[0 0 0])
scatter(ic,bb,50,colors(ic,:),'filled','MarkerFaceAlpha',0.7','jitter','on','jitterAmount',0.15);
set(gca,'FontWeight','bold','FontSize',FntSze,'FontName','Arial')
set(gca,'XTickLabels',[])
set(gca,'YTickLabels',[])
plot(0:0.1:5,ones(size(0:0.1:5)),'k-.','LineWidth',1.3) % y=1 line (EGF condition)

% Plot third metric
subplot(1,3,3); hold on;
patch([x1 fliplr(x2)], [y1 fliplr(y2)], [0.1,0.1,0.1],'FaceAlpha',.15,'LineStyle','none')
boxplot(cc,gg,'ColorGroup',gg,'colors',colors,'symbol','')
ylim([0 2])
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    tt = ids(j);
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(tt,:),'FaceAlpha',.5);
end
h2 = findobj(gca,'tag','Median');
set(h2,'linewidth',3);
set(h2,'Color',[0 0 0])
scatter(ic,cc,50,colors(ic,:),'filled','MarkerFaceAlpha',0.7','jitter','on','jitterAmount',0.15);
set(gca,'FontWeight','bold','FontSize',FntSze,'FontName','Arial')
set(gca,'XTickLabels',[])
set(gca,'YTickLabels',[])
plot(0:0.1:5,ones(size(0:0.1:5)),'k-.','LineWidth',1.3) % y=1 line (EGF condition)





