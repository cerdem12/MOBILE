function [intranksBOTH_r,intranksONLY1_r,intranksONLY2_r,TTcolor] = ...
    TTcolorInts_RA_lvl4_v2(TLas1,TLas2,XX_Annots,XX_IDs,YY_Annots,YY_IDs,toplot)
if nargin<7
    toplot = 1;
end

[~,~,TTcolor] = plotTwoMatCombos(TLas1,TLas2,toplot);
if toplot
    c11 = colorbar;
    set(c11,'YTick',[0.5,1.5,2.5,3.5])
    set(c11,'YTickLabel',{'NONE';'FULL';'LeGO';'BOTH';},'FontSize',12,'FontWeight','bold')
    set(gcf, 'renderer', 'zbuffer');
    axis square
end

[zjONLY1only,ziONLY1only,~] = find(TTcolor==1);
[zjONLY2only,ziONLY2only,~] = find(TTcolor==2);
[zjBOTHonly,ziBOTHonly,~] = find(TTcolor==3);

numBOTHonly = size(ziBOTHonly,1);
intranksBOTH = {};
for qq = 1:numBOTHonly
intranksBOTH(qq,1) = num2cell(ziBOTHonly(qq));
intranksBOTH(qq,2) = num2cell(zjBOTHonly(qq));
intranksBOTH(qq,3) = cellstr(YY_Annots.ensembl_id(YY_IDs(ziBOTHonly(qq)),:));
intranksBOTH(qq,4) = cellstr(YY_Annots.hgnc_id(YY_IDs(ziBOTHonly(qq)),:));
intranksBOTH(qq,5) = cellstr(XX_Annots.ensembl_id(XX_IDs(zjBOTHonly(qq)),:));
intranksBOTH(qq,6) = cellstr(XX_Annots.hgnc_id(XX_IDs(zjBOTHonly(qq)),:));
intranksBOTH(qq,7) = cellstr(strjoin(([strtrim(XX_Annots.hgnc_id(XX_IDs(zjBOTHonly(qq)),:)) '_(' num2str((XX_Annots.dist2TSS{XX_IDs(zjBOTHonly(qq)),:})) ')'])));
intranksBOTH(qq,8) = num2cell(TLas1(ziBOTHonly(qq),zjBOTHonly(qq)));
intranksBOTH(qq,9) = num2cell(abs(TLas1(ziBOTHonly(qq),zjBOTHonly(qq))));
intranksBOTH(qq,10) = num2cell(TLas2(ziBOTHonly(qq),zjBOTHonly(qq)));
intranksBOTH(qq,11) = num2cell(abs(TLas2(ziBOTHonly(qq),zjBOTHonly(qq))));
intranksBOTH(qq,12) = num2cell(abs(TLas1(ziBOTHonly(qq),zjBOTHonly(qq)))+abs(TLas2(ziBOTHonly(qq),zjBOTHonly(qq))));
intranksBOTH(qq,13) = num2cell((TLas1(ziBOTHonly(qq),zjBOTHonly(qq)))./(TLas2(ziBOTHonly(qq),zjBOTHonly(qq))));
intranksBOTH(qq,14) = num2cell(abs(TLas1(ziBOTHonly(qq),zjBOTHonly(qq)))./abs(TLas2(ziBOTHonly(qq),zjBOTHonly(qq))));
intranksBOTH(qq,15) = num2cell((abs(TLas1(ziBOTHonly(qq),zjBOTHonly(qq)))+abs(TLas2(ziBOTHonly(qq),zjBOTHonly(qq))))./2);
intranksBOTH(qq,16) = num2cell(1); % means Y are RNAseq
intranksBOTH(qq,17) = num2cell(4); % means X are ATAC
intranksBOTH(qq,18) = num2cell(3); % means in-both
intranksBOTH(qq,19) = num2cell(max([abs(TLas1(ziBOTHonly(qq),zjBOTHonly(qq))),abs(TLas2(ziBOTHonly(qq),zjBOTHonly(qq)))]));
intranksBOTH(qq,20) = cellstr(XX_Annots.peak(XX_IDs(zjBOTHonly(qq)),:));
end
if size(intranksBOTH,1)
    intranksBOTH_r = flipud(sortrows(intranksBOTH,19));
    aa = find(strcmp('',intranksBOTH_r(:,4)));
    intranksBOTH_r(aa,4) = intranksBOTH_r(aa,3);
    aa2 = find(strcmp('NA',intranksBOTH_r(:,4)));
    intranksBOTH_r(aa2,4) = intranksBOTH_r(aa2,3);
    aa3 = find(strcmp('',intranksBOTH_r(:,6)));
    intranksBOTH_r(aa3,6) = intranksBOTH_r(aa3,5);
    aa4 = find(strcmp('NA',intranksBOTH_r(:,6)));
    intranksBOTH_r(aa4,6) = intranksBOTH_r(aa4,5);
else
    intranksBOTH_r = {};
end

numONLY1only = size(ziONLY1only,1);
intranksONLY1 = {};
for qq = 1:numONLY1only
intranksONLY1(qq,1) = num2cell(ziONLY1only(qq));
intranksONLY1(qq,2) = num2cell(zjONLY1only(qq));
intranksONLY1(qq,3) = cellstr(YY_Annots.ensembl_id(YY_IDs(ziONLY1only(qq)),:));
intranksONLY1(qq,4) = cellstr(YY_Annots.hgnc_id(YY_IDs(ziONLY1only(qq)),:));
intranksONLY1(qq,5) = cellstr(XX_Annots.ensembl_id(XX_IDs(zjONLY1only(qq)),:));
intranksONLY1(qq,6) = cellstr(XX_Annots.hgnc_id(XX_IDs(zjONLY1only(qq)),:));
intranksONLY1(qq,7) = cellstr(strjoin(([strtrim(XX_Annots.hgnc_id(XX_IDs(zjONLY1only(qq)),:)) '_(' num2str((XX_Annots.dist2TSS{XX_IDs(zjONLY1only(qq)),:})) ')'])));
intranksONLY1(qq,8) = num2cell(TLas1(ziONLY1only(qq),zjONLY1only(qq)));
intranksONLY1(qq,9) = num2cell(abs(TLas1(ziONLY1only(qq),zjONLY1only(qq))));
intranksONLY1(qq,10) = num2cell(TLas2(ziONLY1only(qq),zjONLY1only(qq)));
intranksONLY1(qq,11) = num2cell(abs(TLas2(ziONLY1only(qq),zjONLY1only(qq))));
intranksONLY1(qq,12) = num2cell(abs(TLas1(ziONLY1only(qq),zjONLY1only(qq)))+abs(TLas2(ziONLY1only(qq),zjONLY1only(qq))));
intranksONLY1(qq,13) = num2cell((TLas1(ziONLY1only(qq),zjONLY1only(qq)))./(TLas2(ziONLY1only(qq),zjONLY1only(qq))));
intranksONLY1(qq,14) = num2cell(abs(TLas1(ziONLY1only(qq),zjONLY1only(qq)))./abs(TLas2(ziONLY1only(qq),zjONLY1only(qq))));
intranksONLY1(qq,15) = num2cell((abs(TLas1(ziONLY1only(qq),zjONLY1only(qq)))+abs(TLas2(ziONLY1only(qq),zjONLY1only(qq))))./2);
intranksONLY1(qq,16) = num2cell(1);
intranksONLY1(qq,17) = num2cell(4);
intranksONLY1(qq,18) = num2cell(1); % means in-first(FULL)-only
intranksONLY1(qq,19) = num2cell(max([abs(TLas1(ziONLY1only(qq),zjONLY1only(qq))),abs(TLas2(ziONLY1only(qq),zjONLY1only(qq)))]));
intranksONLY1(qq,20) = cellstr(XX_Annots.peak(XX_IDs(zjONLY1only(qq)),:));
end
if size(intranksONLY1,1)
    intranksONLY1_r = flipud(sortrows(intranksONLY1,19));
    bb = find(strcmp('',intranksONLY1_r(:,4)));
    intranksONLY1_r(bb,4) = intranksONLY1_r(bb,3);
    bb2 = find(strcmp('NA',intranksONLY1_r(:,4)));
    intranksONLY1_r(bb2,4) = intranksONLY1_r(bb2,3);
    bb3 = find(strcmp('',intranksONLY1_r(:,6)));
    intranksONLY1_r(bb3,6) = intranksONLY1_r(bb3,5);
    bb4 = find(strcmp('NA',intranksONLY1_r(:,6)));
    intranksONLY1_r(bb4,6) = intranksONLY1_r(bb4,5);
else
    intranksONLY1_r = {};
end

numONLY2only = size(ziONLY2only,1);
intranksONLY2 = {};
for qq = 1:numONLY2only
intranksONLY2(qq,1) = num2cell(ziONLY2only(qq));
intranksONLY2(qq,2) = num2cell(zjONLY2only(qq));
intranksONLY2(qq,3) = cellstr(YY_Annots.ensembl_id(YY_IDs(ziONLY2only(qq)),:));
intranksONLY2(qq,4) = cellstr(YY_Annots.hgnc_id(YY_IDs(ziONLY2only(qq)),:));
intranksONLY2(qq,5) = cellstr(XX_Annots.ensembl_id(XX_IDs(zjONLY2only(qq)),:));
intranksONLY2(qq,6) = cellstr(XX_Annots.hgnc_id(XX_IDs(zjONLY2only(qq)),:));
intranksONLY2(qq,7) = cellstr(strjoin(([strtrim(XX_Annots.hgnc_id(XX_IDs(zjONLY2only(qq)),:)) '_(' num2str((XX_Annots.dist2TSS{XX_IDs(zjONLY2only(qq)),:})) ')'])));
intranksONLY2(qq,8) = num2cell(TLas1(ziONLY2only(qq),zjONLY2only(qq)));
intranksONLY2(qq,9) = num2cell(abs(TLas1(ziONLY2only(qq),zjONLY2only(qq))));
intranksONLY2(qq,10) = num2cell(TLas2(ziONLY2only(qq),zjONLY2only(qq)));
intranksONLY2(qq,11) = num2cell(abs(TLas2(ziONLY2only(qq),zjONLY2only(qq))));
intranksONLY2(qq,12) = num2cell(abs(TLas1(ziONLY2only(qq),zjONLY2only(qq)))+abs(TLas2(ziONLY2only(qq),zjONLY2only(qq))));
intranksONLY2(qq,13) = num2cell((TLas1(ziONLY2only(qq),zjONLY2only(qq)))./(TLas2(ziONLY2only(qq),zjONLY2only(qq))));
intranksONLY2(qq,14) = num2cell(abs(TLas1(ziONLY2only(qq),zjONLY2only(qq)))./abs(TLas2(ziONLY2only(qq),zjONLY2only(qq))));
intranksONLY2(qq,15) = num2cell((abs(TLas1(ziONLY2only(qq),zjONLY2only(qq)))+abs(TLas2(ziONLY2only(qq),zjONLY2only(qq))))./2);
intranksONLY2(qq,16) = num2cell(1);
intranksONLY2(qq,17) = num2cell(4);
intranksONLY2(qq,18) = num2cell(2); % means in-second(LeGO)-only
intranksONLY2(qq,19) = num2cell(max([abs(TLas1(ziONLY2only(qq),zjONLY2only(qq))),abs(TLas2(ziONLY2only(qq),zjONLY2only(qq)))]));
intranksONLY2(qq,20) = cellstr(XX_Annots.peak(XX_IDs(zjONLY2only(qq)),:));
end
if size(intranksONLY2,1)
    intranksONLY2_r = flipud(sortrows(intranksONLY2,19));
    cc = find(strcmp('',intranksONLY2_r(:,4)));
    intranksONLY2_r(cc,4) = intranksONLY2_r(cc,3);
    cc2 = find(strcmp('NA',intranksONLY2_r(:,4)));
    intranksONLY2_r(cc2,4) = intranksONLY2_r(cc2,3);
    cc3 = find(strcmp('',intranksONLY2_r(:,6)));
    intranksONLY2_r(cc3,6) = intranksONLY2_r(cc3,5);
    cc4 = find(strcmp('NA',intranksONLY2_r(:,6)));
    intranksONLY2_r(cc4,6) = intranksONLY2_r(cc4,5);
else
    intranksONLY2_r = {};
end
return
