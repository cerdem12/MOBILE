function [sumTG1, sumTG2, TTcolor] = plotTwoMatCombos(TG1,TG2,toplot)

if nargin<3
    toplot = 1;
end

TG1 = transpose(TG1);
TG1(abs(TG1)>0) = 1;
TG2 = transpose(TG2);
TG2(abs(TG2)>0) = 1;
norows = size(TG1,1);
nocols = size(TG1,2);
%%%%% 0=not in both // 1=in only igf1 // 2=in only ins // 3=in both
TTcolor = zeros(norows,nocols);
for ppp = 1:norows
    for rrr = 1:nocols
        if (TG1(ppp,rrr)==0) && (TG2(ppp,rrr)==0)
            TTcolor(ppp,rrr) = 0;
        elseif (TG1(ppp,rrr)==1) && (TG2(ppp,rrr)==0)
            TTcolor(ppp,rrr) = 1;
        elseif (TG1(ppp,rrr)==0) && (TG2(ppp,rrr)==1)
            TTcolor(ppp,rrr) = 2;   
        elseif (TG1(ppp,rrr)==1) && (TG2(ppp,rrr)==1)
            TTcolor(ppp,rrr) = 3;
        end
    end
end
sumTG1 = sum(sum(TG1));
sumTG2 = sum(sum(TG2));
% keyboard
if toplot
    figure;
    pcolor(([(TTcolor) zeros(size(TTcolor,1),1);zeros(1,size(TTcolor,2)+1)]))
    map = [1, 1, 1; 0, 0, 1; 1 1 0; 1 0 0];
    colormap(map);
    caxis([0, 4])
    ylabel('Mat\_X analytes')
    xlabel('Mat\_Y analytes')
    set(gcf, 'renderer', 'zbuffer');
    colorbar
end
return