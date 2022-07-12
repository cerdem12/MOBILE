function [intranks_r,intranks_ii,numedges,TLasTarget] = runAssocRankerRA(cutoffval,TLasMats,XX_Annots,XX_IDs,YY_Annots,YY_IDs)
%%%% Function that finds Robust Lasso Coefficient Matrix 

%%%% Inputs
% cutoffval: Cut-off (between 0-1) to define robustness. Defines the number of occurences 
%            an association must appear to be coined as robust.
% TLasMats: The list of Lasso coefficient matrices (usually 10000 many)
% XX_Data: The data and annotations of the RHS matrix
% XX_IDs: The row indices for RHS input data
% YY_Data: The data and annotations of the LHS matrix
% YY_IDs: The row indices for LHS input data

%%%% Outputs
% intranks_r: List of associations in the Robust Lasso Coefficients Matrix
% intranks_ii: List of "self" associations (i.e. between same gene products)
% numedges: Number of non-zero associations inferred in each matrix
% TLasTarget: The ensemble median matrix with robust association locations defined -> the Target matrix

%%%% see the function "runAssocRankerRpR" for details below

numnetworks = size(TLasMats,1);
numedges = zeros(size(TLasMats,1),1);
TLasAllsum_1s = zeros(size(TLasMats{1,1}));
TLasAllsum = zeros(size(TLasMats{1,1}));

for qq = 1:numnetworks
    numedges(qq,1) = nnz(TLasMats{qq,1});
    TLastemp = (TLasMats{qq,1});
    TLasAllsum = TLasAllsum + TLastemp;    
    TLastemp(abs(TLastemp)>0) = 1;
    TLasAllsum_1s = TLasAllsum_1s + TLastemp;
end
cutofflimit = cutoffval;
cutofflimit = floor(numnetworks*cutofflimit);
TLas = (TLasAllsum_1s)>=(cutofflimit);

TLasAllmean = TLasAllsum./TLasAllsum_1s;
TLasAllmean(isnan(TLasAllmean)) = 0;
TLasTarget = TLas.*TLasAllmean;
TLasTarget(isnan(TLasTarget)) = 0;

nn = find(cellfun(@isnumeric, XX_Annots.hgnc_id));
XX_Annots.hgnc_id(nn) = XX_Annots.ensembl_id(nn);

[ii,jj,ss] = find(TLasTarget);
intranks = {};
intranks(:,1) = num2cell(jj);
intranks(:,2) = num2cell(ii);
intranks(:,3) = cellstr(XX_Annots.ensembl_id(XX_IDs(jj),:));
intranks(:,4) = cellstr(XX_Annots.hgnc_id(XX_IDs(jj),:));
intranks(:,5) = cellstr(num2str(cell2mat(XX_Annots.entrez_id(XX_IDs(jj)))));
intranks(:,6) = cellstr(YY_Annots.ensembl_id(YY_IDs(ii),:));
intranks(:,7) = cellstr(YY_Annots.hgnc_id(YY_IDs(ii),:));
intranks(:,8) = cellstr(' ');%cellstr(YY_Annots.entrez_gene(YY_IDs(ii),:));
intranks(:,9) = num2cell(ss);
intranks(:,10) = num2cell(abs(ss));
intranks(:,11) = cellstr(XX_Annots.peak(XX_IDs(jj),:));
intranks(:,12) = cellstr(cellstr(num2str(cell2mat(XX_Annots.dist2TSS(XX_IDs(jj),:)))));
intranks(:,13) = cellstr(XX_Annots.annot(XX_IDs(jj),:));

intranks_r = flipud(sortrows(intranks,10));
aa = find(strcmp('',intranks_r(:,4)));
intranks_r(aa,4) = intranks_r(aa,3);
aa = find(strcmp('NA',intranks_r(:,4)));
intranks_r(aa,4) = intranks_r(aa,3);
iinornot = find(strcmp(cellstr(intranks_r(:,4)),cellstr(intranks_r(:,7))));
intranks_ii = intranks_r(iinornot,:);
return
