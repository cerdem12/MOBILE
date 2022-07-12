function [intranks_r,intranks_ii,numedges,TLasTarget] = runAssocRankerRpR(cutoffval,TLasMats,XX_Data,XX_IDs,YY_Data,YY_IDs)
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

numnetworks = size(TLasMats,1); % Number of Lasso matrices
numedges = zeros(size(TLasMats,1),1); % Start empty array to store numbers of associations inferred
TLasAllsum_1s = zeros(size(TLasMats{1,1})); % Matrix to only store association locations
TLasAllsum = zeros(size(TLasMats{1,1})); % Matrix to store association magnitudes

for qq = 1:numnetworks % processing all coefficient matrices
    numedges(qq,1) = nnz(TLasMats{qq,1});
    TLastemp = (TLasMats{qq,1});
    TLasAllsum = TLasAllsum + TLastemp;    
    TLastemp(abs(TLastemp)>0) = 1;
    TLasAllsum_1s = TLasAllsum_1s + TLastemp;
end
cutofflimit = cutoffval;
cutofflimit = floor(numnetworks*cutofflimit);
TLas = (TLasAllsum_1s)>=(cutofflimit); % The target matrix (will be mostly zeros, ones will denote 
%                                        the associations that appeared more than cutoff times)

TLasAllmean = TLasAllsum./TLasAllsum_1s; % Calculate mean value for associations
TLasAllmean(isnan(TLasAllmean)) = 0;
TLasTarget = TLas.*TLasAllmean; % The Robust matrix with the most number of target associations
TLasTarget(isnan(TLasTarget)) = 0;

% Create a list that shows association values, genes, and other annotations
[ii,jj,ss] = find(TLasTarget); 
intranks = {};
intranks(:,1) = num2cell(jj);
intranks(:,2) = num2cell(ii);
intranks(:,3) = cellstr(XX_Data.ensembl_id(XX_IDs(jj),:));
intranks(:,4) = cellstr(XX_Data.hgnc_id(XX_IDs(jj),:));
intranks(:,5) = cellstr(' ');
intranks(:,6) = cellstr(' ');%cellstr(YY_Annots.ensembl_gene_id(YY_IDs(ii),:));
intranks(:,7) = cellstr(YY_Data.hgnc_id(YY_IDs(ii),:));
intranks(:,8) = cellstr(' ');%cellstr(YY_Annots.entrez_gene(YY_IDs(ii),:));
intranks(:,9) = num2cell(ss);
intranks(:,10) = num2cell(abs(ss));
intranks(:,11) = cellstr(YY_Data.antibodyfull(YY_IDs(ii),:));
intranks(:,12) = cellstr(YY_Data.pathway(YY_IDs(ii),:));

intranks_r = flipud(sortrows(intranks,10)); % rank based on magnitudes
aa = find(strcmp('',intranks_r(:,4))); % if no HGNC present, replace to ENSEMBL id below
intranks_r(aa,4) = intranks_r(aa,3);
aa = find(strcmp('NA',intranks_r(:,4))); % if no HGNC present, replace to ENSEMBL id below
intranks_r(aa,4) = intranks_r(aa,3);
iinornot = find(strcmp(cellstr(intranks_r(:,4)),cellstr(intranks_r(:,7)))); % Find the self-gene associations
intranks_ii = intranks_r(iinornot,:);
return
