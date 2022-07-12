function [uqList2,uqList1] = findNetworkGeneWeights(Network2check0,colno,valcolno,num2normlz,cutoffval)
%%%% Function to select association lists based on some criteria

%%%% Inputs:
% Network2check0: Initial network/association list  
% colno: Column order number if want to select based on ligand cases
% valcolno: Column order number if want to select based on association magnitude
% num2normlz: Association weight normalization value
% cutoffval: Magnitude cutoff value

%%%% Outputs:
% uqList2: Summarized list of associations
% uqList1: List of associations

%%
%%% Uncomment if selecting based on ligand condition
% Network2check = Network2check0(find(cell2mat(Network2check0(:,colno))>=1),:); 
% numInts = size(Network2check,1);

Network2check = Network2check0(find(cell2mat(Network2check0(:,valcolno))>cutoffval),:);
numInts = size(Network2check,1);

%% Find and split multiple gene RPPA edges, make duplicates to count for each gene
for qq = 1:numInts
    tempG = Network2check{qq,7};
    if any(isspace(tempG))
       aa = split(tempG);
       Network2check(qq,7) = aa(1);
       for ww = 2:length(aa)
           Network2check(end+1,:) = Network2check(qq,:);
           Network2check(end,7) = aa(ww);
       end
    end
end
%% 
uqList1 = unique([Network2check(:,4);Network2check(:,7)]);
lenL1 = length(uqList1);
numAnlyts = ones(lenL1,1);
for qq = 1:lenL1
    tempG = uqList1(qq,1);
    % Find the gene in RNAseq side first
    aa = find(strcmp(Network2check(:,4),tempG));
    for ww = 1:length(aa)
        uqList1(qq,2*ww) = Network2check(aa(ww),7);
        uqList1(qq,2*ww+1) = Network2check(aa(ww),valcolno);
    end
    % Next, find the gene in RPPA side 
    bb = find(strcmp(Network2check(:,7),tempG));
%     if bb
%         numAnlyts(qq) = length(unique(Network2check(bb,8)));
%     end
    counter = 1;
    for ww = length(aa)+1:length(aa)+length(bb)
        if ~strcmp(Network2check(bb(counter),4),tempG)
            uqList1(qq,2*ww) = Network2check(bb(counter),4);
            uqList1(qq,2*ww+1) = Network2check(bb(counter),valcolno);
        end
        counter = counter + 1;
    end
end

%% Now, sum edge weights for each gene
uqList2 = uqList1(:,1);
lenL2 = length(uqList2);
for qq = 1:lenL2
    uqList2(qq,2) = num2cell(sum(cell2mat(uqList1(qq,3:2:end)))./(num2normlz)./numAnlyts(qq,1));
end
return