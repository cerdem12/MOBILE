function [TTg,BBg,a0Vals,YYpg,dsedevg2] = glmnetRunner(XX,YY,toplot,CVnfolds)
%%%% Function to run customized glmnet simulations
%%%% By Cemal Erdem, Ph.D. (June 2022)

%%%% Inputs
% XX: RHS of the equation (RNAseq or ATACseq data)
% YY: LHS of the equation (RPPA or RNAseq data)
% toplot: Flag for plotting (1=plot, 0=no plot)
% CVnfolds: Number of cross-validation folds in glmnet

%%%% Outputs
% TTg: Lasso coefficients matrix instance 
% BBg: glmnet package output structure
% a0Vals: Array of y-axis intercepts, usually an all-zero array (val<1e-6)
% YYpg: Calculated LHS array
% dsedevg2: Inference error value

% glmnet options
nosamples = size(YY,1);
opts = glmnetSet;
opts.standardize = 0;
options = glmnetSet(opts);
count2 = 0;

if nargin < 4 % Set cross-validation fold number if not provided as input
    CVnfolds = 3;  
end

for qq = 1:nosamples
    %%%% cvglmnet default inputs:(x,y,family,options,type,nfolds,foldid,parallel,keep,grouped)
    [BBg] = cvglmnet(XX',transpose(YY(qq,:)),'gaussian',options,'default',CVnfolds,[],1);
    Lming = BBg.lambda_min; % Find lowest prediction error lambda value
    Lindexg = find(BBg.lambda==Lming); % Find index of lambda_min
    Bmin = BBg.glmnet_fit.beta(:,Lindexg); % Retrieve the Lasso coefficients array
    if any(Bmin) % If the coefficients array is not empty, finish the iteration
        TTg(qq,:) = transpose(Bmin); % A new row of Lasso coefficients matrix is found
        a0Vals(qq,1) = BBg.glmnet_fit.a0(Lindexg); % the y-intercept value
    else % if the beta array from lambda_min is empty:
        Lming2 = BBg.lambda_1se; % get the lambda value/index pair for which fitting error is within 1 std-error
        Lindexg2 = find(BBg.lambda==Lming2);
        Bmin2 = BBg.glmnet_fit.beta(:,Lindexg2); % Coeff. array at 1 std-err lambda
        if Bmin~=Bmin2 % if the array is not empty:
            TTg(qq,:) = transpose(Bmin);
            a0Vals(qq,1) = BBg.glmnet_fit.a0(Lindexg2);
        else % finally, if the second beta array is also empty, find the beta array with most consistent number of non-zero coefficients
            count2 = count2 + 1;
            DFuniqs = unique(BBg.glmnet_fit.df);
            DFreps = histc(BBg.glmnet_fit.df,DFuniqs);
            [~, Dfmaxid] = max(DFreps);
            Dfmaxid2 = DFuniqs(Dfmaxid);
            Dfmaxid3 = find(BBg.glmnet_fit.df==Dfmaxid2);
            Bmaxid = Dfmaxid3(end); % or Dfmaxid(1);
            Bmin3 = BBg.glmnet_fit.beta(:,Bmaxid); % Latest beta coefficient array candidate
            TTg(qq,:) = transpose(Bmin3); % Lasso coefficients matrix row 
            a0Vals(qq,1) = BBg.glmnet_fit.a0(Bmaxid); % the y-intercept value
        end
    end
end 
devg2 = TTg*XX-YY; % Find the prediction error
dsedevg2 = (sum(sum(devg2.^2))); % Prediciton error value
YYpg = TTg*XX; % Fitting
if toplot % If plotting flag is ON:
    figure; hold on;
    plot(1:nosamples,YY(:,1),'b','LineWidth',2)
    plot(1:nosamples,YYpg(:,1),'r:','LineWidth',2)
    % ylim([-5 5])
    xlim([0 nosamples+2])
end
return