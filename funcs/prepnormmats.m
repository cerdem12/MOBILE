function [Xn] = prepnormmats(X,rc,B)
%%%% Function to process/normalize columns and rows of matrices
%%%% By Cemal Erdem, Ph.D. (June 2022)

%%%% Inputs
% X: Matrix to process normalization for 
% rc: Condition flag to define the normalization type
% B: Parameter to set matrix length

%%%% Output
% Xn: Processes/normalized matrix

rowsnum = size(X,1);
colmnsnum = size(X,2);
Xn = X;
if rc == 1 % row center
    Xr = Xn;
    for qq = 1:rowsnum
        muj = mean(Xr(qq,:));
        Xr(qq,:) = (Xr(qq,:) - muj);
    end
    Xn = Xr;
elseif  rc == 2 % row len = 1
    Xr = Xn;
    for qq = 1:rowsnum
        norm2j = norm(Xr(qq,:));
        Xr(qq,:) = (Xr(qq,:))./norm2j;
    end
    Xn = Xr;
elseif  rc == 3 % column center
    Xc = Xn;
    for qq = 1:colmnsnum
        muj = mean(Xc(:,qq));
        Xc(:,qq) = (Xc(:,qq) - muj);
    end
    Xn = Xc;
elseif  rc == 4 % column len = 1
    Xc = Xn;
    for qq = 1:colmnsnum
        norm2j = norm(Xc(:,qq));
        Xc(:,qq) = (Xc(:,qq))./norm2j;
    end
    Xn = Xc;
elseif  rc == 5 % mat len = A
    Xc = Xn;
    A0 = sqrt(sum(sum(Xc.^2)));
    Xn = sqrt(B).*((Xc)./(A0)); 
elseif  rc == 6 % column-center, row-center, row-len=1
    Xc = Xn;
    for qq = 1:colmnsnum
        muj = mean(Xc(:,qq));
        Xc(:,qq) = (Xc(:,qq) - muj);
    end
    Xn = Xc;
    Xr = Xn;
    for qq = 1:rowsnum
        muj = mean(Xr(qq,:));
        Xr(qq,:) = (Xr(qq,:) - muj);
        norm2j = norm(Xr(qq,:));
        Xr(qq,:) = (Xr(qq,:))./norm2j;
    end
    Xn = Xr;    
elseif  rc == 7 % column-center, row-center, row-std=1
    Xc = Xn;
    for qq = 1:colmnsnum
        muj = mean(Xc(:,qq));
        Xc(:,qq) = (Xc(:,qq) - muj);
    end
    Xn = Xc;
    Xr = Xn;
    for qq = 1:rowsnum
        muj = mean(Xr(qq,:));
        Xr(qq,:) = (Xr(qq,:) - muj);
        stdj = std(Xr(qq,:),1);
        Xr(qq,:) = (Xr(qq,:))./stdj;
    end
    Xn = Xr;
end
return
