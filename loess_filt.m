function [zo, rmse, b, c, sk, nu] = loess_filt(X,z,sigma_z,lx, porder);
% return loess filter estimate at t=0
%
% [zo, rmse, b, c, sk, nu] = loess_filt(X,z,sigma_z,lx, porder);
%
% Input
%   X are independent variables, centered at prediction location
%        is nxd, where d is dimension
%        do NOT include constant in any collumns of X
%   z are observations, nx1
%      IMPORTANT, z are assumed to be somehow detrended.  
%      Large mean values over the width, lx, confuse the error analysis
%   sigma_z are the expected observation errors on z
%      if errors unknown, then sigma_z=0 should be used
%   lx are correlation scales, dx1
%      ok, I did synthetic tests with white spectrum:
%      -Linear    filter: (3 dx) < lx < (Lpreserve/5)
%      -Quadratic filter: (6 dx) < lx < (Lpreserve/2)
%      where dx is sample spacing 
%      and Lpreserve suffers less than 10% loss in variance
%   porder is the order of the local polynomial, 
%      porder=1 is linear in all parameters, porder=2 is quadratic, ...
%      porder>2 is not supported at the moment
%
% Output
%    zo is the estimate at x=0
%    rmse is the estimated error
%    b are coeficients (excluding zo)
%    c are coeficients variances (excluding zo)
%    sk are model skill of regression
%    nu are the effective dof for each regression

% check the porder parameter validity
if(porder>2)
   fprintf('porder>2 is not supported\n')
   return
end

% initialize output
[n,d] = size(X);
zo = nan;
rmse = nan;
sk=nan;
nu=0;
b = repmat(nan,d,1);
c = repmat(nan,d,1);

% get weights
% Loess, in this case
w = loess_wt(X,lx(:));

% should exit here if weights are zero
if(max(w)==0)
   return
end

% weight by variance
% don't worry about magnitude of weights, that is managed in the regression code
w = w./sigma_z(:);

% build the input
m = d+1;

if(porder==1)
   % build linear input
   X=[ones(n,1),X];
elseif(porder==2)
   % get quadratic terms
   % there will be (m^2 + m)/2 entries (including the constant)
   q = 0.5*m*(m+1);
   X=[ones(n,1),X, zeros(n,q-m)];
   for i = 2:(d+1)
      % get the quadratic self-terms
      for j = 2:i
         m  = m+1;
         X(:,m) = (X(:,i).*X(:,j));
      end
   end
end

% effective dof
nu = sum(w)/max(w);

% check to see if we are underdetermined
if (fix(nu)<=m)
   return
end

% send to regression, with weights
[b,c,sk, nu, msz] = regr_xzw(X,z,w);

% put the primary output in its place
zo = b(1);
rmse = c(1);

% get the 1st order terms
b = [b(2:(d+1))];
c = [c(2:(d+1))];

% now, if mse is negative
if (imag(rmse)~=0)
   rmse = 0;
   return
end