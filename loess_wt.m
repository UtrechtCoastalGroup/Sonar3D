function x = loess_wt(x,lx)
%
% [w] = loess_wt(x,lx)
%
% Input
%   x are nxm inputs, weights will be centered on x=0
%   lx is mx1 length scales
%
% Output
%   w are weights 0<=w<=1

[n,m] = size(x);
lx = lx(:);

% normalize by correlation scale
x = x./repmat(lx',n,1);

% weights for least-sq. minimization
x = x.^2;

% sum accross all inputs, if m>1
if(m>1)
   x = x*ones(m,1);
end

% the loess weighting function (zero if x>1)
x = ((1-(x.^3)).^3).*(x<1);  % used by greenslade
% x = ((1-(x.^1.5)).^3).*(x<1); % used by schlax and chelton
