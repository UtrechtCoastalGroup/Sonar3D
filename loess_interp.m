function [zi, rmse, bi, ci, sk, Ni, mso_dev] = ...
   loess_interp(x,z,sigma_z,xi,lx, model_order);
% filter/interpolate data z(x) to xi using filter characterized by lx; 
% use method of Schlax and Chelton, 1992
%
% [zi, rmse, bi, ci, sk, Ni, mso_dev] = loess_interp(x,z,sigma_z,xi,lx, model_order);
% Input
%  x, the nxm observation location matrix, m is the number of Dimenisions
%  z, the nx1 observation array
%  sigma_z are the expected rms of observation errors on z
%     sigma_z could be constant (equal weight) or variable
%  xi, the pxm interpolation location matrix
%  lx are correlation scales, dx1
%      ok, I did synthetic tests with white spectrum, which suggest that:
%      -Linear    filter: (3 dx) < lx < (Lpreserve/5)
%      -Quadratic filter: (6 dx) < lx < (Lpreserve/2)
%         where dx is sample spacing ...
%         and Lpreserve is wavelength that suffers less than 10% loss in variance
%   model_order=1 is linear in all parameters, porder=2 is quadratic, ...
%      model_order>2 is not supported at the moment
%
% Output
%    zi, are the estimated values at xi
%    rmse, are the estimated rms errors for the estimate
%    bi, are the estimated model coefficients
%    ci, are the estimated model coefficients rmse
%    sk, are the skill estimates of the weighted regression
%    Ni, are the number of effective dof each time
%    mso_dev, is the variance og detrended the input data

% measure input
[n,m] = size(x);
[ni, mi] = size(xi);
if (mi~=m)
   fprintf('dimension of location matrices differ\n');
   return;
end

% fix up sigma_z, if constant
[ns,ms]=size(sigma_z);
if(max([ns,ms])==1)
   sigma_z = repmat(sigma_z, n,1);
else
   sigma_z = sigma_z(:);
end

% deal with nans
id = find(isfinite(sum([z,x,sigma_z]')));
n = length(id);
x=x(id,:);
z = z(id);
sigma_z = sigma_z(id);

% find out if any parameters are constant
xmean = mean(x);
xsdev = std(x);
idp = find(xsdev>0);
d=length(idp);

% center and normalize field area
x = x(:,idp)-repmat(xmean(idp),n,1);
x = x(:,idp)./repmat(xsdev(idp),n,1);
xi = xi(:,idp)-repmat(xmean(idp),ni,1);
xi = xi(:,idp)./repmat(xsdev(idp),ni,1);
lx = lx(:);
lx = lx(idp)./xsdev(idp)';

% remove global linear (for stability) trend
zmean = mean(z);
zsdev = std(z);
z = z - zmean;
XX = (x'*x)/n;
XZ = (z'*x)/n;
XX_inv = inv(XX);
btrend = XX_inv*XZ'; 
ztrend = x*btrend;
zi_trend = zmean + xi*btrend;

% compute deviations from trend
zdev = z - ztrend;

% compute global variance
mso_dev = zdev'*zdev/n;
fprintf('standard deviation of obs. is %.3f\n', sqrt(mso_dev))

% normalize the obs. errors with the total variance
% rationale is errors near zero are mapped to some constant, and spatial interpolation proceeds as expected
% large errors lead to weights that increase as 1/sigma_z
sigma_z = sqrt(1+(sigma_z.^2)/mso_dev); 

% pass each xi to model
zi = repmat(nan,ni,1);
rmse = zi;
sk = zi;
Ni = zeros(ni,1);
bi = repmat(nan,ni,m);
ci = repmat(nan,ni,m);
kx = 1./lx(:);
for i = 1:ni
    % center on xi
    t = x-repmat(xi(i,:),n,1);
    
    % don't send all data
    id = find((abs(t)*kx)<d);
    if(length(id)>m)
        % interpolate one
        [zi(i), rmse(i), b, c, sk(i), Ni(i)] = ...
            loess_filt(t(id,:),zdev(id),sigma_z(id),lx, model_order);
        bi(i,idp) = (b'+ btrend')./xsdev(idp);
        ci(i,idp) = c'./xsdev(idp);
    end
end
zi = zi + zi_trend;