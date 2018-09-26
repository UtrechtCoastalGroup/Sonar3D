function [Xi,Yi,Zi, Ei] = loess_grid2dh(x,y,z, xi,yi, lx,ly);
% [Xi,Yi,Zi, Ei] = loess_grid2dh(x,y,z, xi,yi, lx,ly);
%
% take advantage of some optimization steps in order to interpolate random data
% 
% Input
%   x,y,z, the data
%   xi, yi, the grid locations in x,y
%   lx,ly, the smoothing scales
%
% Output
%   Xi, Yi, the grid location matrices
%   Zi, the interpolated surface
%   Ei, the error estimates
%
% Note, interpolation model is quadratic
model_order = 2;

% first subsample the data
% subsample spacing set so that loess filter will not alias
% the filter will alias variation if dx is too large, so sample at dx=lx/4, or smaller
fprintf('subsampling the data\n')
tic
DXsub = [lx, ly]/8;
[Xbin, zbin, sbin, nbin] = subsample_data([x, y],z,DXsub);
toc

% the filter will fully remove features shorter than lx/2 (quadratic model)
% if this is the nyquist wavenumber, then we need to sample at dx=lx/4
dxi = lx/4;
dyi = ly/4;
xr=(xi(1)-dxi):dxi:(xi(end)+dxi); 
yr=(yi(1)-dyi):dyi:(yi(end)+dyi); 
[Xr, Yr] = meshgrid(xr,yr);
[Xi,Yi] = meshgrid(xi(:),yi(:));

% check that subsample is not over sampling
if(length(Xr(:))<length(Xi(:)))
    
    
    % interpolate
    fprintf('interpolating the data\n')
    tic
    [zr, rmse] = loess_interp(Xbin,zbin, sbin, [Xr(:),Yr(:)],[lx,ly], model_order);
    toc
    
    % reconstruct the subsampled grid
    Zi = reshape(zr,length(yr),length(xr));
    Ei = reshape(rmse,length(yr),length(xr));
    
    % now, do fast linear interpolation to get to desired grid
    Zi = interp2(Xr,Yr,Zi,Xi,Yi);
    Ei = interp2(Xr,Yr,Ei,Xi,Yi);
    
else
    
    fprintf('interpolating the data\n')
    tic
    [zr, rmse] = loess_interp(Xbin,zbin, sbin, [Xi(:),Yi(:)],[lx,ly], model_order);
    toc
    
    % reconstruct the subsampled grid
    Zi = reshape(zr,length(yi),length(xi));
    Ei = reshape(rmse,length(yi),length(xi));
    
end
