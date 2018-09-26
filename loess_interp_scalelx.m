function [zi, ei] = loess_interp_scalelx(x,z,sigma_z,xi,lx_target);
% filter/interpolate data z(x) to xi using filter characterized by lx; 
% allow for scale to depend on sample spacing
%
% [zi, ei] = loess_interp_scalelx(x,z,sigma_z,xi,lx_target);
%
% Input
%  x, the nxm observation location matrix, m is the number of Dimenisions
%  z, the nx1 observation array
%  sigma_z are the expected rms of observation errors on z
%     sigma_z could be constant (equal weight) or variable
%  xi, the pxm interpolation location matrix
%  lx_target are correlation scales, qxm,
%     where each of the q ROWS are different smoothing scales 
%     to be applied sequentially, e.g., [lxbig, lybig; lxsmall, lysmall; ...]
%
% Output
%    zi, are the estimated values at xi
%    ei, are the estimated rms errors for the estimate
model_order = 2;

% check input
[n,m] = size(x);
[ni, mi] = size(xi);
if (mi~=m)
    fprintf('dimension of location matrices differ\n');
    break;
end
if (m>3)
    fprintf('cannot do D>3, sorry\n');
    break;
end
if (length(sigma_z)==1)
    sigma_z = repmat(sigma_z,n,1);
end

% set final scales
lx = lx_target(end,:);
dx = lx/4;

% make a common grid
Mingrid = min([xi;x]); % grid must span data and output grid!!!
Maxgrid = max([xi;x]);

% assume 1-d for now
xg = [Mingrid(1)-dx(1), Mingrid(1):dx(1):Maxgrid(1), Maxgrid(1)+dx(1)];
Nx = length(xg);

% find out which cells have data
idx = 1+fix((x(:,1)-xg(1))/dx(1));

% expand to n-d
if(m>1)
    yg = [Mingrid(2)-dx(2), Mingrid(2):dx(2):Maxgrid(2), Maxgrid(2)+dx(2)];
    [Xg,Yg] = meshgrid(xg,yg);
    Ny = length(yg);
    idy = 1+fix((x(:,2)-yg(1))/dx(2));
end
if(m>2)
    tg = [Mingrid(3)-dx(3), Mingrid(3):dx(3):Maxgrid(3), Maxgrid(3)+dx(3)];
    [Xg,Yg,Tg] = meshgrid(xg,yg,tg);
    Nt = length(tg);
    idt = 1+fix((x(:,3)-tg(1))/dx(3));
end

% make a binary index
data_true = 0*Xg;
if(m==1)
    data_true(idx) = ones(idx);
elseif(m==2)
    data_true(idy+Ny*(idx-1))=ones(n,1);
elseif(m==3)
    data_true(idy+Ny*(idx-1)+(idt-1)*Ny*Nx)=ones(n,1);
end

% scroll through scales
cnt = length(lx_target(:,1));

% initialize
Zg = nan*zeros(size(Xg));
Eg = Zg;
if(m==1)
    Ag = Xg(:);
elseif(m==2)
    Ag = [Xg(:), Yg(:)];
elseif(m==3)
    Ag = [Xg(:), Yg(:), Tg(:)];
end

% store the final output here
zi = zeros(ni,1);
ei = zi;

% begin
for i=1:cnt
    % use appropriate grid spacing
    di = max([fix((lx_target(i,:)/4)./dx); ones(1,m)]);
    idx = [1,2:di(1):(Nx-2),(Nx-1),Nx]';
    Nix = length(idx);
    if(m>1)
        idy = [1,2:di(2):(Ny-2),(Ny-1),Ny]';
        Niy = length(idy);
    end     
    if(m>2)
        idt = [1,2:di(3):(Nt-2),(Nt-1),Nt]';
        Nit = length(idt);
    end
    
    % avoid oversampling of raw data
    [xs,zs, sigma_z] = subsample_data(x,z,(di.*dx)/4);
    
    % tack on gridded data (must include endpoints)
    if(m==1)
        Ag = col(Xg(idx));
        Zg = col(Zg(idx));
        Eg = col(Eg(idx));
        id = find(isfinite(Zg) & ~col(data_true(idx)));   
    elseif(m==2)
        Ag = [col(Xg(idy, idx)), col(Yg(idy, idx))];;
        Zg = col(Zg(idy, idx));
        Eg = col(Eg(idy, idx));
        id = find(isfinite(Zg) & ~col(data_true(idy,idx)));   
    elseif(m==3)
        Ag = [col(Xg(idy, idx, idt)), col(Yg(idy, idx, idt)), col(Tg(idy, idx, idt))];;
        Zg = col(Zg(idy, idx, idt));
        Eg = col(Eg(idy, idx, idt));
        id = find(isfinite(Zg) & ~col(data_true(idy,idx,idt)));   
    end
    xdata = [xs; Ag(id,:)];
    zdata = [zs; Zg(id)];
    edata = [sigma_z; Eg(id)];     
    
    % do it
    tic
    [Zg,Eg,bis,cis,skis,nis,mso]=loess_interp(xdata,zdata,edata,Ag,lx_target(i,:),model_order);
    toc
    if(m==2)
        Zg = reshape(Zg, Niy,Nix);
        Eg = reshape(Eg, Niy,Nix);
    elseif(m==3)
        Zg = reshape(Zg, Niy,Nix,Nit);
        Eg = reshape(Eg, Niy,Nix,Nit);
    end
    
    % watch   
    if(m>1)
        figure
        subplot 211
        imagesc(xg(idx), yg(idy), Zg(:,:,1)); axis xy
        shading flat
        caxis([-1 1]*sqrt(mso))
        colorbar
        title(num2str(lx_target(i,:)))
        subplot 212
        imagesc(xg(idx), yg(idy), Eg(:,:,1)); axis xy
        hold on
        plot(xdata(:,1), xdata(:,2),'.m')
        hold off
        shading flat
        caxis([0 1]*sqrt(mso))
        colorbar
        drawnow
    end
    
    % add this field to summation of z,e
    intmethod = 'linear';
    if(m==1)
        zi = zi + interp1(xg(idx), Zg, xi(:,1), intmethod);
        ei = ei + interp1(xg(idx), Eg, xi(:,1), intmethod);
        zfilt = interp1(xg(idx), Zg, x(:,1), intmethod);
    elseif(m==2)
        zi = zi + interp2(xg(idx),yg(idy), Zg, xi(:,1),xi(:,2), intmethod);
        
        % the interpolation error variance is cumulative
        ei = ei + interp2(xg(idx),yg(idy), Eg, xi(:,1),xi(:,2), intmethod);
        
        % need this to get deviations
        zfilt = interp2(xg(idx),yg(idy), Zg, x(:,1),x(:,2), intmethod);
        
    elseif(m==3)
        zi = zi + interp3(xg(idx),yg(idy),tg(idt), Zg, xi(:,1),xi(:,2),xi(:,3), intmethod);
        ei = ei + interp3(xg(idx),yg(idy),tg(idt), Eg, xi(:,1),xi(:,2),xi(:,3), intmethod);
        zfilt =   interp3(xg(idx),yg(idy),tg(idt), Zg, x(:,1),x(:,2),x(:,3), intmethod);
    end
    
    % get deviations
    z = z - zfilt;
    
    % variance removed:
    msd = max([0, mso - var(zfilt)]);
    
    % weight the gridded data by interp error plus removed scales
    Eg = sqrt(Eg.^2 + msd);
    
    % and fill out the entire grid for to supplement data
    Zg = 0*Xg; % deviations are what we are gridding  
    if(m==1)
        Eg = interp1(xg(idx),Eg, Xg, intmethod);
    elseif(m==2)
        Eg = interp2(xg(idx),yg(idy), Eg, Xg,Yg, intmethod);
    elseif(m==3)
        Eg = interp3(xg(idx),yg(idy), tg(idt), Eg, Xg,Yg,Tg, intmethod);
    end
    drawnow;
end

