function [Xi,zi, si, ni, Ji, Jmax, X0] = subsample_data(X,z,DX);
% [Xi,zi, si, ni, Ji, Jmax, X0] = subsample_data(X,z,DX);
%
% interpolate data into regular sample bins using boxcar window
%
% Input
%   X, an NxM set of coordinates
%   z, an Nx1 array of observations
%   DX, an Mx1 array of scale parameters, indicating the step size in each dimension
%
% Output
%   Xi, the mean position of the data in this cell
%   zi, the mean value at each interp. cell
%   si, the standard error (=std. dev./sqrt(n)) (or, if ni<3, insert average value)
%   ni, the number of observations going into this cell
%   Ji, the array of indices into each cell
%   Jmax, the array of number of cells in each dimension
%   X0, the location of the first grid point
%      % to make a quick grid of the data:
%      [XX,YY] = meshgrid([1:Jmax(1)]*DX(1)+X0(1),[1:Jmax(2)]*DX(2)+X0(2));
%      ZZ = repmat(nan,Jmax(1),Jmax(2)); % careful, read in flipped
%      EE = ZZ;
%      ZZ(Ji) = zi; ZZ = ZZ'; % flip to usual orientation for matlab
%      EE(Ji) = si; EE = EE'; 
%      pcolor(XX,YY,ZZ);

[N,M] = size(X);
DX = DX(:)'; % need row vector

% map data to scaled points
% J = 1,1...,1 is location X0(1,1,...,1)
X0 = floor(min(X)./DX).*DX; % make nice integer values
J = round(1+(X-repmat(X0,N,1))./repmat(DX,N,1));

% map these to index into array of unique indices
Jmax = max(J);
Ji=ones(N,1);
Jprod = 1;
for i=1:M
   Ji = Ji+(J(:,i)-1)*Jprod;
   Jprod = Jprod*Jmax(i);
end

% initialize output arrays that are as large as the largest index
Ni = max(Ji);
zi = repmat(0,Ni,1);
si = zi;
ni = zi;
Xi = repmat(0,Ni,M);

% insert values
for i=1:N
   ni(Ji(i)) = ni(Ji(i))+1;
   zi(Ji(i)) = zi(Ji(i))+z(i);
   si(Ji(i))= si(Ji(i))+z(i).^2;
   
   % keep track of cell locations used
   % Xi(Ji(i),:) = (J(i,:)-1).*DX + X0;
   
   % keep track of mean data location within cell
   Xi(Ji(i),:) = Xi(Ji(i),:) + X(i,:);
end

% get mean values at each cell
Ji = find(ni>0);
Ni = length(Ji);
ni = ni(Ji);
zi = zi(Ji)./ni;
% Xi = Xi(Ji,:); % use this for cell location center
Xi = Xi(Ji,:)./repmat(ni, 1,M); % use this for mean data location in cell

% standard deviation
%si = real( sqrt( (si(Ji)./ni)-(zi.^2) ) ); % real is used to deal with roundoff error
% standard error
si = real( sqrt( (si(Ji)./ni)-(zi.^2) ) )./sqrt(ni); % real is used to deal with roundoff error

% replace bogus std
id = find(ni<3 | si==0);
if(length(id)<Ni)
   % we will pad the error estimate with the mean value
   idg = find(ni>=3);
   stot = mean(si(idg));
else
   % all of them are missing estimates, put in data std.dev.
   stot = std(z)/sqrt(N);
end
si(id) = repmat(stot, length(id),1);