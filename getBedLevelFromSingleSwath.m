function [xBed,zBed] = getBedLevelFromSingleSwath(Xswath,Zswath,bS,minZ,belowBedZ)

% GETBEDLEVELFROMSINGLESWATH(Xswath,Zswath,BS,BEGINID,NOISEIDS) determines the bed
% level line from a single swath of the 3D Sonar.
%
% INPUT
%   Xswath, [nPings nSamples] matrix of x values
%   Zswath, [nPings nSamples] matrix of z values
%   bS, [nPings nSamples] matrix of backscatter
%   minZ, bins with a Zswath > minZ are not used. Make sure that minZ is not too close
%         to the actual bed.
%   belowBedZ, bins with a Zswath < belowBedZ are considered to be below bed.
%         Used for noise floor. Make sure belowBedZ is not too close to the
%         actual bed.
%
% OUTPUT
%   xBed, x-values of bed level line
%   zBed, z-values of bed level line
%      When zBed is a NaN, xBed is a dummy value! When interpoloating to a
%      regular grid afterwards, these dummy x ensure that gaps remain gaps.
%
% During the estimation of xBed and zBed, the algorithm will first detect the
% maximum backscatter with corresponding xBed and zBed. A second-order polynomial
% is then fitted centered around the maximum to go subgrid. This fit is based on 11 or more points.
%
% v1, 24 March 2015, Gerben Ruessink, modified from getBedLevelFromSRPSScan

% number of pings and samples
[nPings,nSamples] = size(Xswath);

% intialize output
xBed = NaN(nPings,1);
zBed = NaN(nPings,1);

% noise level (below bed)
threshold = max(max(bS((Zswath < belowBedZ))));
bS(bS <= threshold) = 0;

% loop through each ping
for ping = 1:nPings
    % ids of interest
    beginId = find(Zswath(ping,:)<minZ);
    idInt = beginId(1):nSamples;
    
    if sum(bS(ping,idInt))~= 0,  % there is a signal
        
        % find maximum
        [~,id] = max(bS(ping,idInt));   % find maximum
        
        % Go subgrid by fitting 2nd polynomial over 11 points centered at
        % maximum value. This provides a slighthly smoother bed level line
        % than just taking the position of the maximum backscatter in each
        % scanline. Sometimes this width is too small, yielding an
        % unrealistic bed level. Then, width is increased by two bins until
        % a more realistic value is obtained.
        width = 5;
        while 1
            idRange = idInt(id)-width:idInt(id)+width;
            idRange(idRange < 1) = [];
            idRange(idRange > nSamples) = [];
            if length(idRange) < 3,
                continue;
            end;
            P = polyfit(idRange,bS(ping,idRange),2);
            idSubGrid = -P(2)/(2*P(1));
            if idSubGrid >= idRange(1) && idSubGrid <= idRange(end),
                break;
            end;
            width = width + 1;
            if width > 15; % too wide, bed will not be found
                break
            end
        end;
        xBed(ping) = interp1(idInt,Xswath(ping,idInt),idSubGrid);
        zBed(ping) = interp1(idInt,Zswath(ping,idInt),idSubGrid);
        
    end; % if
    
    % if no maximum was found, xBed and zBed are NaN for this ping
end;

% detect and remove outliers
Zcrit = 3.5;   % trial and error based value
while 1
    zNoTrend = zBed;
    zNoTrend(~isnan(zBed)) = poly_detrend(xBed(~isnan(zBed)),zBed(~isnan(zBed)),2) + nanmean(zBed);
    Z = (zNoTrend - nanmean(zNoTrend)) / nanstd(zNoTrend);
    if any(abs(Z) > Zcrit),
        zBed(abs(Z) > Zcrit) = NaN;
    else
        break;
    end;
end;
% 
% % x may contain NaNs; remove and adjust x-axis
% % indices of good data points; add 0 at begin and end to detect NaN
% % at first and last point
% goodPoints = find(~isnan([0; xBed; 0]));
% 
% % do we need to do anything?
% if (length(goodPoints)-2) < length(xBed),  % -2 because 0 was added twice above
%     
%     % length of NaN parts
%     d = diff(goodPoints);
%     d = d - 1;
%     lengthGaps = d(d > 0);
%     if exist('idGaps','var'),
%         clear('idGaps');
%     end;
%     idGaps(:,1) = goodPoints(d > 0); % not +1 because 0 was added above before correctedV;
%     idGaps(:,2) = idGaps(:,1) + lengthGaps - 1;
%     for gap = 1:size(idGaps,1);
%         
%         if idGaps(gap,1) == 1,  % er zit een gat aan het begin van x
%             
%             %            % zet dummy x-waarden aan het begin van x. Dit komt
%             %            % later bij het interpoleren naar een regelmatige x-as
%             %            % weer goed.
%             %            stapgrootte = (100-x(idGaps(gap,2)+1))/(idGaps(gap,2)-idGaps(gap,1)+1);
%             %            xDummy = interp1([100 x(idGaps(gap,2)+1)],[100 x(idGaps(gap,2)+1)],...
%             %                             [100:-stapgrootte:x(idGaps(gap,2)+1)]');
%             %            x(idGaps(gap,1):idGaps(gap,2)) = xDummy(1:end-1);
%             
%             continue;
%             
%         elseif idGaps(gap,2) == length(xBed),  % er zit een gat aan het eind van x
%             
%             %             stapgrootte = (-100-x(idGaps(gap,1)-1)) / (idGaps(gap,2)-idGaps(gap,1)+1);
%             %             xDummy = interp1([x(idGaps(gap,1)-1) -100], [x(idGaps(gap,1)-1) -100],...
%             %                              [x(idGaps(gap,1)-1):stapgrootte:-100]');
%             %             x(idGaps(gap,1):idGaps(gap,2)) = xDummy(2:end);
%             
%             continue;
%             
%         else % het gat zit ergens in het midden
%             
%             stapgrootte = (xBed(idGaps(gap,1)-1) - xBed(idGaps(gap,2)+1)) / (idGaps(gap,2)-idGaps(gap,1)+2);
%             xDummy = interp1([xBed(idGaps(gap,1)-1) xBed(idGaps(gap,2)+1)],[xBed(idGaps(gap,1)-1) xBed(idGaps(gap,2)+1)],...
%                 [xBed(idGaps(gap,1)-1):-stapgrootte:xBed(idGaps(gap,2)+1)]');
%             xBed(idGaps(gap,1):idGaps(gap,2)) = xDummy(2:end-1);
%             
%         end;
%         
%     end;
%     
% end;
% 
% % remove NaNs in x (i.e., gaps at begin and/or end)
% % id = find(isnan(x));
% % x(id) = []; z(id) = [];
% 
% % remove identical x-values
% [xBed,ID] = unique(xBed);
% zBed = zBed(ID);
% 
% % sort x
% [xBed,id] = sort(xBed);
% zBed = zBed(id);
% 
% % % remove outliers (again)
% % % Zcrit = 3.5; % trial and error based value
% % iter = 0; maxIter = 3;
% % while 1
% %    iter = iter+1;
% %    q = gradient(z,x);
% %    qNoTrend = q;
% %    qNoTrend(~isnan(q)) = poly_detrend(x(~isnan(q)),q(~isnan(q)),2) + nanmean(q);
% %    Z = (qNoTrend - nanmean(qNoTrend)) / nanstd(qNoTrend);
% %    id = find(Z > Zcrit);
% %    if isempty(id) || iter == maxIter,
% %        break;
% %    end;
% %    for i = 1:length(id)
% %        if id(i)+2 > length(q),
% %            z(id(i)) = NaN;
% %            continue;
% %        end;
% %        if abs(Z(id(i)+2)) > Zcrit,
% %           z(id(i)+1) = NaN;
% %        else
% %           z(id(i)) = NaN;
% %        end;
% %    end;
% % end;

% ready
return
