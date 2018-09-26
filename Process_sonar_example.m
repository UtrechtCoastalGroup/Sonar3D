% Script file to process a single data file collected by the 3D Ripple Profiling
% Logging Sonar during the DVA campaign

clc;
clear;
close all;
warning('off','all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BEGIN INVULSECTIE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showResults = 0;   % plot intermediate results
plotting = 1;      % plot final results?

outputPath = 'd:\SONAR analyse\AZG\Frame 4\';
campaign = 'AZG';
frame = 'F4P5';

% location of the data
workPath = 'w:\SeaWad\KG2Database\trunk\morphology\rippleprofiler\raw\2017_09_ameland\seawad\frame04\';
mapinfo = dir([workPath,'*2017']);
dates = size(mapinfo,1);

% instrument data, copy from metadata file
SONAR3D(1).timeIN = datenum([2017,8,29,17,55,0]);     % net voor de eerste inzet
SONAR3D(1).timeOUT = datenum([2017,9,13,13,29,59]);   % net voor servicen
SONAR3D(1).zeroOrientation = 283;                     % = offset x-axis with respect to AQD + angle between AQD and North
SONAR3D(1).z = 0.977;                                 % Bottom head center ==> not sampling volume yet!
SONAR3D(1).depth = 5;                                 % approximate depth
SONAR3D(1).location = [165276, 611043];               % RDx,RDy

SONAR3D(2).timeIN = datenum([2017,9,13,13,30,0]);     % net na servicen
SONAR3D(2).timeOUT = datenum([2017,9,19,11,46,59]);   % net na recovery
SONAR3D(2).zeroOrientation = 267.1;                   % = offset x-axis with respect to AQD + angle between AQD and North
SONAR3D(2).z = SONAR3D(1).z;                          % Bottom head center ==> not sampling volume yet!
SONAR3D(2).depth = SONAR3D(1).depth;                  % approximate depth
SONAR3D(2).location = SONAR3D(1).location;            % RDx,RDy

SONAR3D(3).timeIN = datenum([2017,9,19,11,47,0]);     % net na servicen
SONAR3D(3).timeOUT = datenum([2017,10,9,17,50,0]);    % net na recovery
SONAR3D(3).zeroOrientation = 321;                     % = offset x-axis with respect to AQD + angle between AQD and North
SONAR3D(3).z = SONAR3D(1).z;                          % Bottom head center ==> not sampling volume yet!
SONAR3D(3).depth = SONAR3D(1).depth;                  % approximate depth
SONAR3D(3).location = SONAR3D(1).location;            % RDx,RDy

smallscalefilter = 0.05;
largescalefilter = 0.10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EIND INVULSECTIE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
tel = 0;
% read the data
for date = 1:dates;
    fileInfo = dir([workPath, mapinfo(date).name,'\','*.RW2']);
    nRawFiles = size(fileInfo,1);
    
    for h = 1:nRawFiles
        tel = tel+1;
        fileName = fileInfo(h).name;
        fprintf(1,'Filename: %s\n',fileInfo(h).name);
        year = mapinfo(date).name(5:8);
        month = mapinfo(date).name(3:4);
        day = mapinfo(date).name(1:2);
        hour = fileName(1:4);
        [header, data, range] = readSonar3DRW2([workPath,mapinfo(date).name,'\',fileName]);
        if datenum(header.when)<SONAR3D(1).timeIN||datenum(header.when)>SONAR3D(end).timeOUT;
            continue
        end
        
        timevec(tel,:) = header.when;
        % The header contains relevant information on how the 3D Sonar was
        % programmed, for example, the arcwidth centered around the vertical axis,
        % the swathstep, and rotationstep. Part of the information within the
        % header has already been used to compute the range.
        
        % The data are the signal amplitudes organized as [Nsamples Npings
        % NSwaths]. That is, a single swath contains Npings (profiles), each from
        % range(1) to range(end), with each profile comprising Nsamples. In total
        % NSwaths were performed to cover a full circle. Nsamples equals the length
        % of range, Npings equals header.arc.
        [Nsamples, Npings, Nswaths] = size(data);
        
        % The unit stepsize in a swath and a rotation step is 1 in the 400-degree
        % system
        stepSize = 1*360/400;
        
        % The processing is carried out in the following steps:
        % (1) Determine bed level for each swath
        % (2) Put into (x,y,z)
        % (3) Rotate to have y north and x east
        % (4) Outlier detection, probably redundant
        % (5) Interpolate to a regular grid, and remove trend to reveal ripples as
        %     perturbations
        
        
        %% STEP 1
        
        % We are going to make a "virtual" x-z axis for each swath, irrespective of
        % the rotation angle. This is going to result in a bed profile from each
        % swath. Unfortunately, the manual is not very specific about the angles
        % within the swath, other than that the arc is centered around the
        % vertical. So, 0 should be pointing straight down. The Npings is even,
        % which suggests that there are Npings/2 < 0 and Npings/2 > 0. Curiously,
        % for the Bardex II data, the viewer that comes with the sonar produces a
        % marked difference between the first and last swath in terms of bed
        % elevation. This is not due to some overall roll or tilt, but indicates
        % 'uncertainties' with how the angles within the arc are defined. This
        % difference between first and last swath is not found for occasions when
        % the sonar is externally triggered. In Bardex II, it was triggered on-line
        % from a laptop. Could that matter?
        
        % Interestingly, we do not know whether the Sonar3D samples from negative
        % to positive angles, or from positive to negative. We will assume here
        % that it does so from negative to positive. In the swath with rotation =
        % 0, the positive x goes out from the 0 rotation angle of the sonar head.
        % Or, said for Bardex II, the 0 pointed approximately to wave paddle. So,
        % positive x for rotation = 0 is pointing toward the wave paddle and
        % negative x to the swash zone.
        
        slotoffset = 0;
        THdeg = (-header.Arc/2:header.SwathStep:header.Arc/2-1).*stepSize + slotoffset*stepSize;
        THrad = deg2rad(THdeg);
        
        % Note: for externally triggered systems it appears that slotoffset should be 0.
        % And curiously, some Bardex II scans required step to be -2. For AZG frame
        % 5, slotoffset was -1. For DVA1, several values were tested and 0 is correct.
        
        % Turn range and TH into Cartesian coordinates. This is a coordinate scheme
        % within this swath, with X the horizontal component and Z the vertical
        % component.
        THETA = ones(length(range),1)*THrad;
        RHO = ones(length(THrad),1)*range;
        [Zswath,Xswath] = pol2cart(THETA',RHO);       % Convert swaths from polar to cartesian
        
        % Note Zswath and Xswath are now Npings x Nsamples, which is different from the order
        % in data. The bed level detection algorithm called below assumes Npings x
        % Nsamples. There is a transpose there for each swath in data.
        %
        % And, pol2cart has Zs,Xs as output as THETA is defined with respect to the
        % x-axis, which is here vertically down, and hence Z.
        
        % change sign of Zswath to have below transducer head as negative
        Zswath = -Zswath;
        
        % Prepare output. For every ping we are going to find a bed level at (xBed,zBed)
        xBed = NaN(Npings,Nswaths);
        zBed = xBed;
        
        % The bed is expected between minZ and belowBedZ. The signal beneath belowBedZ is used to quantify noise.
        % At least for the UU frame, the height of the Sonar was ~0.9 m, so we take
        % minZ as sensorheight/2 and belowBedZ as sensorheight*2. However, if the
        % bed level changes are large, this might need to be changed!
        minZ = -0.45;
        belowBedZ = -1.8;
        for i = 1:Nswaths          % for every swath
            
            bS = data(:,:,i)';     % input signal amplitudes, with transpose to go Npings x Nsamples
            [xBed(:,i), zBed(:,i)] = getBedLevelFromSingleSwath(Xswath,Zswath,bS,minZ,belowBedZ);  % get bed level estimates for this swath
            
            if showResults,   % if 1, show intermediate results
                figure(1);
                pcolor(Xswath,Zswath,bS);
                shading flat;
                hold on;
                plot(xBed(:,i),zBed(:,i),'*-k');
                hold off
                title(['Swath ' num2str(i) ' of ' num2str(Nswaths)])
                pause(0.1);
            end;
            
        end
        
        % figure();
        % plot(xBed,zBed)
        % xlabel('x (m)')
        % ylabel('z (m)')
        %% STEP 2
        
        % Put all xBed and zBed into a proper (x,y,z) matrix
        THrot = (0:header.RotateStep:200-header.RotateStep)*stepSize;
        THrotrad = deg2rad(THrot);
        
        THETArot = ones(Npings,1)*THrotrad;
        [X,Y] = pol2cart(THETArot,xBed);
        
        % Turn these into columns and remove NaNs
        XCol = X(:);
        YCol = Y(:);
        ZCol = zBed(:);
        nandat = isnan(XCol)|isnan(YCol)|isnan(ZCol);
        XCol(nandat) = [];
        YCol(nandat) = [];
        ZCol(nandat) = [];
        
        % figure();
        % scatter(YCol,XCol,2,ZCol);
        % xlabel('x (m)')
        % ylabel('y (m)')
        % title('Depth below sensor')
        % hcb=colorbar;
        % title(hcb,'(m)')
        %% STEP 3
        
        % rotate X and Y such that X becomes East-West axis and Y becomes
        % South-North axis
        if datenum(header.when)<SONAR3D(2).timeIN   % voor het servicen
            rotation = SONAR3D(1).zeroOrientation;
        elseif length(SONAR3D)>2 && datenum(header.when)>SONAR3D(end).timeIN                                        % na het servicen
            rotation = SONAR3D(3).zeroOrientation;
        else
            rotation = SONAR3D(2).zeroOrientation;
        end
        
        phiRad = deg2rad(rotation);
        coord = [XCol, YCol];
        rot = [cos(phiRad) -sin(phiRad); sin(phiRad) cos(phiRad)];
        coordRot = (rot*coord')';
        XCol = coordRot(:,1);
        YCol = coordRot(:,2);
        %
        %% STEP 4
        
        % Outliers are removed based on a deviation from the median. This is
        % presumably a redundant step, especially when only scans with no-wave
        % action are analyzed. No harm done because then number of outliers
        % is 0.
        zOutRemoved = NaN(size(ZCol));
        for ix = -3.5:0.25:3.5
            for iy = -3.5:0.25:3.5
                windowloc = find(XCol>ix & XCol<=ix+0.25 & YCol>iy & YCol<=iy+0.25);
                zwindow = ZCol(windowloc);
                n = removeOutlierMedian(zwindow,0.07,NaN);
                zOutRemoved(windowloc) = n.data;
            end
        end
        fprintf('%d outliers removed\n',sum(sum(isnan(zOutRemoved))));
        ZCol = zOutRemoved;
        XCol(isnan(ZCol)) = [];
        YCol(isnan(ZCol)) = [];
        ZCol(isnan(ZCol)) = [];
        
        % figure();
        % scatter(YCol,XCol,2,ZCol,'filled')
        % xlabel('x (m)')
        % ylabel('y (m)')
        % title('Depth below sensor')
        % hcb=colorbar;
        % title(hcb,'(m)')
        
        %% STEP 5
        
        % Interpolate to a regular grid, from -2.5 to 2.5 in both x and y with a
        % 0.01 m step size. Use lx and ly to fill gaps; note that this larger than
        % for the perturbations from the central line. Remove zGrid values when the error
        % exceeds 0.01 m. This removes points that were extrapolated.
        
        xGrid = -2.5:0.01:2.5;
        yGrid = -2.5:0.01:2.5;
        
        lxs = smallscalefilter;   % small
        lys = lxs;%04;
        [xi,yi,zGrids,eGrids] = loess_grid2dh(XCol,YCol,ZCol,xGrid,yGrid,lxs,lys);
        zGrids(eGrids > 0.01) = NaN;  % remove points at the border
        
        lxl = largescalefilter;   % small
        lyl = lxl;%04;
        [xi,yi,zGridl,eGridl] = loess_grid2dh(XCol,YCol,ZCol,xGrid,yGrid,lxl,lyl);
        zGridl(eGridl > 0.01) = NaN;  % remove points at the border
        
        % figure();
        % pcolor(yi,xi,zGrids);
        % shading flat
        % title('Depth with lx=ly=0.05 m')
        % hcb=colorbar;
        % title(hcb,'(m)')
        % caxis([-1.083 -1.005])
        % xlabel('x (m)')
        % ylabel('y (m)')
        %
        % figure();
        % pcolor(yi,xi,zGridl);
        % shading flat
        % title('Depth with lx=ly=0.10 m')
        % hcb=colorbar;
        % title(hcb,'(m)')
        % caxis([-1.083 -1.005])
        % xlabel('x (m)')
        % ylabel('y (m)')
        
        %%
        % remove NaNs once more
        XCol_s = xi(:);
        YCol_s = yi(:);
        ZCol_s = zGrids(:);
        XCol_s(isnan(ZCol_s)) = [];
        YCol_s(isnan(ZCol_s)) = [];
        ZCol_s(isnan(ZCol_s)) = [];
        
        Zmean_s = mean(ZCol_s);
        Zmeans(tel) = Zmean_s;
        
        XCol_l = xi(:);
        YCol_l = yi(:);
        ZCol_l = zGridl(:);
        XCol_l(isnan(ZCol_l)) = [];
        YCol_l(isnan(ZCol_l)) = [];
        ZCol_l(isnan(ZCol_l)) = [];
        
        Zmean_l = mean(ZCol_l);
        Zmeanl(tel) = Zmean_l;
        
        zGridPert_s = zGrids - Zmean_s;
        zGridPert_l = zGridl - Zmean_l;
        
        
        %% Plotting
        if plotting
            figure();
            pcolor(yi,xi,zGridPert_s); shading('flat');
            set(gca,'fontsize',12,'fontname','times','clim',[-0.09 0.09]);
            hcb=colorbar;
            title(hcb,'[m]')
            hold on
            text(1,0.9,'1','fontname','times')
            xlabel('x (m)','fontsize',12,'fontname','times');
            ylabel('y (m)','fontsize',12,'fontname','times');
            runn = [day '-' month '-' year ' ' hour(1:2) ':' hour(3:4)];
            title(runn,'fontname','times');
            axis([-1 1 -1 1])
            axis equal
            
            figure();
            pcolor(yi,xi,zGridPert_l); shading('flat');
            set(gca,'fontsize',12,'fontname','times','clim',[-0.09 0.09]);
            hcb=colorbar;
            title(hcb,'[m]')
            hold on
            text(2,2,'1','fontname','times')
            xlabel('x (m)','fontsize',12,'fontname','times');
            ylabel('y (m)','fontsize',12,'fontname','times');
            runn = [day '-' month '-' year ' ' hour(1:2) ':' hour(3:4)];
            title(runn,'fontname','times');
            axis([-2.5 2.5 -2.5 2.5])
            axis equal
        end
        
        %% store data
        data05.x = yi;
        data05.y = xi;
        data05.z = zGridPert_s;
        data05.e = eGrids;
        data05.lx = lxs;
        data05.ly = lys;
        data05.z_mean = Zmean_s;
        data05.header = header;
        data05.SONARdata = SONAR3D;
        data05.flag = 1;
        dat = data05;
        save([outputPath,'Matfiles\',campaign,'_',frame,'_Depth',num2str(SONAR3D(1).depth),'_Filter05_',year,month,day,hour,'.mat'],'dat');
        
        data10.x = yi;
        data10.y = xi;
        data10.z = zGridPert_l;
        data10.e = eGridl;
        data10.lx = lxl;
        data10.ly = lyl;
        data10.z_mean = Zmean_l;
        data10.header = header;
        data10.SONARdata = SONAR3D;
        data10.flag = 1;
        dat = data10;
        save([outputPath,'Matfiles\',campaign,'_',frame,'_Depth',num2str(SONAR3D(1).depth),'_Filter10_',year,month,day,hour,'.mat'],'dat');
        %
        print([outputPath,'Figures\',campaign,'_',frame,'_Depth',num2str(SONAR3D(1).depth),'_Filter05_',year,month,day,hour,'.png'],'-dpng','-r600','-f1');
        print([outputPath,'Figures\',campaign,'_',frame,'_Depth',num2str(SONAR3D(1).depth),'_Filter10_',year,month,day,hour,'.png'],'-dpng','-r600','-f2');
        close all
    end
end
