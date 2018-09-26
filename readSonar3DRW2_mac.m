function [header, data, range] = readSonar3DRW2(fullSonar3DFilename)

% [header, data, range] = readSonar3DRW2(fullSonar3DFilename)
%
% RW2 files (Marine Electronics) contain an ASCII header and a binary data
% block. readHeader reads and returns the header information as a structure
% called "header". The actual data is contain in the matrices "data" and
% "range".
%
% INPUT
%   fullSonar3DFilename: path+fillname of the RW2 file
% OUTPUT
%   header: structure with settings of the sonar, necessary for later
%           processing
%   data: backscatter data for all swaths
%   range: range [m]
%   
% v1, Gerben Ruessink, 19/April/2013
% v2, Gerben Ruessink, 7/November/2013: updated to read entire file, rather
%                      than just the header
% v3, Joost Brinkkemper, 6/November/2014: made more general, for BARDEX and
%                      Zandmotor data
%
% TO DO: update help once I understand again how the binary file is built
% up.

% open file; return error if file does not exist
fidSonar3D = fopen(fullSonar3DFilename,'r');
if fidSonar3D == -1,
    error('File not found\n');
end;

% line 1: 3D_Profiler_RW2
header.fileType = fgetl(fidSonar3D);

% line 2: date time
dateTime = fgetl(fidSonar3D);
if ~isempty(strfind(dateTime,'-'))              % Added for Bardex data
    dateTime = strrep(dateTime, '-', '/');
end

header.whenStarted = datevec(dateTime,'dd/mm/yy HH:MM:SS');


% line 3: blank
header.comment = fgetl(fidSonar3D);

% lines 4-26
for line = 4:26
    
    % read line
    information = fgetl(fidSonar3D);
    
    % find location of =
    idIS = strfind(information,'=');
    
    % split information into string and number
    informationString = information(1:idIS-1);
    informationString = strrep(informationString,' ','');           % remove any blanks
    informationString = matlab.lang.makeValidName(informationString,'ReplacementStyle','delete');
    informationString = strrep(informationString,'in09steps','');
    informationString = strrep(informationString,'insec','');
    informationNumber = str2double(information(idIS+1:end));        % number
    header.(informationString) = informationNumber;                 % add2header
    
end;

% line 27: blank
fgetl(fidSonar3D);

% read the binary data
data = NaN(header.NoofSamples,floor(header.Arc/header.SwathStep),floor(200/header.RotateStep));
for rot = 1:floor(200/header.RotateStep)
    for swath = 1:floor(header.Arc/header.SwathStep)
        data(:,swath,rot) = fread(fidSonar3D,header.NoofSamples,'uint8');
    end;
end;

% and close the file
fclose(fidSonar3D);

% compute the range
range = (header.Delay:header.SampleInterval:header.Delay+(header.NoofSamples-1)*header.SampleInterval)...
        *1e-6*header.CaptureVoS/2;

% done
end
