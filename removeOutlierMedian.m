function n=removeOutlierMedian(data,val,dummy)

%REMOVEOUTLIERMedian removes outliers on basis of deviaton from the median
% Input
%   data is the data to be filtered
%   val is deviation from median when datapoint is considered an outlier
%   dummy is the replacement for outliers
% Output
%   n is the filtered data

n.data=data;
n.dummy=dummy;
n.val=val;

tt=data;
tt=tt(:);
n.med=median(tt(tt~=dummy));
n.data(abs(data-n.med)>val)=dummy;