function Y = poly_detrend(t,y,p);

% function Y = poly_detrend(t,y,p);
% poly_detrend verwijdert een trend met een polynome vorm van orde p
% standaard is een 2de orde polynoom
%
% INPUT
%   t is x-as met tijdseenheden
%   y is tijdserie waarvan de trend moet worden verwijderd
%   p is orde van de polynoom (default = 2)
% OUTPUT
%   Y is de gedetrende tijdserie

if nargin < 3, p = 2; end;

t = t(:);
y = y(:);

s = polyfit(t,y,p);
fit = polyval(s,t);
Y = y - fit;
