function f=fldnrmCDF(x,mu,sigma)
%FLDNRMPDF  computes the probability density function of the folded normal
% distribution
%
% REFERENCES
% [1] FC Leone, LS Nelson, RB Nottingham, Technometrics 3 (1961) 543.
% [2] RC Elandt, Technometrics 3 (1961) 551.
%
% 2011-07-07 | Eric Payton, Ruhr-Universitaet Bochum (payton.28[at]osu.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=sigma*sqrt(2);
F=0.5*(erf((x+mu)/a)+erf((x-mu)/a));