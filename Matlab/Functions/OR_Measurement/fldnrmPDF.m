function f=fldnrmpdf(x,mu,sigma)
%FLDNRMPDF  computes the probability density function of the folded normal
% distribution
%
% REFERENCES
% [1] FC Leone, LS Nelson, RB Nottingham, Technometrics 3 (1961) 543.
% [2] RC Elandt, Technometrics 3 (1961) 551.
%
% 2011-07-07 | Eric Payton, Ruhr-Universitaet Bochum (payton.28[at]osu.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a=2*sigma^2;
f=(1/(sqrt(2*pi)*sigma))*(exp(-(x-mu).^2./a)+exp(-(x+mu).^2./a));
f(x<0)=1e-10;