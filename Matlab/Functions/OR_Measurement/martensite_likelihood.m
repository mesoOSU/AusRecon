 function [likelihood] = martensite_likelihood(martensite,ksi,austenite,halfwidth,CS_R,CS_T)
%evaluates the likelihood of measuring marteniste orientations for a given
%set of parameters

aus_odf=unimodalODF(austenite,austenite.CS,austenite.SS,'halfwidth',halfwidth);

[T2R flag]=calc_T2R(ksi,CS_R,CS_T);

pot_aus=symmetrise(martensite)*T2R;

l=eval(aus_odf,pot_aus);
l(l<1e-6)=1e-6;
l=reshape(l,24,length(l)/24)';

likelihood=sum(l,2)*0.057; %factor is to convert to proper probability distribution from tiem uniform random odf


end

