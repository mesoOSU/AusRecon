function [p] = posterior_pdf_fminunc(samples,prior_pars,martensite)

ksi1=abs(samples(1));
ksi2=abs(samples(2));
ksi3=abs(samples(3));
halfwidth=abs(samples(4))*degree;

[T2R flag]=calc_T2R([ksi1,ksi2,ksi2],prior_pars.CS_A,prior_pars.CS_M);

num_martensite_orientations=length(martensite);

austenite=symmetrise(martensite)*T2R;

[austenite_proposal] = global_pole_figure_estimation(austenite,austenite.CS,austenite.SS,1);

for kk=1:length(austenite_proposal)
    temp(kk)=martensite_posterior_log_likelihood(martensite,[ksi1,ksi2,ksi3],halfwidth,austenite_proposal(kk),prior_pars);
end
id=find(temp==max(temp));
p=temp(id);%ma

%p=martensite_posterior_log_likelihood(martensite,[ksi1,ksi2,ksi3],halfwidth,austenite,prior_pars);
end

