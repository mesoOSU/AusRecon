function [log_like] = martensite_posterior_log_likelihood(martensite,ksi,halfwidth,austenite,prior_pars)

%calculate data likelihood
data_log_likelihood=sum(log(martensite_likelihood(martensite,ksi,austenite,halfwidth,prior_pars.CS_A,prior_pars.CS_M)));

ksi_log_prior_probability=sum(log(ksi_prior(ksi,prior_pars.ksi_prior_mu,prior_pars.ksi_prior_sigma)));
halfwidth_log_prior_probability=log(halfwidth_prior(halfwidth,prior_pars.halfwidth_prior_mu,prior_pars.halfwidth_prior_sigma));
%austenite_log_prior_probability=log(austenite_prior(austenite,prior_pars.austenite_prior_odf));

log_like=data_log_likelihood+ksi_log_prior_probability+halfwidth_log_prior_probability;%+austenite_log_prior_probability;

% if log_like<-500
%     log_like=-500;
% else if isnan(log_like)
%         log_like=-500;
%     end
% end

end
        
    


