function p=ksi_prior(ksi_sample,mu,sigma)

%prior distribution on the orientation relationship.

% inputs
%ksi_sample=[ksi1,ksi2,ksi3]
%mu and sigma and 1x3 column vectors containing the means and standard
%deviations of the individual ksi terms

% outputs
% p=[p1,p2,p3] collumn vector of the prior probabilities
%
if numel(find(ksi_sample<0))>0
    p=[1,1,1]*1e-100;
else if ksi_sample(1)>ksi_sample(2)
        p=[1,1,1]*1e-100;
    else
        
        
        p1=fldnrmPDF(ksi_sample(1),mu(1),sigma(1));
        p2=fldnrmPDF(ksi_sample(2),mu(2),sigma(2));
        p3=fldnrmPDF(ksi_sample(3),mu(3),sigma(3));
        
        p=[p1,p2,p3];
    end
end