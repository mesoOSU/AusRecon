function q = sample_VMF(mean,kappa,num_samples)
% backdoor way to sample from a Von Mises Fisher distribution using MTEX

if kappa>=500
    kappa=499;
end
psi = vonMisesFisherKernel(kappa);
odf_pertube = unimodalODF(mean,psi);
orientations = calcOrientations(odf_pertube,num_samples);
q = quaternion(orientations);