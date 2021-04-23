function orientations = sample_halfwidth(mean,eps,num_samples,CS,SS)

odf_pertube = unimodalODF(mean,'halfwidth',eps);
orientations = calcOrientations(odf_pertube,num_samples);
%q = quaternion(orientations);