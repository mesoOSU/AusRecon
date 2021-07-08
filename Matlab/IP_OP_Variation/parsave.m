function par_run(fname, input)

inputs(i).myEBSD = Initial_with_MODF;
[inputs(i).myEBSD,inputs(i).Parent,inputs(i).Twin] = Call_Recon(...
    inputs(i).myEBSD,inputs(i).IP,inputs(i).OP,inputs(i).iters);
inputs(i).Original = inputs(i).myEBSD.Ebsd;
inputs(i).Reconstruction = inputs(i).myEBSD.Recon.Ebsd;
inputs(i).Likelyhood = inputs(i).myEBSD.Recon.Likelihood;
parsave([fileparts(ang_filename) '/R' Writeout_name],inputs(1));
inputs(i).myEBSD = [];
inputs(i).Original = [];
inputs(i).Reconstruction = [];
inputs(i).Likelyhood = [];

Original = input.Original;
Reconstruction = input.Reconstruction;
Likelyhood = input.Likelyhood;
save(fname,'Original','Reconstruction','Likelyhood');
end