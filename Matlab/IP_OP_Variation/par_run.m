function par_run(input,Initial_with_MODF,ang_filename)
myEBSD = Initial_with_MODF;
iters = [];
[myEBSD,~,~] = Call_Recon(myEBSD,input.IP,input.OP,iters);
Original = myEBSD.Ebsd;
Reconstruction = myEBSD.Recon.Ebsd;
Likelyhood = myEBSD.Recon.Likelihood;
name = input.Writeout_name
save([fileparts(ang_filename) '/R' name],'Original','Reconstruction','Likelyhood');
end


