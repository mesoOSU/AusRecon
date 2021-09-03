function par_run(input,Initial_with_MODF,ang_filename)
%cd ../mtex-5.1.1
%startup_mtex
%cd ../IP_OP_Variation
%ls
name = input.Writeout_name;
% disp('aaaaaaa')
% disp(input)
% disp('aaaaaaa')
% disp(Initial_with_MODF)
% disp('aaaaaaa')
% disp(ang_filename)
% disp('aaaaaaa')
save_name = [fileparts(ang_filename) '\R' name]
% disp('aaaaaaa')
% getappdata(0,'mtex')
if  ~isempty(dir([save_name '.mat']))
    disp([save_name '.mat already exists, skipping...'])
else
    myEBSD = Initial_with_MODF;
    iters = [];
    [myEBSD,~,~] = Call_Recon(myEBSD,input.IP,input.OP,iters);
    Original = myEBSD.Ebsd;
    Reconstruction = myEBSD.Recon.Ebsd;
    Likelyhood = myEBSD.Recon.Likelihood;
    save(save_name,'Original','Reconstruction','Likelyhood');
end
end


