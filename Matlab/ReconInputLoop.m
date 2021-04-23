%% Run Recontruction Algorithm
% Script that runs the reconstruction algorithm. For the input to the
% function 'RunRecon', the first is the filename; second is the OR if it 
% has already been measured (automated doesn't work as of now); 
% OPTIONAL INPUTS: 3rd/4th is the set of interface parameters 
% (regularization); 3d/4th is a structure that will bypass redundant 
% calculations, (such as the MODF), for optimization of the interface 
% parameters. The order of these doesn't matter.

% Jump into mtex folder (alter as needed based on version) and start it up
% NOTE: THIS IS NOT IN THE REPO BY DEFAULT. If you need MTEX, download the 
% it from the mtex website: https://mtex-toolbox.github.io/download
%cd mtex-5.1.1
%startup_mtex
%cd ..
addpath(genpath('Functions'));
%% Access Reconstruction Algorithm
% Whatever the folder name is. For the example below, a DISTRO A copy of the
% data can be found at either of the below locations:
%AF Research Lab link: 
%https://drive.google.com/drive/folders/1cw21aKhuaXX_yEq_0JA0GJTm7L0Oi91o?usp=sharing
%ELSZ Storage File location:
%sftp://sshfs.rdte.afrl.dren.mil/   location: /project/RXCM/DISTRO_A_Datasets_AusRecon
Foldername = 'EBSD_Data/angFiles';
% This searches through the folder names and collects them into a single
% array
Fnames = dir(Foldername);
% The AF96 orientation relationship. The algorithm can now automatically
% determine it if need be.
OR = [3.09,8.10,8.48];      % AF96 Sample

% CHECK YOURSELF! On mine, the first 2 structures are '.'1 and '..', so the
% file names begin at 3. This could be different for you!
for ii = 3:length(Fnames)
%    try
    % Path and file name
    fname = [Foldername,'/',Fnames(ii).name]
    % Remove '.ang' from filename to save plots as the filename.png
    Pltname = Fnames(ii).name;
    Pltname = strrep(Pltname,'.ang','');

    % Run the algorithm. Remove OR to automatically determine it.
    RunRecon(fname,OR)
    
    % All plots are 4x6 inches on standard paper size, and should be
    % represented by the same number of pixels.
    %% Plot Grain Boundaries
    figure; plot(AusGrains.grains.boundary,'MicronBar','off')
    set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
    set(gcf, 'PaperPositionMode', 'auto')
    Plt1name = [Foldername,'/',Pltname,'_','Grain_Boundaries','.png'];
    saveas(gcf,Plt1name)

    %% Plot Sub-Blocks
    PltSubBlocks(myEBSD,1)
    set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
    set(gcf, 'PaperPositionMode', 'auto')
    Plt2name = [Foldername,'/',Pltname,'_','Variants','.png'];
    saveas(gcf,Plt2name)
%    catch
%    end
    
end