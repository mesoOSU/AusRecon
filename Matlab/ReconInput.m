%% Run Recontruction Algorithm
% Script that runs the reconstruction algorithm. For the input to the
% function 'RunRecon', the first is the filename; second is the OR if it 
% has already been measured (automated doesn't work as of now); 
% OPTIONAL INPUTS: 3rd/4th is the set of interface parameters 
% (regularization); 3rd/4th is a structure that will bypass redundant 
% calculations, (such as the MODF), for optimization of the interface 
% parameters. The order of these doesn't matter.

% Jump into mtex folder (alter as needed based on version) and start it up
% NOTE: THIS IS NOT IN THE REPO BY DEFAULT. If you need MTEX, download the 
% it from the mtex website: https://mtex-toolbox.github.io/download
cd ..
cd ..
cd mtex-5.1.1
startup_mtex
cd ..
cd AusRecon
cd Matlab
addpath(genpath('Functions'));
%% Access Reconstruction Algorithm
% Folder/filename
% fname = 'EBSD_Data/Timkin1.ctf';

%fname = 'C:\Users\paytonej\Desktop\AusRecon-main\AusRecon-main\EBSD_Data\AF96_example.ang'

fnames = {'C:\Users\paytonej\Desktop\EBSD Data-20210707T175853Z-001\EBSD Data\Ultrafort 6355\Fine (860 deg C, 4 cycles)\S500_A3_RA2-Scan2\S500_A3_RA2-Scan2.ang',...
'C:\Users\paytonej\Desktop\EBSD Data-20210707T175853Z-001\EBSD Data\Wrought AF9628\Fine (850 deg C, 4 Cycles)\2B-RA5\2B-RA5.ang',...
'C:\Users\paytonej\Desktop\EBSD Data-20210707T175853Z-001\EBSD Data\P675 Steel\Fine Microstructure\P675_P-RA1-UL\P675_P-RA1-UL.ang',...
'..\EBSD_Data\AF96_example.ang'};

for fname_iter = 1:length(fnames)
    close all;
    %clear('myEBSD')
    fname = fnames{fname_iter};
% Possibilities after a single run
% OR = myEBSD.OR; % You can type this myEBSD.OR into the workspace and then
% OR = [3.09,8.10,8.48];      % AF96 Sample [0.915447390857388,9.245004817251438,9.258634331312445]
% use OR = [...] (whatever the values are) 
% intparams = [...] % Whatever you want them to be
% myEBSD is the structure that returns everything. If you add myEBSD to the
% function input, then you bypass modf calculation, etc. and just get to
% the reconstruction.

% When you are optimizing the results, run:

RunRecon(fname, [3.09, 8.10, 8.48])
% RunRecon(fname)
% You can still add OR, and/or intparams; myEBSD to the RunRecon1 function 
% input as well! This will perform reconstruction on a truncated portion of
% the dataset to save time.

% Run for the entire dataset!
% RunRecon(fname,OR,intparams)

%% Plot Initial Transformation Microstructure
plotEBSD(myEBSD.Ebsd,myEBSD)

%% Plot Austenite Reconstruction With Overlying Grain Boundaries

% Plot microstructure. Input is the Ebsd dataset desired to be plotted, the
% and the myEBSD structure (for phase information). If desired, add a third
% variable which will remove the micronbar from the plot. Final line will
% outline the determined PAG boundaries in black.
% figure; plot(AusRecon_Ebsd,AusRecon_Ebsd.orientations,'MicronBar','off')
mbar = 1;
plotEBSD(myEBSD.Recon.Ebsd,myEBSD,mbar)
hold on
plot(myEBSD.AusGrains.grains.boundary)

%% Plot Austenite Likelihood Plot

% Plot likelihood plot of chosen austenite orientations transforming the
% observed martensitic microstructure
figure; plot(myEBSD.Recon.Ebsd,myEBSD.Recon.Likelihood,'MicronBar','off')

%% Plot Grain Structure With Labels

% Plot grain structure or merged twins with the corresponding numbering of
% each grain (if set to 1; otherwise, it will just plot the grain boundary
% structure without labeling
 plotGrains(myEBSD.AusGrains,1);

%% Save results

tmp1 = split(fname, filesep);
tmp2 = join(tmp1(1:end-1), filesep);
tmp3 = split(tmp1(end), '.');
outputfilename = strcat(tmp2{1}, filesep, tmp3{1}, '.mat');
save(outputfilename)

end % loop over files