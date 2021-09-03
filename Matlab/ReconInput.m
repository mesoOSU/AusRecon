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
% cd mtex-5.1.1
% startup_mtex
% cd ..
addpath(genpath('Functions'));
%% Access Reconstruction Algorithm
% Folder/filename
% fname = 'EBSD_Data/Timkin1.ctf';

fname = '../EBSD_Data/AF96_example.ang';

% Possibilities after a single run
% OR = myEBSD.OR; % You can type this myEBSD.OR into the workspace and then
OR = [3.09,8.10,8.48];      % AF96 Sample [0.915447390857388,9.245004817251438,9.258634331312445]
% use OR = [...] (whatever the values are) 
% intparams = [...] % Whatever you want them to be
% myEBSD is the structure that returns everything. If you add myEBSD to the
% function input, then you bypass modf calculation, etc. and just get to
% the reconstruction.

% When you are optimizing the results, run:

RunRecon(fname)
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

%% Plot Segmented Packets Over Entire Microstructure 

% Plot the packets for generated austenite grains for the entire
% microstructure. Again, mbar gives the choice of removing the micron bar 
% from the image and the coloring scheme uses the hsv color map.
mbar = 1;
PltAllPackets(myEBSD,mbar);
hold on
plot(myEBSD.AusGrains.grains.boundary)

%% Plot Segmented Blocks Over Entire Microstructure
% Same concept for plotting blocks except the block colors adhere to the 
% custom colormap and both grain boundaries (white) and packet boundaries
% (black) are overlayed over the microstructure

PltAllBlocks(myEBSD,mbar);
hold on
plot(myEBSD.AusGrains.grains.boundary,'FaceColor','white')

%% Plot Sub-Blocks

% Plots a map of sub-blocks from 1 to 24 (twin variants are reassigned to
% the parent numbering scheme). Block boudnaries are overlayed in 'black'
% and the PAG boundaries in 'white.'

PltSubBlocks(myEBSD,mbar)
hold on
plot(myEBSD.AusGrains.grains.boundary,'FaceColor','white')

% You don't need to run anything below HERE!!!

%% Plot Single Austenite Grain And Segmented Packets

% Plot single austenite grain packets for individual analysis of a grain.
% The same packet coloring scheme is used from above. Output is the 
% martensite variants from within the chosen PAG; packet boundaries; and 
% the corresponding weights.
grnId = 7;
PltSngGrnPackets(myEBSD,myEBSD.AusGrnPackets,grnId)

%% Plot Packet Orientation Pole Figures (Both Theoretical and Segmented)

% Plots the pole figures that relate the segmented martensite variants from 
% the selected grain to what should be theoretically expected. Colormap is
% consistent with Packets.

% Directions (can choose up to three)
Dir = [0,0,1];
PltPacketPDFs(myEBSD,myEBSD.AusGrnPackets,grnId,Dir)

%% Plot Single Austenite Grain And Segmented Blocks

% Plots the block boundaries over a chosen PAG  using the same custom
% colormap used to plot block boundaries over the entire microstructure.
% Same output as the packet case except substitute blocks for packets.
grnId = 7;
PltSngGrnBlocks(myEBSD,myEBSD.AusGrnPackets,grnId)

%% Plot Block Orientation Pole Figures (Theoretical and Segmented)
% Directions (can choose up to three)
grnId = 7;

Dir = [0,0,1];
PltBlockPDFs(myEBSD,myEBSD.AusGrnPackets,grnId,Dir)

%% Plot Paired Bain Variants (Even And Odd)

% Plots a map of the closest matching sub-block variants across the entire
% microstructure (separates 'even' and 'odd' labeled variants, as these
% correspond to sub-blocks). Block boundaries are overlayed in 'black,' and
% PAG boundaries in 'white.'
PltPairedVars(myEBSD,mbar)
hold on
plot(myEBSD.AusGrains.grains.boundary,'FaceColor','white')