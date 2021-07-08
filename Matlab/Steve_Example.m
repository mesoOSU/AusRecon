% Run Recontruction Algorithm

addpath(genpath('Functions'));
Foldername = '../EBSD_Data/bad_ang';
% This searches through the folder names and collects them into a single
% array
fname = [Foldername '\' 'AF_005.ang'];
% The AF96 orientation relationship. The algorithm can now automatically
% determine it if need be.
OR = [3.09,8.10,8.48];      % AF96 Sample

Pltname = fname(1:end-4)
RunRecon(fname,OR,[4,6])

% Make some plots
% All plots are 4x6 inches on standard paper size, and should be
% represented by the same number of pixels.

% Plot 1: EBSD Map
plotEBSD(myEBSD.Ebsd,myEBSD)
set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
set(gcf, 'PaperPositionMode', 'auto')
Plt1name = [Foldername,'/',Pltname,'_EBSD.png'];
saveas(gcf,Plt1name)

% Plot 2: Prior Austenite Grains
mbar = 1;
plotEBSD(myEBSD.AusRecon_Ebsd,myEBSD,mbar)
hold on
plot(myEBSD.AusGrains.grains.boundary)
set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
set(gcf, 'PaperPositionMode', 'auto')
Plt1name = [Foldername,'/',Pltname,'_Austenite.png'];
saveas(gcf,Plt1name)

% Plot 3: Sub Blocks
PltSubBlocks(myEBSD,mbar)
hold on
plot(myEBSD.AusGrains.grains.boundary,'FaceColor','white')
set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
set(gcf, 'PaperPositionMode', 'auto')
Plt1name = [Foldername,'/',Pltname,'_Variants.png'];
saveas(gcf,Plt1name)

% Plot 4: Likelihood plots
figure; plot(myEBSD.AusRecon_Ebsd,myEBSD.AusRecon_Likelihood,'MicronBar','off')
hold on
plot(myEBSD.AusGrains.grains.boundary,'lineColor','red','linewidth',2)
set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
set(gcf, 'PaperPositionMode', 'auto')
Plt1name = [Foldername,'/',Pltname,'_Likelihood.png'];
saveas(gcf,Plt1name)
