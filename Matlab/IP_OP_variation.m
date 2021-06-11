% Run a few bad files with multiple IP/OP choices
addpath(genpath('Functions'));
Foldername = '../EBSD_Data/bad_ang';
Fnames = dir([Foldername '\AF' '*.ang']);
OR = [3.09,8.10,8.48];      % AF96 Sample
IP_Weight = [2,3,4,5,6];    %4
IP_Scale = [4,5,6,7,8];     %6
OP_Weight = [1,2];          %2
OP_Scale = [0.175];% [0,0.175,0.25] %0.175
ang_tol = [1,2,5];% 2

for ii = 1:3% length(Fnames)
    fname = [Foldername,'/',Fnames(ii).name];
    
    
    Initial_Recon_Dataset.Material = 'Steel';
    Initial_Recon_Dataset.Phase.Name{1} = 'Martensite';
    Initial_Recon_Dataset.Phase.Name{2} = 'Austenite';
    Initial_Recon_Dataset.Celldims{1} = [2.87, 2.87, 2.87];
    Initial_Recon_Dataset.Celldims{2} = [3.65, 3.65, 3.65];
    Initial_Recon_Dataset.rec_space = 'Mixed'; 
    [Initial_Recon_Dataset] = import_EBSD_Alex(fname,Initial_Recon_Dataset);
    Initial_Recon_Dataset.Ebsd.y = flipud(Initial_Recon_Dataset.Ebsd.y);
    Initial_Recon_Dataset.OR  = [3.09,8.10,8.48];
    Initial_Recon_Dataset.noise.halfwidth = 1.7*degree;
    % [myEBSD] = AutoOR_estimation(myEBSD,0,0,2000,0);
    [Initial_Recon_Dataset] = calcMODF(Initial_Recon_Dataset);
    [Initial_Recon_Dataset] = DataQuadrants(Initial_Recon_Dataset);
    [Initial_Recon_Dataset] = initialize_recon(Initial_Recon_Dataset);

    for a = 1:length(IP_Weight)
        for b = 1:length(IP_Scale)
            IP = [IP_Weight(a) IP_Scale(b)];
            for c =1:length(OP_Weight)
                d =1;
                %                for d = 1:length(OP_Scale)
                OP = [OP_Weight(c) OP_Scale(d)];
                for e = 1:length(ang_tol)
                    try
                        Writeout_name = [Fnames(ii).name(1:end-4), '_IP', int2str(IP(1)),'-',int2str(IP(2)), '_OP', int2str(OP(1)),'-',int2str(OP(2)), '_ang',int2str(ang_tol(e))]
                        RunRecon_IP_OP(Initial_Recon_Dataset,fname,IP,OP,ang_tol(e));
                        
                        % Make some plots
                        % All plots are 4x6 inches on standard paper size, and should be
                        % represented by the same number of pixels.
                        
                        % Plot 1: EBSD Map
                        plotEBSD(myEBSD.Ebsd,myEBSD)
                        set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
                        set(gcf, 'PaperPositionMode', 'auto')
                        Plt1name = [Foldername,'/',Writeout_name,'_EBSD.png'];
                        saveas(gcf,Plt1name)
                        
                        % Plot 2: Prior Austenite Grains
                        mbar = 1;
                        plotEBSD(myEBSD.AusRecon_Ebsd,myEBSD,mbar)
                        hold on
                        plot(myEBSD.AusGrains.grains.boundary)
                        set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
                        set(gcf, 'PaperPositionMode', 'auto')
                        Plt1name = [Foldername,'/',Writeout_name,'_Austenite.png'];
                        saveas(gcf,Plt1name)
                        
                        % Plot 3: Sub Blocks
                        PltSubBlocks(myEBSD,mbar)
                        hold on
                        plot(myEBSD.AusGrains.grains.boundary,'FaceColor','white')
                        set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
                        set(gcf, 'PaperPositionMode', 'auto')
                        Plt1name = [Foldername,'/',Writeout_name,'_Variants.png'];
                        saveas(gcf,Plt1name)
                        
                        % Plot 4: Likelihood plots
                        figure; plot(myEBSD.AusRecon_Ebsd,myEBSD.AusRecon_Likelihood,'MicronBar','off')
                        hold on
                        plot(myEBSD.AusGrains.grains.boundary,'lineColor','red','linewidth',2)
                        set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
                        set(gcf, 'PaperPositionMode', 'auto')
                        Plt1name = [Foldername,'/',Writeout_name,'_Likelihood.png'];
                        saveas(gcf,Plt1name)
                        
                        mat_name = [Foldername,'/',Writeout_name,'_Recon.mat'];
                        save(mat_name,'myEBSD')
                        
                    catch
                    end
                end
            end
        end
    end
end
%end