%MTEX setup stuff
clear all
close all
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
addpath(genpath('Functions'));
OR = [3.09,8.10,8.48];                % AF96 Sample

% Choose some EBSD scans to test out
Foldername = '../EBSD_Data/bad_ang';
Fnames = dir([Foldername '\AF' '*.ang']);


%Choose IP and OP weights to iterate over
%IP_Weight = [3.0,3.5,4.0,4.5,5.0];        %4
%IP_Scale =  [5.0,5.5,6.0,6.5,7.0];        %6
%OP_Weight = [1.0,1.5,2.0,2.5,3.0];        %2
%OP_Scale  = [0.1,0.15,0.175,0.2,0.225];       %0.175

IP_Weight = [30,40,50];        %4
IP_Scale =  [50,60,70];        %6
OP_Weight = [10,20,30];        %2
OP_Scale  = [100,175,200];     %0.175


for ii = 1:length(Fnames)
    fname = [Foldername,'/',Fnames(ii).name];
    Initial_Recon_Dataset.Material = 'Steel';
    Initial_Recon_Dataset.Phase.Name{1} = 'Martensite';
    Initial_Recon_Dataset.Phase.Name{2} = 'Austenite';
    Initial_Recon_Dataset.Celldims{1} = [2.87, 2.87, 2.87];
    Initial_Recon_Dataset.Celldims{2} = [3.65, 3.65, 3.65];
    Initial_Recon_Dataset.rec_space = 'Mixed';
    [Initial_Recon_Dataset] = import_EBSD_Alex(fname,Initial_Recon_Dataset);
    % Initial_Recon_Dataset.Ebsd.y = flipud(Initial_Recon_Dataset.Ebsd.y);
    Initial_Recon_Dataset.OR  = [3.09,8.10,8.48];
    Initial_Recon_Dataset.noise.halfwidth = 1.7*degree;
    [Initial_Recon_Dataset] = calcMODF(Initial_Recon_Dataset);
    [Initial_Recon_Dataset] = DataQuadrants(Initial_Recon_Dataset);
    [Initial_Recon_Dataset] = initialize_recon(Initial_Recon_Dataset);
    CS = Initial_Recon_Dataset.origEbsd('Austenite').CS;
    
    % define S3 and S9. (Used later for plotting grains)
    %S3 is defined using u1v1 method from Steve
    u1 = Miller( 1, 2, 1,CS);
    v1 = Miller(-1,-2,-1,CS);
    u2 = Miller(-1, 0, 1,CS);
    v2 = Miller( 1, 0,-1,CS);
    S3_twin_miso = orientation('map',u1,v1,u2,v2);
    % S9 was done from 1984 Hans Grimmer Paper
    S9_angle = 38.94*degree;
    S9_axis = vector3d(1,1,0);
    S9_twin_miso = orientation('axis',S9_axis,'angle',S9_angle,CS);
    S9_twin_miso.SS = S9_twin_miso.CS;
    
    a = 2;
    %    for a = 1:length(IP_Weight)
    for b = 1:length(IP_Scale)
        IP = [IP_Weight(a)/10 IP_Scale(b)/10];
        for c =1:length(OP_Weight)
            for d = 1:length(OP_Scale)
                OP = [OP_Weight(c)/10 OP_Scale(d)/1000];
                try
                    IP_name = ['_IP',int2str(IP_Weight(a)),'-',int2str(IP_Scale(b))];
                    OP_name = ['_IP',int2str(OP_Weight(c)),'-',int2str(OP_Scale(d))];
                    Writeout_name = [Fnames(ii).name(1:end-4), IP_name,OP_name];
                    clear IP_name OP_name
                    iters = [];
                    
                    % Run the Reconstruction
                    myEBSD = Initial_Recon_Dataset;
                    [myEBSD,Parent,Twin] = Call_Recon(myEBSD,IP,OP,iters);
                    
                    % Extract the outputs and save them
                    Original = myEBSD.Ebsd;
                    Reconstruction = myEBSD.Recon.Ebsd;
                    Likelyhood = myEBSD.Recon.Likelihood;
                    save_na = split(fname,'/');
                    save_nam = string(join(save_na(1:end-1),'/'));
                    save_name = strcat(save_nam, '//', Writeout_name);
                    save(save_name,'Original','Reconstruction','Likelyhood');
                    clear save_na save_nam save_name my_EBSD
                    
                    % The reconstruction is done, the rest of this is
                    % Austin's methods for plotting what we want
                    
                    Mart = Original('Martensite');
                    Aus  = Reconstruction('Austenite');
                    
                    % segment
                    [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(Aus.fill,'angle',1*degree);
                    gB = grains.boundary;
                    gB_Aus = gB('Austenite','Austenite');
                    
                    % Find s3 and s9 twins
                    S3_isTwinning = angle(gB_Aus.misorientation,S3_twin_miso) < 0.5*degree;
                    S9_isTwinning = angle(gB_Aus.misorientation,S9_twin_miso) < 0.5*degree;
                    S3_TwinBoundaries = gB_Aus(S3_isTwinning);
                    S9_TwinBoundaries = gB_Aus(S9_isTwinning);
                    All_Twin_Boundaries = [S9_TwinBoundaries S3_TwinBoundaries];
                    [mergedGrains,parentId] = merge(grains,All_Twin_Boundaries);
                    
                    close all
                    figure(3)
                    plot(Mart.fill,Mart.fill.orientations,'MicronBar','off')
                    hold on
                    plot(All_Twin_Boundaries,'linecolor','w','linewidth',4.5)
                    plot(S3_TwinBoundaries,'linecolor','r','linewidth',2.5,'displayName','S3')
                    % plot(S9_TwinBoundaries,'linecolor','g','linewidth',2.5,'linestyle','-','displayName','S9')
                    plot(mergedGrains.boundary,'linecolor','k','linewidth',3,'displayName','Aus')
                    set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
                    set(gcf, 'PaperPositionMode', 'auto')
                    Plt1name = [Foldername,'/',Writeout_name,'_Likelihood.png'];
                    saveas(gcf,Plt1name)
                    hold off
                    
                    figure(4)
                    plot(Mart,Likelyhood,'MicronBar','off')
                    hold on
                    plot(All_Twin_Boundaries,'linecolor','w','linewidth',4.5)
                    plot(S3_TwinBoundaries,'linecolor','r','linewidth',2.5,'displayName','S3')
                    % plot(S9_TwinBoundaries,'linecolor','g','linewidth',2.5,'linestyle','-','displayName','S9')
                    plot(mergedGrains.boundary,'linecolor','k','linewidth',3,'displayName','Aus')
                    set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
                    set(gcf, 'PaperPositionMode', 'auto')
                    Plt1name = [Foldername,'/',Writeout_name,'_Parent_Grains.png'];
                    saveas(gcf,Plt1name)
                    hold off
                    
                catch
                end
            end
        end
    end
end