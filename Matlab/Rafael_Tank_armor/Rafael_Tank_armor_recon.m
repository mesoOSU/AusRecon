%MTEX setup stuff
clear all
close all

%cd ../mtex-5.1.1
%startup_mtex
%cd ../Rafael_Tank_armor

setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
addpath(genpath('../Functions'));

% Identify the .ang files for testing
ang_filenames = dir('Data/test.ang');
% If you had a bunch of files, you could do a "for" loop here, but there is 
% just one, so inseadI set ii = 1 and move on
ii = 1;
ang_filename = [ang_filenames(ii).folder '/' ang_filenames(ii).name]
%When AusRecon gets refactored, all this setup will be automated, but for
%now need to set some variables the hard way.
clear tank_EBSD
tank_EBSD.Material = 'Steel';
tank_EBSD.Phase.Name{1} = 'Martensite';
tank_EBSD.Phase.Name{2} = 'Austenite';
tank_EBSD.Celldims{1} = [2.87, 2.87, 2.87];
tank_EBSD.Celldims{2} = [3.65, 3.65, 3.65];
tank_EBSD.rec_space = 'Mixed';
[tank_EBSD] = import_EBSD_Alex(ang_filename,tank_EBSD);
CS_A = tank_EBSD.origEbsd('Austenite').CS;
CS_M = tank_EBSD.origEbsd('Martensite').CS;
%Instead of KS or NW, this code determines the correct Aus/Mart Orientation
%relationship (OR)using a method developed by Vitoria Yardley and coded up by
%Eric payton. A bit confusing but much more accurate. check the function below
% for a link to the paper. First time you run a new sample, you should determine
%the correct OR for that sample, as it can change depending on composition and 
%heat treat.
%====================
%[tank_EBSD] = AutoOR_estimation(tank_EBSD,1,1,2000,1);
%====================
%However, AutoOR_estimation takes some time to run, so instead i ran it for
%you and saved the result below. feel free to comment out this line and
%uncomment the one above to do the OR-determination on the fly.
tank_EBSD.OR  = [0.6201,    4.9984,    4.9981];
% At this point, all the pre-Reconstruction work is done. If you would like
% to see what your sample looks like loaded into MTEX pre-construction,
% here is the commend you need:
%plot(tank_EBSD.origEbsd('Martensite'),tank_EBSD.origEbsd('Martensite').orientations)



% 225 tests per file
IP_Weight = [20,30,40,50,60];        %4
IP_Scale =  [40,50,60,70,80];        %6
OP_Weight = [10,20,30];        %2
OP_Scale  = [100,200,300];     %0.175
MODF_Kernel_noise_halfwidth = [17];

% do a bunch of nested for loops to make structs of input
for e = 1:length(MODF_Kernel_noise_halfwidth)
    MODF_HW = MODF_Kernel_noise_halfwidth(e);
    % add MODF to initial Recon Structure
    clear Initial_with_MODF
    Initial_with_MODF = tank_EBSD;
    Initial_with_MODF.noise.halfwidth = 0.1*MODF_HW*degree;
    [Initial_with_MODF] = calcMODF(Initial_with_MODF);
    [Initial_with_MODF] = DataQuadrants(Initial_with_MODF);
    [Initial_with_MODF] = initialize_recon(Initial_with_MODF);
    
    %This is the last point where we modify data. everything else is
    %parameter modification, so this is where we create an input list
    %for looping over
    clear inputs
    i=0;
    for a = 1:length(IP_Weight)
        for b = 1:length(IP_Scale)
            IP = [IP_Weight(a)/10 IP_Scale(b)/10];
            for c =1:length(OP_Weight)
                for d = 1:length(OP_Scale)
                    OP = [OP_Weight(c)/10 OP_Scale(d)/1000];
                    IP_name =['_IP',int2str(IP_Weight(a)),'-',int2str(IP_Scale(b))];
                    OP_name = ['_OP',int2str(OP_Weight(c)),'-',int2str(OP_Scale(d))];
                    HW_name = ['_HW',int2str(OP_Weight(c))];
                    Writeout_name = [ang_filenames(ii).name(1:end-4),HW_name,IP_name,OP_name];
                    clear IP_name OP_name HW_name
                    i = i+1;
                    inputs(i).IP = IP;
                    inputs(i).OP = OP;
                    inputs(i).Writeout_name = Writeout_name;
                    inputs(i).MODF_HW = MODF_HW;
                end
            end
        end
    end
    total_tests = i;
    disp(['Total of ' ,int2str(i),' inputs in this parfor loop. Starting now...'])
    save('Initial_with_MODF','Initial_with_MODF')
    for i =1:total_tests
        par_run(inputs(i),Initial_with_MODF,ang_filename);
    end
    disp('Parfor Loop finished. beginning next one...')
end
disp('======================================================')
disp('All loops finished for this file, moving to next one')
disp('======================================================')


%
%                     % The reconstruction is done, the rest of this is
%                     % Austin's methods for plotting what we want
%
%                     Mart = Original('Martensite');
%                     Aus  = Reconstruction('Austenite');
%
%                     % segment
%                     [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(Aus.fill,'angle',1*degree);
%                     gB = grains.boundary;
%                     gB_Aus = gB('Austenite','Austenite');
%
%                     % Find s3 and s9 twins
%                     S3_isTwinning = angle(gB_Aus.misorientation,S3_twin_miso) < 0.5*degree;
%                     S9_isTwinning = angle(gB_Aus.misorientation,S9_twin_miso) < 0.5*degree;
%                     S3_TwinBoundaries = gB_Aus(S3_isTwinning);
%                     S9_TwinBoundaries = gB_Aus(S9_isTwinning);
%                     All_Twin_Boundaries = [S9_TwinBoundaries S3_TwinBoundaries];
%                     [mergedGrains,parentId] = merge(grains,All_Twin_Boundaries);
%
%                     close all
%                     figure(1)
%                     plot(Mart.fill,Mart.fill.orientations,'MicronBar','off')
%                     hold on
%                     plot(All_Twin_Boundaries,'linecolor','w','linewidth',4.5)
%                     plot(S3_TwinBoundaries,'linecolor','r','linewidth',2.5,'displayName','S3')
%                     % plot(S9_TwinBoundaries,'linecolor','g','linewidth',2.5,'linestyle','-','displayName','S9')
%                     plot(mergedGrains.boundary,'linecolor','k','linewidth',3,'displayName','Aus')
%                     set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
%                     set(gcf, 'PaperPositionMode', 'auto')
%                     Plt1name = [Foldername,'/',Writeout_name,'_Likelihood.png'];
%                     saveas(gcf,Plt1name)
%                     hold off
%
%                     figure(2)
%                     plot(Mart,Likelyhood,'MicronBar','off')
%                     hold on
%                     plot(All_Twin_Boundaries,'linecolor','w','linewidth',4.5)
%                     plot(S3_TwinBoundaries,'linecolor','r','linewidth',2.5,'displayName','S3')
%                     % plot(S9_TwinBoundaries,'linecolor','g','linewidth',2.5,'linestyle','-','displayName','S9')
%                     plot(mergedGrains.boundary,'linecolor','k','linewidth',3,'displayName','Aus')
%                     set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
%                     set(gcf, 'PaperPositionMode', 'auto')
%                     Plt1name = [Foldername,'/',Writeout_name,'_Parent_Grains.png'];
%                     saveas(gcf,Plt1name)
%                     hold off
%                 catch
%                 end
%             end
%         end
%     end
% end



% Grab the saved CS, use it for defining S3 and S9
% save('Austenite_CS','CS') <save function here if you need to overwrite CS
load('Austenite_CS','CS')
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
