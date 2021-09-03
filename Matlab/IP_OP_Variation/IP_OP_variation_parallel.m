%MTEX setup stuff
clear all
close all

% cd ../mtex-5.1.1
% startup_mtex
% cd ../IP_OP_Variation

parpool(5)

setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
addpath(genpath('../Functions'));

% Identify the .ang files for testing
ang_filenames = dir('EBSD/AF_00*/*AF_00*.ang');

%Choose IP and OP weights to iterate over
%IP_Weight = [3.0,3.5,4.0,4.5,5.0];        %4
%IP_Scale =  [5.0,5.5,6.0,6.5,7.0];        %6
%OP_Weight = [1.0,1.5,2.0,2.5,3.0];        %2
%OP_Scale  = [0.1,0.15,0.175,0.2,0.225];       %0.175

% IP_Weight = [30,40,50];        %4
% IP_Scale =  [50,60,70];        %6
% OP_Weight = [10,20,30];        %2
% OP_Scale  = [100,175,200];     %0.175

% 735 tests per file
IP_Weight = [20,30,35,40,45,50,60];        %4
IP_Scale =  [40,50,55,60,65,70,80];        %6
OP_Weight = [10,15,20,25,30];        %2
OP_Scale  = [100,175,200];     %0.175
MODF_Kernel_noise_halfwidth = [17];

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


for ii = 1:length(ang_filenames)
    ang_filename = [ang_filenames(ii).folder '/' ang_filenames(ii).name];
    clear Initial_Recon_Dataset;
    Initial_Recon_Dataset.Material = 'Steel';
    Initial_Recon_Dataset.Phase.Name{1} = 'Martensite';
    Initial_Recon_Dataset.Phase.Name{2} = 'Austenite';
    Initial_Recon_Dataset.Celldims{1} = [2.87, 2.87, 2.87];
    Initial_Recon_Dataset.Celldims{2} = [3.65, 3.65, 3.65];
    Initial_Recon_Dataset.rec_space = 'Mixed';
    [Initial_Recon_Dataset] = import_EBSD_Alex(ang_filename,Initial_Recon_Dataset);
    Initial_Recon_Dataset.OR  = [3.09,8.10,8.48];
    % CS = Initial_Recon_Dataset.origEbsd('Austenite').CS;
    
    % do a bunch of nested for loops to make structs of input
    for e = 1:length(MODF_Kernel_noise_halfwidth)
        MODF_HW = MODF_Kernel_noise_halfwidth(e);
        % add MODF to initial Recon Structure
        clear Initial_with_MODF
        Initial_with_MODF = Initial_Recon_Dataset;
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
                        Writeout_name = [ang_filenames(ii).name(1:end-4),HW_name,IP_name,OP_name]
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
        parfor i =1:total_tests
            par_run(inputs(i),Initial_with_MODF,ang_filename);
        end
        disp('Parfor Loop finished. beginning next one...')
    end
    disp('======================================================')
    disp('All loops finished for this file, moving to next one')
    disp('======================================================')
end

%% Now make some default images and save num
% firsd define S3 and S9 criteria
close all
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
addpath(genpath('../Functions'));

dirlist =  dir('EBSD/AF_00*/R*.mat');
skipped = {'skipped_filenames'}

IP_weight = [];
IP_scale  = [];
OP_weight = [];
OP_scale  = [];

for i = 1:size(dirlist,1)
    close all;
    fprintf(' ------ Iter: %d ------ \n', i)
    name = [dirlist(i).folder '\' dirlist(i).name];
    load(name);
    Writeout_name = name(1:end-4);
    if length(dir([name(1:end-4) '*.png'])) <3
        try
            [Bounds,Aus_grains,Aus_parentId] = Boundary_Finder(Reconstruction,Original);
            % plot with twins on Aus
            figure
            plot(Aus_grains('Austenite'),Aus_grains('Austenite').meanOrientation,'MicronBar','off');
            hold on
            plot(Bounds.All,'linecolor','w','linewidth',4)
            try
                plot(Bounds.S3_Twin,'linecolor','#6b6b6b','linewidth',2.5)
                plot(Bounds.S9_Twin,'linecolor','#3b3b3b','linewidth',4)
            catch
            end
            plot(Bounds.All,'linecolor','k','linewidth',2)
            set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
            set(gcf, 'PaperPositionMode', 'auto')
            Pltname = [Writeout_name,'_Aus.png'];
            saveas(gcf,Pltname)
            hold off
            
            % same but on Mart
            figure
            plot(Original('Martensite'),Original('Martensite').orientations,'MicronBar','off');
            hold on
            plot(Bounds.All,'linecolor','w','linewidth',4)
            try
                plot(Bounds.S3_Twin,'linecolor','#6b6b6b','linewidth',2.5)
                plot(Bounds.S9_Twin,'linecolor','#3b3b3b','linewidth',4)
            catch
            end
            plot(Bounds.All,'linecolor','k','linewidth',2)
            set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
            set(gcf, 'PaperPositionMode', 'auto')
            Pltname = [Writeout_name,'_Mart.png'];
            saveas(gcf,Pltname)
            hold off
            
            % same but on Likelyhood
            figure
            plot(Reconstruction,Likelyhood(1:end,1),'MicronBar','off');
            hold on
            plot(Bounds.All,'linecolor','w','linewidth',4)
            try
                plot(Bounds.S3_Twin,'linecolor','#6b6b6b','linewidth',2.5)
                plot(Bounds.S9_Twin,'linecolor','#3b3b3b','linewidth',4)
            catch
            end
            plot(Bounds.All,'linecolor','k','linewidth',2)
            set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,3]);
            set(gcf, 'PaperPositionMode', 'auto')
            Pltname = [Writeout_name,'_L.png'];
            saveas(gcf,Pltname)
            hold off
        catch
            fprintf('\n\n Beans!!!! \n\n');
            skipped{size(skipped,2)+1}=name;
        end
        
    else
        fprintf('skipping...');
    end
    % grab the metadata from the filename and add to some variables to
    % make generating heatmaps easier later
    
    % the fact that i have to wrap "split" to get a single line to assign
    % multiple variables values from a cell will never NOT annoy me.
    % this for loop and preceeding empty variable nonsense should jsut be
    % one line of code
    [~,scan_ID,HW,IP,OP] =feval(@(x) x{:},split(dirlist(i).name,"_",2));
    [a,b] =feval(@(x) x{:},split(IP(3:end),"-",2));
    [c,d] =feval(@(x) x{:},split(OP(3:end-4),"-",2));
    IP_weight(i) = str2num(a);
    IP_scale(i)  = str2num(b);
    OP_weight(i) = str2num(c);
    OP_scale(i)  = str2num(d);
    
    %Now time to try out a few different ways of quantifying a "good"
    %Likelyhood plot
    
    sorted_L = sort(Likelyhood);
    L_10(i) = sorted_L(round(size(Likelyhood,1)*0.1));
    L_20(i) = sorted_L(round(size(Likelyhood,1)*0.2));
    L_80(i) = sorted_L(round(size(Likelyhood,1)*0.8));
    L_90(i) = sorted_L(round(size(Likelyhood,1)*0.9));
    L_ave(i) = mean(Likelyhood);
    %Median Filters
    L = zeros(321*321,1);
    L(Reconstruction.id)=Likelyhood;
    L = reshape(L,[321,321]);
    L_M2(i) = mean(mean(medfilt2(L)));
    L_M5(i) = mean(mean(medfilt2(L,[5,5])));
    % Looks like the above are all just variants of a mean filter, so lets
    % try something new: threshold by grain
    [grains,Reconstruction.grainId] = calcGrains(Reconstruction);
    G = grains(grains.grainSize>10).id;
    M_L = [];
    G_S = [];
    LL = Likelyhood;
    for ii = 1:length(G)
        M_L(ii) = mean(Likelyhood(Reconstruction.grainId == G(ii)));
        G_S(ii)  = sum(Reconstruction.grainId == G(ii));
        LL(Reconstruction.grainId == G(ii)) = M_L(ii);
    end
    WT50(i) = mean(LL> (mean(Likelyhood)*0.50));
    WT75(i) = mean(LL> (mean(Likelyhood)*0.75));
    WT90(i) = mean(LL> (mean(Likelyhood)*0.90));
    GT50(i) = mean(M_L>(mean(Likelyhood)*0.50));
    GT75(i) = mean(M_L>(mean(Likelyhood)*0.75));
    GT90(i) = mean(M_L>(mean(Likelyhood)*0.90));
end
close all

L_10  = L_10(L_10 >0);
L_20  = L_20(L_20 >0);
L_80  = L_80(L_80 >0);
L_90  = L_90(L_90 >0);
L_ave = L_ave(L_ave >0);
L_M2  = L_M2(L_M2 >0);
L_M5  = L_M5(L_M5 >0);

figure()
hold on
plot(L_10/mean(L_10))
plot(L_20/mean(L_20))
plot(L_80/mean(L_80))
plot(L_90/mean(L_90))
plot(L_ave/mean(L_ave))
plot(L_M2/mean(L_M2))
plot(L_M5/mean(L_M5))
hold off

edges = 0:0.2:2;
h1 = histcounts(L_10/mean(L_10),edges);
h2 = histcounts(L_20/mean(L_20),edges);
h3 = histcounts(L_80/mean(L_80),edges);
h4 = histcounts(L_90/mean(L_90),edges);
h5 = histcounts(L_ave/mean(L_ave),edges);
h6 = histcounts(L_M2/mean(L_M2),edges);
h7 = histcounts(L_M5/mean(L_M5),edges);
figure()
bar(edges(1:end-1),[h2; h3; h5; h6; h7]')

Best_A = dirlist(L_10  == max(L_10));
Best_B = dirlist(L_20  == max(L_20));
Best_C = dirlist(L_80  == max(L_80));
Best_D = dirlist(L_90  == max(L_90));
Best_E = dirlist(L_ave == max(L_ave));
Best_F = dirlist(L_M2  == max(L_M2));
Best_G = dirlist(L_M5  == max(L_M5));


% figure()
% histogram(L_10)
% figure()
% histogram(L_20)
% figure()
% histogram(L_80)
% figure()
% histogram(L_90)
% figure()
% histogram(L_ave)


donezo = 'True';
%% Save out here, so dirlist,IP/OP, and metrics all match


Predictors = [L_10;L_20;L_80;L_90;L_ave;L_M2;L_M5;WT50;WT75;WT90;GT50;GT75;GT90];
    
weights = unique(IP_weight);
scales = unique(IP_scale);
sz = size(Predictors,1);
IP_map = zeros(length(weights)+1,length(scales)+1,sz);
for i = 1:length(weights)
    for j = 1:length(scales)
        mask = logical(repmat((IP_weight == weights(i)) .* (IP_scale == scales(j)),sz,1));
        data = Predictors(mask);
        data = reshape(data,sz,length(data)/sz);
        IP_map(i,j,1:end) = mean(data');        
    end
end
IP_map = IP_map./max(max(IP_map));
IP_map(IP_map == 0) = 1 ; 

weights = unique(OP_weight);
scales = unique(OP_scale);
OP_map = zeros(length(weights)+1,length(scales)+1,sz);
for i = 1:length(weights)
    for j = 1:length(scales)
        mask = logical(repmat((OP_weight == weights(i)) .* (OP_scale == scales(j)),sz,1));
        data = Predictors(mask);
        data = reshape(data,sz,length(data)/sz);
        OP_map(i,j,1:end) = mean(data');        
    end
end
OP_map = OP_map./max(max(OP_map));
OP_map(OP_map == 0) = 1;
close all
%%
for i = 1:12
    figure()
    pcolor([15,25,32.5,37.5,42.5,47.5,55,65]+20,[15,25,32.5,37.5,42.5,47.5,55,65],IP_map(1:end,1:end,i));
    hold on
    plot(60,40)
end
%%
for i = 3:3:12
    figure()
    pcolor([125,162.5,187.5,212.5],(7.5:5:33),OP_map(1:end,1:end,i)); 
end

% pcolor([125,162.5,187.5,212.5],(7.5:5:33),OP_map(1:end,1:end,1));
  %%      
Predictors = [L_10;L_20;L_80;L_90;L_ave;L_M2;L_M5;WT50;WT75;WT90;GT50;GT75;GT90];

    sorted_L = sort(Likelyhood);
    L_10(i) = sorted_L(round(size(Likelyhood,1)*0.1));
    L_20(i) = sorted_L(round(size(Likelyhood,1)*0.2));
    L_80(i) = sorted_L(round(size(Likelyhood,1)*0.8));
    L_90(i) = sorted_L(round(size(Likelyhood,1)*0.9));
    L_ave(i) = mean(Likelyhood);
    %Median Filters
    L = zeros(321*321,1);
    L(Reconstruction.id)=Likelyhood;
    L = reshape(L,[321,321]);
    L_M2(i) = mean(mean(medfilt2(L)));
    L_M5(i) = mean(mean(medfilt2(L,[5,5])));


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
