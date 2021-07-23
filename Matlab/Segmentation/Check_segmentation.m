clear all
close all

%Load up data from saved myEBSD file.
load('EBSD/AF_001_Recon.mat','myEBSD')
%load('../../EBSD_Data/AF96_example_Recon.mat','myEBSD')

% Grab the data we need from myEBSD: The original EBSD, the Austenite
% Reconstruction, the CS and SS objects, and the OR
Aus  = myEBSD.Recon.FullEbsd('Austenite');
Mart = flipud(myEBSD.origEbsd('Martensite'));
CS_A = Aus.CS;
CS_M = Mart.CS;
SS   = myEBSD.SS;
OR = myEBSD.OR;
Alex_Var = rem(myEBSD.Variants.IDs-1,24)+1;

% At some point Aus.prop got modified, either to allow side by side
% viewing, or just as a side effect of a calculation. This fixes that.
Mart.prop.y = Mart.prop.y- min(Mart.prop.y);

% What we want to Emulate:
plot(Aus,Aus.orientations)
Plot_Variants( myEBSD,'cbar',0,'mbar',0);
Plot_Blocks( myEBSD,'cbar',0,'mbar',0);
Plot_Packets( myEBSD,'cbar',0,'mbar',0);


% to better understand what is happening, I want to plot the S3, S9, and
% non-special grain boundaries as well, but I want all that in another
% function for mental simplicity.
[Bounds,~,~] = Map_Sigma_Boundaries(Aus,CS_A,CS_M);
clear Bounds

% =========================================================== %
% EVERYTHING ABOVE HERE IS UNCHANGING. Start below this point
% =========================================================== %
% 
% [grains,Aus.grainId,~] = calcGrains(Aus,'angle',5*degree);
% small_grains = grains(grains.grainSize<6);
% Aus(small_grains) = [];
% [grains,Aus.grainId,~] = calcGrains(Aus,'angle',5*degree);
% 
% 
% 
% grains(grains.grainSize<6) = [];
% %[grains,Aus.grainId,~] = calcGrains(Aus,'angle',5*degree);
% 
% 
% Small_grains = grains(grains.grainSize<6);
% Large_grains = grains(grains.grainSize>5);
% Small_grains_points = Aus.grainId(Aus.grainId == Small_grains);
% 
% 
% [grains,Aus.grainId,~] = calcGrains(Aus,'angle',5*degree);
% 
% 
% % remove two pixel grains
% Aus(grains(grains.grainSize<=6)) = [];
% [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',5*degree);
% 
% % smooth them
% grains = grains.smooth(5);
% 
% % visualize the grains
% plot(grains,grains.meanOrientation)
% 
% % store crystal symmetry of Magnesium
% CS = grains.CS;
% 
% 
% 
% 
% 
% 
% 





% Aus and Mart are different sizes when there is retained Austenite in the
% sample. for now, we are going to just fill those points in (this is how
% the current code works) and continue. the plotting function will color
% these pixels white in the final image.
Aus = Aus.fill;
Mart=Mart.fill;


% Get the 24 3x3 rotation matrices from YardleyVarians
YV = YardleyVariants(OR);
% these need to be transposed and converted into MTEX rotations
clear YV_rot YV_ori
for i =1:24
    YV_rot(i) = rotation('matrix',transpose(YV{i}));
end
YV_ori = orientation(YV_rot,CS_M);
% Verify these look right
figure()
plotPDF(YV_ori,Miller(0,0,1,CS_M))

figure()
cmap = myEBSD.Variants.RGB;
caxis([0,25]);
colormap(cmap);
for ii = 1:24
    plotPDF(YV_ori(ii),Miller(0,0,1,CS_M),'MarkerColor',cmap(ii+1,:)) 
    hold on
end

figure()
cmap = myEBSD.Packets.RGB;
caxis([0,5]);
colormap(cmap);
for ii = 1:24
    plotPDF(YV_ori(ii),Miller(0,0,1,CS_M),'MarkerColor',cmap(ceil((ii)/6)+1,:)) 
    hold on
end



%% Setup
Aus_ori = Aus.orientations;
Mart_ori = Mart.orientations;

Martx24 = repmat(Mart_ori,1,24); 
Ausx24 = repmat(Aus_ori,1,24);
fakeEBSD=myEBSD;

%% Method 1 - min(angle(Mart, (Aus * rot)))
choices_1 = Aus_ori * YV_rot;
miso_1    = angle(choices_1, Martx24);
[~,Vars_1]= min(miso_1,[],2);


Vars = flipud(reshape(Vars_1,[321,321]));
Vars(1,1) = 0;
fakeEBSD.Variants.IDs = Vars;
fakeEBSD.Packets.IDs = ceil(Vars/6);
Plot_Variants(fakeEBSD,'mbar',0,'cbar',0);
Plot_Packets(fakeEBSD,'mbar',0,'cbar',0);

%% Method 2 - min(angle(Mart, (Aus * ori)))
choices_2 = Aus_ori * YV_ori;
miso_2    = angle(choices_2, Martx24);
[~,Vars_2]= min(miso_2,[],2);


Vars = flipud(reshape(Vars_2,[321,321]));
Vars(1,1) = 0;
fakeEBSD.Variants.IDs = Vars;
fakeEBSD.Packets.IDs = ceil(Vars/6);
Plot_Variants(fakeEBSD,'mbar',0,'cbar',0);
Plot_Packets(fakeEBSD,'mbar',0,'cbar',0);

%% Method 3 - min(angle(Aus, (Mart * rot)))
choices_3 = Mart_ori * YV_rot;
miso_3    = angle(choices_3, Ausx24);
[~,Vars_3]= min(miso_3,[],2);


Vars = flipud(reshape(Vars_3,[321,321]));
Vars(1,1) = 0;
fakeEBSD.Variants.IDs = Vars;
fakeEBSD.Packets.IDs = ceil(Vars/6);
Plot_Variants(fakeEBSD,'mbar',0,'cbar',0);
Plot_Packets(fakeEBSD,'mbar',0,'cbar',0);

%% Method 4 - min(angle(Aus, (ori*Mart)))
choices_4 = transpose(YV_ori*Mart_ori);
miso_4    = angle(choices_4, Ausx24);
[~,Vars_4]= min(miso_4,[],2);


Vars = flipud(reshape(Vars_4,[321,321]));
Vars(1,1) = 0;
fakeEBSD.Variants.IDs = Vars;
fakeEBSD.Packets.IDs = ceil(Vars/6);
Plot_Variants(fakeEBSD,'mbar',0,'cbar',0);
Plot_Packets(fakeEBSD,'mbar',0,'cbar',0);

%% Method 4 - min(angle(Aus, (ori*Mart)))
choices_4 = transpose(YV_ori*Mart_ori);
miso_4    = angle(choices_4, Ausx24);
[~,Vars_4]= min(miso_4,[],2);


Vars = flipud(reshape(Vars_4,[321,321]));
Vars(1,1) = 0;
fakeEBSD.Variants.IDs = Vars;
fakeEBSD.Packets.IDs = ceil(Vars/6);
Plot_Variants(fakeEBSD,'mbar',0,'cbar',0);
Plot_Packets(fakeEBSD,'mbar',0,'cbar',0);

%% Methods 9-12: try for loop
    V_A2M_1 = zeros(size(Aus,1),24);
    V_A2M_2 = zeros(size(Aus,1),24);
    V_M2A_1 = zeros(size(Aus,1),24);
    V_M2A_2 = zeros(size(Aus,1),24);

for i = 1:size(YV_ori,2)
    A2M_1 = YV_ori(i) * Aus_ori;
    A2M_2 = Aus_ori * YV_ori(i);
    M2A_1 = YV_ori(i) * Mart_ori;
    M2A_2 = Mart_ori * YV_ori(i);
    
    V_A2M_1(:,i) = angle(A2M_1,Mart_ori)./degree;   
    V_A2M_2(:,i) = angle(A2M_2,Mart_ori)./degree;   
    V_M2A_1(:,i) = angle(M2A_1,Aus_ori)./degree;   
    V_M2A_2(:,i) = angle(M2A_2,Aus_ori)./degree;   
end
%Pick the lowest misorientation choice (again, both ways)
[~,V_A] = min(V_A2M_1,[],2);
[~,V_B] = min(V_A2M_2,[],2);
[~,V_C] = min(V_M2A_1,[],2);
[~,V_D] = min(V_M2A_2,[],2);


Vars = flipud(reshape(V_A,[321,321]));
fakeEBSD.Variants.IDs = Vars;
fakeEBSD.Packets.IDs = ceil(Vars/6);
Plot_Variants(fakeEBSD,'mbar',0,'cbar',0);
%Plot_Packets(fakeEBSD,'mbar',0,'cbar',0);


Vars = flipud(reshape(V_B,[321,321]));
fakeEBSD.Variants.IDs = Vars;
fakeEBSD.Packets.IDs = ceil(Vars/6);
Plot_Variants(fakeEBSD,'mbar',0,'cbar',0);
%Plot_Packets(fakeEBSD,'mbar',0,'cbar',0);


Vars = flipud(reshape(V_C,[321,321]));
fakeEBSD.Variants.IDs = Vars;
fakeEBSD.Packets.IDs = ceil(Vars/6);
Plot_Variants(fakeEBSD,'mbar',0,'cbar',0);
%Plot_Packets(fakeEBSD,'mbar',0,'cbar',0);


Vars = flipud(reshape(V_D,[321,321]));
fakeEBSD.Variants.IDs = Vars;
fakeEBSD.Packets.IDs = ceil(Vars/6);
Plot_Variants(fakeEBSD,'mbar',0,'cbar',0);
%Plot_Packets(fakeEBSD,'mbar',0,'cbar',0);

%%

%EVERYTHING ABOVE HERE WORKS NOW







Mart_options = Mart.orientations*YVar_ori;
Ausx24 = repmat(Aus.orientations,1,24);
miso = angle(Mart_options,Ausx24);
[~,VarsA]= min(miso,[],2);
VarsA = flipud(reshape(VarsA,[321,321]));

fakeEBSD=myEBSD;
fakeEBSD.Variants.IDs = VarsA;
Plot_Variants(fakeEBSD)





% Make array of MartxYV, and Aus repeated, then find angles
%Possible_Marts = transpose(YVar_ori*Mart.orientations);
Possible_Marts = Mart.orientations*YVar_ori;
Ausx24 = repmat(Aus.orientations,1,24);
miso = angle(Possible_Marts,Ausx24);
[~,VarsB]= min(miso,[],2);
VarsB(1) = 0;
VarsB = reshape(VarsB,[321,321]);

% Make array of MartxYV, and Aus repeated, then find angles
%Possible_Austs = transpose(YVar_ori*Aus.orientations);
Possible_Austs = transpose(Aus.orientations*YVar_ori);
Martx24 = repmat(Mart.orientations,1,24);
miso = angle(Possible_Austs,Martx24);
[~,VarsC]= min(miso,[],2);
VarsC = reshape(VarsC,[321,321]);


%compare to correct answer
figure(100)
image(reshape(Alex_Var,[321,321]),'CDataMapping','scaled')
figure()
image(VarsA,'CDataMapping','scaled')
figure()
image(VarsB,'CDataMapping','scaled')
figure()
image(VarsC,'CDataMapping','scaled')

% Maybe got Variant transform backwards? try it both ways
%Variant_miso_A2M = zeros(size(Aus_f,1),24);
%Variant_miso_M2A = zeros(size(Aus_f,1),24);
    V_A2M_1 = zeros(size(Aus_f,1),24);
    V_A2M_2 = zeros(size(Aus_f,1),24);
    V_M2A_1 = zeros(size(Aus_f,1),24);
    V_M2A_2 = zeros(size(Aus_f,1),24);

for i = 1:size(Y_Variants,1)
%    Variant_ori_A2M = (rotation('matrix',Y_Variants{i,1}).inv)*Aus_f.orientations;
%    Variant_ori_M2A = (rotation('matrix',Y_Variants{i,1}).inv)*Mart_f.orientations;
%    Variant_ori_A2M = (Y_Variants{i}.inv)*Aus_f.orientations;
%    Variant_ori_M2A = (Y_Variants{i}.inv)*Mart_f.orientations;
%    Variant_miso_A2M(:,i) = angle(Variant_ori_A2M,Mart_f.orientations)./degree;   
%    Variant_miso_M2A(:,i) = angle(Variant_ori_M2A,Aus_f.orientations)./degree;   
    A2M_1 = (Y_Variants(i).inv)*Aus_f.orientations;
    A2M_2 = (Y_Variants(i))*Aus_f.orientations;
    M2A_1 = (Y_Variants(i).inv)*Mart_f.orientations;
    M2A_2 = (Y_Variants(i))*Mart_f.orientations;
    
    V_A2M_1(:,i) = angle(A2M_1,Mart_f.orientations)./degree;   
    V_A2M_2(:,i) = angle(A2M_2,Mart_f.orientations)./degree;   
    V_M2A_1(:,i) = angle(M2A_1,Aus_f.orientations)./degree;   
    V_M2A_2(:,i) = angle(M2A_2,Aus_f.orientations)./degree;   
end

%Pick the lowest misorientation choice (again, both ways)
[~,Variants_A2M] = min(Variant_miso_A2M,[],2);
[~,Variants_M2A] = min(Variant_miso_M2A,[],2);
x_size = size(unique(Aus_f.prop.x),1);
y_size = size(unique(Aus_f.prop.y),1);
sq_var_A2M = reshape(Variants_A2M,x_size,y_size);
sq_var_M2A = reshape(Variants_M2A,x_size,y_size);
 
X = [0:0.5:160];
Y = [0:0.5:160];
%EVERYTHING ABOVE HERE WORKS NOW

%   ===================================================================  %
% % ERIC: Start looking here.
% Figure 9 is an example of just the packet segmentation, no coloration.
% Segmentation looks good to me (as in, Variant boundaries are in the right
% spot)
%   ===================================================================  %
figure(9)
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
hold on
% var_plot = pcolor(X,Y,sq_var_A2M);
var_plot = pcolor(X,Y,sq_var_A2M);
var_plot.EdgeColor = 'none';
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
mtexColorbar
hold off

% Figure 10 and 11 SHOULD be the variants grouped by Packets, but neither
% lines up correctly. Not sure why. One goes from A-to-M, the other M-to-A.
% Noteably, they are both wrong in different ways, which is counter to what
% I thought would happen.
figure(10)
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
hold on
% var_plot = pcolor(X,Y,sq_var_A2M);
var_plot = pcolor(X,Y,floor(sq_var_A2M/6)*6);
var_plot.EdgeColor = 'none';
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
mtexColorbar
hold off

figure(11)
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
hold on
% var_plot = pcolor(X,Y,sq_var_M2A);
var_plot = pcolor(X,Y,floor(sq_var_M2A/6)*6);
var_plot.EdgeColor = 'none';
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
mtexColorbar
hold off

%ERIC: here is where i started coloring things. Ignore the Container.Maps,
%I dont use them anymore, but i'm leaving them there in case i revert back
%for some reason

%Austin note to future Austin: brackets mean cell. semi-colons mean column
% start of thing that should be function
Packet_RGB  = {
    [1 0 0];[1 0 0];[1 0 0];[1 0 0];[1 0 0];[1 0 0];
    [0 1 0];[0 1 0];[0 1 0];[0 1 0];[0 1 0];[0 1 0];
    [0 0 1];[0 0 1];[0 0 1];[0 0 1];[0 0 1];[0 0 1];
    [1 1 0];[1 1 0];[1 1 0];[1 1 0];[1 1 0];[1 1 0];
    [0 0 0]};
Packet_RGB = cellfun(@(x) x./1,Packet_RGB,'UniformOutput',false);
Packet_map = containers.Map([0:24],Packet_RGB);

Block_RGB = {
    [252 122 120]; [252 122 120]; [250 4   0  ]; [250 4   0  ];
    [102 0   0  ]; [102 0   0  ]; [120 252 122]; [120 252 122];
    [4   250 0  ]; [4   250 0  ]; [0   125 2  ]; [0   125 2  ];
    [130 125 253]; [130 125 253]; [2   2   250]; [2   2   250];
    [1   1   100]; [1   1   100]; [220 220 170]; [220 220 170];
    [250 250 2  ]; [250 250 2  ]; [227 176 10 ]; [227 176 10 ];
    [0 0 0];};
Block_RGB = cellfun(@(x) x./255,Block_RGB,'UniformOutput',false);
Block_map = containers.Map([0:24],Block_RGB);


Variant_RGB = {
    [227,97 ,95 ];[255,147,145];[225,4  ,0  ];[255,4  ,0  ];
    [77 ,0  ,0  ];[127,0  ,0  ];[95 ,227,97 ];[145,255,147];
    [4  ,225,0  ];[4  ,255,0  ];[0  ,100,2  ];[0  ,150,2  ];
    [105,100,228];[155,150,255];[2  ,2  ,225];[2  ,2  ,255];
    [1  ,1  ,75 ];[1  ,1  ,125];[195,195,145];[245,245,195];
    [225,225,2  ];[255,255,2  ];[202,151,10 ];[252,201, 10]
    [0 0 0];};
Variant_RGB = cellfun(@(x) x./255,Variant_RGB,'UniformOutput',false);
Variant_map = containers.Map([0:24],Variant_RGB);

PRGB = cell2mat(Packet_RGB);
BRGB = cell2mat(Block_RGB);
VRGB = cell2mat(Variant_RGB);
%end of thing that should be function

blk = (floor((sq_var_A2M-1)/2)+1)*2;
pkt = (floor((sq_var_A2M-1)/6)+1)*6;

%ERIC: Plots below should be maps colored by variant, block, and packet,
%but again, they are wrong. Not sure why anymore.
figure(12)
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
hold on
var_plot = pcolor(X,Y,blk);
var_plot.EdgeColor = 'none';
colormap(BRGB)
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
mtexColorbar
hold off

figure(13)
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
hold on
var_plot = pcolor(X,Y,pkt);
var_plot.EdgeColor = 'none';
colormap(PRGB)
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
mtexColorbar
hold off

figure(14)
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
hold on
var_plot = pcolor(X,Y,sq_var_A2M);
var_plot.EdgeColor = 'none';
colormap(VRGB)
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
mtexColorbar
hold off




sq_var_A2M(sq_var_A2M>6)=0;
figure(14)
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
hold on
var_plot = pcolor(X,Y,sq_var_M2A);
var_plot.EdgeColor = 'none';
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2)
mtexColorbar
hold off


% % Compute variants and corresponding groupoid from euler angles
% ksi=myEBSD.OR
% [Vtmp,~] = YardleyVariants(ksi);
% [Gtmp,~] = GroupoidVariantsCubic(Vtmp);

% Now concatenate the variant and groupoid (misorientation)
% vectors to account for the existance of any twins
%V = vertcat(V,Vtmp);
%G = vertcat(G,Gtmp);

% % Eric Payton (CIV)9:34 AM
% % Otto & Payton J Mater Sci 2011
% % Eric Payton (CIV)9:39 AM
% % DOI:10.1107/S0108767384000246
% % Eric Payton (CIV)9:44 AM
% % 38.94 degrees @<110>
% % Eric Payton (CIV)9:47 AM
% % ori = orientation.byEuler(30*degree,50*degree,10*degree,cs)
% % (26.56,83.62,26.56)



% determine parent grain ori with twins merged in
% This requires fixing, leave be for now
%L = length(Aus_gb_merged.id);
%merged_ori = orientation('Euler',ones(L,1),ones(L,1),ones(L,1),CS_A,CS_A);
%for i = 1:size(Aus_gb_merged,1)
%    [~,index] = max(Aus_grains(Aus_parentId== i).grainSize);
%    merged_ori(i) = Aus_grains(Aus_parentId== i).meanOrientation(index); 
%end
%mergedGrains.meanOrientation = merged_ori;
%clear i index merged_ori
%figure(8)
%plot(mergedGrains,mergedGrains.meanOrientation)

% %%Austin Note: I am pretty positive this works and does simplify things
% %% a bit, but I am ignoring for now to avoid symmeterizing stuff
% %% Steve calc_m2A
% %Since defining explicitly as misorientation and giving symm, only need 1 to get all 24
% function [A2M,flag]=calc_A2M(ksi,CS_A,CS_M);
% [OR,flag]=YardleyVariants(ksi);%[3.3,8.5,8.9]);
% A2M=orientation(‘matrix’,(OR{1}’),CS_M,CS_A);
% %M2A=orientation(‘matrix’,OR{1},CS_A,CS_M);
% end
% 
% function [M2A,flag]=calc_M2A(ksi,CS_A,CS_M);
% [OR,flag]=YardleyVariants(ksi);%[3.3,8.5,8.9]);
% %A2M=orientation(‘matrix’,(OR{1}’),CS_M,CS_A);
% M2A=orientation(‘matrix’,OR{1},CS_A,CS_M);
% end
% 
% %%
% 
% 
%Not relevant, but helpful reminder of how to do inline function def
%Normalize_block = @(x) x./255
%Block_map = cellfun(Normalize_block,Block_RGB,'UniformOutput',false)
