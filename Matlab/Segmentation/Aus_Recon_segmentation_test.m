clear all
close all

%Load up data from saved myEBSD file.
load('EBSD/AF_001_Recon.mat','myEBSD')

%grab the reconstructed map of the prior orientations (no identificataion
%of twin vs grain, just "what was i RIGHT before the transform")
Aus  = myEBSD.Recon.Ebsd('Austenite');
Mart = flipud(myEBSD.origEbsd('Martensite'));
CS_A = Aus.CS;
CS_M = Mart.CS;
SS   = myEBSD.SS;
OR = myEBSD.OR;
% clear myEBSD

%segment the Aus map and create a GB structure (this step treats Twin and
%prior grians the same way, S3 and s9 ID come later
[Aus_grains,~,~] = calcGrains(Aus,'angle',5*degree);
Aus_grains = Aus_grains.smooth(5);
Aus_gb = Aus_grains.boundary('Austenite','Austenite');



%Plot prior Aus
figure(1)
plot(Aus,Aus.orientations);
%Plot prior Aus filled
figure(2)
plot(Aus.fill,Aus.fill.orientations);
% visualize the grains
% figure(3)
% plot(Aus_grains,Aus_grains.meanOrientation)
%plot GB alone, colored by misorientation angle
% figure(4)
% hold on
% plot(Aus_gb,Aus_gb.misorientation.angle./degree,'linewidth',2)
% mtexColorbar
% hold off
% %Histogram maker (not needed
% figure(5)
% hold on
% histogram(Aus_gb.misorientation.angle./degree,40)
% xlabel('misorientation angle (degree)')
% hold off
% code to show location of S3 and s9 in rodruiguez space
ind = Aus_gb.misorientation.angle>59.5*degree & Aus_gb.misorientation.angle<60.5*degree;
mori = Aus_gb.misorientation(ind);
%determine the mean of the cluster
% mori_mean = mean(mori,'robust')
% %determine the closest special orientation relation ship
% round2Miller(mori_mean)
% figure(6)
% hold on
% scatter(mori)
% hold off


% define S3 and S9. 
%S3 is defined using u1v1 method from Steve
u1 = Miller( 1, 2, 1,CS_A);
v1 = Miller(-1,-2,-1,CS_A);
u2 = Miller(-1, 0, 1,CS_A);
v2 = Miller( 1, 0,-1,CS_A);
S3_twin_miso = orientation('map',u1,v1,u2,v2);
% S9 was done from 1984 Hans Grimmer Paper
S9_angle = 38.94*degree;
S9_axis = vector3d(1,1,0);
S9_twin_miso = orientation('axis',S9_axis,'angle',S9_angle,CS_A);
S9_twin_miso.SS = S9_twin_miso.CS;

% restrict to twinnings with threshold 0.5 degree
S3_isTwinning = angle(Aus_gb.misorientation,S3_twin_miso) < 0.5*degree;
S9_isTwinning = angle(Aus_gb.misorientation,S9_twin_miso) < 0.5*degree;

S3_TwinBoundaries = Aus_gb(S3_isTwinning);
S9_TwinBoundaries = Aus_gb(S9_isTwinning);
All_Twin_Boundaries = [S9_TwinBoundaries S3_TwinBoundaries];

[Aus_gb_merged,Aus_parentId] = merge(Aus_grains,All_Twin_Boundaries);

% plot the twinning boundaries
figure(7)
plot(Aus_grains,Aus_grains.meanOrientation)
hold on
plot(All_Twin_Boundaries,'linecolor','w','linewidth',6)
plot(S3_TwinBoundaries,'linecolor','r','linewidth',2.5,'displayName','S3 Twin')
plot(S9_TwinBoundaries,'linecolor','b','linewidth',2.5,'displayName','S9 Twin')
plot(Aus_gb_merged.boundary,'linecolor','k','linewidth',2,'displayName','Prior')
hold off

Aus_f  = Aus.fill;
Mart_f = Mart.fill;
Mart_f.prop.y = Mart_f.prop.y*-1;
prop = Mart_f.prop;
Aus_f.prop = prop;

% Get the Variant rotations based on OR from Yardley function
%threebythree_Y_Var = YardleyVariants(OR);
%for i =1:24
%    y_var[i]=rotation('matrix',threebythree_Y_Var
%end
%Y_Variants


YV = YardleyVariants(OR);

clear Y_Variants
for i = 1:24
    Y_Variants(i) = rotation('matrix',YV{i});
end

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
% A2M=orientation(???matrix???,(OR{1}???),CS_M,CS_A);
% %M2A=orientation(???matrix???,OR{1},CS_A,CS_M);
% end
% 
% function [M2A,flag]=calc_M2A(ksi,CS_A,CS_M);
% [OR,flag]=YardleyVariants(ksi);%[3.3,8.5,8.9]);
% %A2M=orientation(???matrix???,(OR{1}???),CS_M,CS_A);
% M2A=orientation(???matrix???,OR{1},CS_A,CS_M);
% end
% 
% %%
% 
% 
%Not relevant, but helpful reminder of how to do inline function def
%Normalize_block = @(x) x./255
%Block_map = cellfun(Normalize_block,Block_RGB,'UniformOutput',false)
