clear all
close all

%load in data from a half finished run of AF96_Example.ang
myEBSD = load('Line_228.mat',['myEBSD']);
Parent = load('Line_228.mat',['Parent']);
Twin = load('Line_228.mat',['Twin']);

% Fix naming
myEBSD=myEBSD.myEBSD;
Parent = Parent.Parent;
Twin = Twin.Twin;

% grab JUST the reconstructed prior EBSD map
Recon_ebsd = myEBSD.Recon.Ebsd;
Recon_ebsd =Recon_ebsd('Austenite');
plot(Recon_ebsd,Recon_ebsd.orientations);

% segment grains,smooth them, and create GB structure
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(Recon_ebsd,'angle',5*degree);
grains = grains.smooth(5);
CS = grains.CS;
gB = grains.boundary;
gB_Aus = gB('Austenite','Austenite');

% visualize the grains
%figure(1)
%plot(grains,grains.meanOrientation)

%plot GB alone, colored by misorientation angle
figure(2)
hold on
plot(gB_Aus,gB_Aus.misorientation.angle./degree,'linewidth',2)
mtexColorbar
hold off

%Histogram maker (not needed
% figure(3)
% hold on
% histogram(gB_Aus.misorientation.angle./degree,40)
% xlabel('misorientation angle (degree)')
% hold off

% code to show location of S3 and s9 in rodruiguez space
% ind = gB_Aus.misorientation.angle>59.5*degree & gB_Aus.misorientation.angle<60.5*degree;
% mori = gB_Aus.misorientation(ind);
% determine the mean of the cluster
% mori_mean = mean(mori,'robust')
% determine the closest special orientation relation ship
% round2Miller(mori_mean)
% figure(4)
% hold on
% scatter(mori)
% hold off


% define S3 and S9. 
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

% restrict to twinnings with threshold 0.5 degree
S3_isTwinning = angle(gB_Aus.misorientation,S3_twin_miso) < 0.5*degree;
S9_isTwinning = angle(gB_Aus.misorientation,S9_twin_miso) < 0.5*degree;

S3_TwinBoundaries = gB_Aus(S3_isTwinning);
S9_TwinBoundaries = gB_Aus(S9_isTwinning);
All_Twin_Boundaries = [S9_TwinBoundaries S3_TwinBoundaries];

[mergedGrains,parentId] = merge(grains,All_Twin_Boundaries);

% plot the twinning boundaries
figure(3)
plot(grains,grains.meanOrientation)
hold on
plot(All_Twin_Boundaries,'linecolor','w','linewidth',6)
plot(S3_TwinBoundaries,'linecolor','r','linewidth',2.5,'displayName','S3 Twin')
plot(S9_TwinBoundaries,'linecolor','b','linewidth',2.5,'displayName','S9 Twin')
plot(mergedGrains.boundary,'linecolor','k','linewidth',2,'displayName','Prior')
hold off

% determine parent grain ori with twins merged in
merged_ori = orientation('Euler',ones(22,1),ones(22,1),ones(22,1),CS,CS);
for i = 1:size(mergedGrains,1)
    [~,index] = max(grains(parentId== i).grainSize);
    merged_ori(i) = grains(parentId== i).meanOrientation(index); 
end
mergedGrains.meanOrientation = merged_ori;
clear i index merged_ori

%figure(6)
%plot(mergedGrains,mergedGrains.meanOrientation)


orig_EBSD = flipud(myEBSD.origEbsd).fill;
Recon_EBSD = myEBSD.Recon.Ebsd.fill;
orig_EBSD.prop.y = orig_EBSD.prop.y*-1;
prop = orig_EBSD.prop;
Recon_EBSD.prop = prop;
Aus = Recon_EBSD('Austenite').fill.orientations;
Mart = orig_EBSD('Martensite').fill.orientations;


% Recon_ebsd = myEBSD.Recon.Ebsd.fill
Y_Variants = YardleyVariants(myEBSD.OR);

% Maybe got Variant transform backwards? try it both ways

OR = YardleyVariants(myEBSD.OR);
A2M = orientation('matrix',OR{1}',CS,CS);
id_Mart_variants = A2M*symmetrise(orientation(idquaternion,CS));



Variant_miso_A2M = zeros(size(Aus,1),24);
Variant_miso_M2A = zeros(size(Aus,1),24);
for i = 1:size(Y_Variants,1)
    Variant_ori_A2M = transpose(rotation('matrix',Y_Variants{i,1})*Aus);
    Variant_ori_M2A = transpose(rotation('matrix',Y_Variants{i,1})*Mart);

    Variant_miso_A2M(:,i) = angle(Variant_ori_A2M,Mart)./degree;   
    Variant_miso_M2A(:,i) = angle(Variant_ori_M2A,Aus)./degree;   
end


[~,Variants_A2M] = min(Variant_miso_A2M,[],2);
[~,Variants_M2A] = min(Variant_miso_M2A,[],2);
x_size = size(unique(Recon_ebsd.prop.x),1);
y_size = size(unique(Recon_ebsd.prop.y),1);
sq_var_A2M = reshape(Variants_A2M,x_size,y_size);
sq_var_M2A = reshape(Variants_M2A,x_size,y_size);

Variant_RGB = zeros(size(Aus,1),24,3);

figure(4)
image(sq_var_A2M,'CDataMapping','scaled')
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
% Packet_RGB  = [[1 0 0];[1 0 0];[1 0 0];[1 0 0];[1 0 0];[1 0 0];
%                [0 1 0];[0 1 0];[0 1 0];[0 1 0];[0 1 0];[0 1 0];
%                [0 0 1];[0 0 1];[0 0 1];[0 0 1];[0 0 1];[0 0 1];
%                [1 1 0];[1 1 0];[1 1 0];[1 1 0];[1 1 0];[1 1 0]]';
%            
% Packet_map = containers.Map([1:24],Packet_RGB);
% 
% 
% 
% 
% 
% 
% 
% 
% Block_map = containers.Map
% 
% Block_RGB             ={[252 122 120],[252 122 120],
%                         [250 4   0  ],[250 4   0  ],
%                         [102 0   0  ],[102 0   0  ],
%                         [120 252 122],[120 252 122],
%                         [4   250 0  ],[4   250 0  ],
%                         [0   125 2  ],[0   125 2  ],
%                         [130 125 253],[130 125 253],
%                         [2   2   250],[2   2   250],
%                         [1   1   100],[1   1   100],
%                         [220 220 170],[220 220 170],
%                         [250 250 2  ],[250 250 2  ],
%                         [227 176 10],[227 176 10],};
% Block_RGB = reshape(Block_RGB,[1,24]);
% Block_RGB = cellfun(@(x) x./255,Block_RGB,'UniformOutput',false);
% Block_map = containers.Map([1:24],Block_RGB);
% 
% % Normalize_block = @(x) x./255
% % Block_map = cellfun(Normalize_block,Block_RGB,'UniformOutput',false)
% 
% 
% % Compute variants and corresponding groupoid from euler angles
% ksi=myEBSD.OR
% [Vtmp,~] = YardleyVariants(ksi);
% [Gtmp,~] = GroupoidVariantsCubic(Vtmp);
% 
% % Now concatenate the variant and groupoid (misorientation)
% % vectors to account for the existance of any twins
% %V = vertcat(V,Vtmp);
% %G = vertcat(G,Gtmp);
% 
% 
% 
% 
% % Eric Payton (CIV)9:34 AM
% % Otto & Payton J Mater Sci 2011
% % Eric Payton (CIV)9:39 AM
% % DOI:10.1107/S0108767384000246
% % Eric Payton (CIV)9:44 AM
% % 38.94 degrees @<110>
% % Eric Payton (CIV)9:47 AM
% % ori = orientation.byEuler(30*degree,50*degree,10*degree,cs)
% % (26.56,83.62,26.56)