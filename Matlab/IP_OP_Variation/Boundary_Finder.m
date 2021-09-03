function [Bounds,Aus_grains,Aus_parentId] = Boundary_Finder(Aus,Mart)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[Aus_grains,~] = calcGrains(Aus,'angle',5*degree);
Aus_grains = Aus_grains.smooth(5);
Aus_gb = Aus_grains.boundary('Austenite','Austenite');
%CS_A = Aus('Austenite').CS;
CS_M = Mart.CSList{1,3};
CS_A = Mart.CSList{1,2};
% define S3 and S9. 
%S3 is defined using u1v1 method from Steve
u1 = Miller( 1, 2, 1,CS_A);
v1 = Miller(-1,-2,-1,CS_A);
u2 = Miller(-1, 0, 1,CS_M);
v2 = Miller( 1, 0,-1,CS_M);
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

Bounds = [];
Bounds.S3_Twin = S3_TwinBoundaries;
Bounds.S9_Twin = S9_TwinBoundaries;
Bounds.All = Aus_gb_merged.boundary;

end

