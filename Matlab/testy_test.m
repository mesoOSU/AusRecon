% Load example data
mtexdata twins

%segment that data
[grains, ebsd.grainId,ebsd.mis2mean]=calcGrains(ebsd('indexed'),'angle',5*degree);

%remove tiny grains
ebsd(grains(grains.grainSize<3)) = [];
[grains,ebsd.id,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'ang;e',5*degree);

%smooth them
grains = grains.smooth(5);
plot(grains,grains.meanOrientation)

CS = grains.CS;
gB =grains.boundary;

gB_MgMg = gB('Magnesium','Magnesium');
plot(gB_MgMg,gB_MgMg.misorientation.angle./degree,'linewidth',2)

ind = gB_MgMg.misorientation.angle>85*degree & gB_MgMg.misorientation.angle<87*degree;
mori = gB_MgMg.misorientation(ind);

scatter(mori)

% determine the mean of the cluster
mori_mean = mean(mori,'robust')
% determine the closest special orientation relation ship
round2Miller(mori_mean)
