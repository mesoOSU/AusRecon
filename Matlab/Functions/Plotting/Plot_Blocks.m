function Plot_Blocks(myEBSD,varargin)
% Plot Blocks. Possible options:
% mbar      : True/False, turns Micronbar on and off
% key       : True/False, turns color key on and off
% Austenite : True/False, turns Austenite border drawing on and off
% Packet    : True/False, turns Packet border drawing on and off

% set all options to true by default, then update from varargin
p = inputParser;
addRequired(p,'myEBSD')
addOptional(p,'mbar',1)
addOptional(p,'cbar',1)
addOptional(p,'Austenite',1)
addOptional(p,'Packet',1)
parse(p,myEBSD,varargin{:})

mbar = p.Results.mbar;
cbar = p.Results.cbar;
Aus_border = p.Results.Austenite;
Packet_border = p.Results.Packet;

% Grab everything we need from the original structure: The Block IDs,
% the values for which phase is the Pre and Post transformation, the
% EBSD map of the Martensite points, the the RGB colorkey
Blocks = myEBSD.Blocks.IDs;
ReconID = myEBSD.Phase.ID{2};
TransID = myEBSD.Phase.ID{1};
Martensite = myEBSD.Ebsd(myEBSD.Ebsd.phase==TransID);
Block_cmap = myEBSD.Blocks.RGB;

% Instead of ALL the Blocks, we want just the Blocks that correspond to
% voxels of the transformed (Martensite) phase.
Blocks(myEBSD.Ebsd.phase==ReconID) = [];
assert(length(Blocks) ==length(Martensite),...
    "Blocks is not the same size as the list of martensite points. Check that the Transformation and Reconstruction Phases have been correctly identified")

% futureproofing in case Blocks return to being Twin-specific.
Blocks = rem(Blocks-1,12)+1;

% create a key to explicitly convert orientations to RGB values.We will
% write over these colors later, but we have to plot something
% initially before recoloring, and adding this in helps with error
% checking. Also, I think explicitly plotting RGB values flips some
% internal MTEX switch to interpret 3xN stringas as RGB, which is what
% we want.
key = ipfHSVKey(Martensite); % IPF2RGB key
key.inversePoleFigureDirection=zvector; % set view direction
mp = key.orientation2color(Martensite.orientations); %ori, but in Nx3 RGB

% Make initial plot
figure;
if mbar == 1
    plot(Martensite,mp,mbar); %EBSD plot
else
    plot(Martensite,mp,'micronbar','off')
end

hold on
caxis([0,13]);
colormap(Block_cmap);

% Assign FaceColor to each respective variant
for ii = 1:13
    plot(Martensite(Blocks == ii-1),mp(Blocks == ii-1),'FaceColor',Block_cmap(ii,:));
end

% ================ %
%  Optional Stuff
% ================ %
% Display colorbar
if cbar == 1
    cm = colorbar;
    set(cm, 'Ticks', 0.5:12.5, 'TickLabels',[0:12]);
end
% Add Packet Border
if Packet_border == 1
    plot(myEBSD.Packets.Boundaries.boundary,'FaceColor','black')
end
% Add Prior Austenite Border
if Aus_border == 1
    plot(myEBSD.AusGrains.grains.boundary,'FaceColor','white')
end
end
