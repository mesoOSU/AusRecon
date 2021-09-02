function Plot_Variants(myEBSD,varargin)
% Plot Variants. Possible options:
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

% Grab everything we need from the original structure: The Variant IDs,
% the values for which phase is the Pre and Post transformation, the
% EBSD map of the Martensite points, the the RGB colorkey
Vars = myEBSD.Variants.IDs;
ReconID = myEBSD.Phase.ID{2};
TransID = myEBSD.Phase.ID{1};
Martensite = myEBSD.Ebsd(myEBSD.Ebsd.phase==TransID);
Var_cmap = myEBSD.Variants.RGB;

% Instead of ALL the Vars, we want just the Vars that correspond to
% voxels of the transformed (Martensite) phase.
Vars(myEBSD.Ebsd.phase==ReconID) = [];
assert(length(Vars) ==length(Martensite),...
    "Vars is not the same size as the list of martensite points. Check that the Transformation and Reconstruction Phases have been correctly identified")

% AusRecon assigns Var ID uniquely for different twins. Floor divide to fix.
Vars = rem(Vars-1,24)+1;

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
caxis([0,25]);
colormap(Var_cmap);

% Assign FaceColor to each respective variant
for ii = 1:25
    plot(Martensite(Vars == ii-1),mp(Vars == ii-1),'FaceColor',Var_cmap(ii,:));
end

% ================ %
%  Optional Stuff
% ================ %
% Display colorbar
if cbar == 1
    cm = colorbar;
    set(cm, 'Ticks', 0.5:24.5, 'TickLabels',[0:24]);
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

