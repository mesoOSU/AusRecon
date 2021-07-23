function Plot_Microstructure_Hierarchy(myEBSD,varargin)
% Plot Variants. Possible options:
% mbar      : True/False, turns Micronbar on/off
% key       : True/False, turns color key on/off
% A_outline : True/False, turns Austenite border drawing on/off
% P_outline : True/False, turns Packet border drawing on/off
% B_outline : True/False, turns Block border drawing on/off
% V_outline : True/False, turns Variant border drawing on/off
% P         : True/False, makes Packet Plot
% B         : True/False, makes Block Plot
% V         : True/False, makes Variant Plot

% set all options to true by default, then update from varargin
p = inputParser;
addRequired(p,'myEBSD')
addOptional(p,'mbar',1)
addOptional(p,'cbar',1)
addOptional(p,'A_outline',1)
addOptional(p,'P_outline',1)
addOptional(p,'B_outline',0)
addOptional(p,'V_outline',0)
addOptional(p,'A',1)
addOptional(p,'B',1)
addOptional(p,'P',1)
addOptional(p,'V',1)
parse(p,myEBSD,varargin{:})

% Grab data that requires altering before viewing
Vars = myEBSD.Variants.IDs;
Martensite = myEBSD.Ebsd(myEBSD.Ebsd.phase==myEBSD.Phase.ID{1});
% Instead of ALL the Vars, we want just the Vars that correspond to
% voxels of the transformed (Martensite) phase.
Vars(myEBSD.Ebsd.phase==myEBSD.Phase.ID{2}) = [];
assert(length(Vars) ==length(Martensite),...
    "Vars is not the same size as the list of martensite points. Check that the Transformation and Reconstruction Phases have been correctly identified")
% AusRecon assigns Var ID uniquely for different twins. Floor divide to fix.
Vars = rem(Vars-1,24)+1;
% create a key to explicitly convert orientations to RGB values.We will
% overwrite these later, but this is the fastest way to initialize the plot
key = ipfHSVKey(Martensite); % IPF2RGB key
key.inversePoleFigureDirection=zvector; % set view direction
mp = key.orientation2color(Martensite.orientations); %ori, but in Nx3 RGB
%Variant Plot
if p.Results.V == true
    figure();
    plot(Martensite,mp,'micronbar',p.Results.mbar)
    hold on
    caxis([0,25]);
    colormap(myEBSD.Variants.RGB);
    for ii = 1:25
        plot(Martensite(Vars == ii-1),mp(Vars == ii-1),'FaceColor',myEBSD.Variants.RGB(ii,:));
    end
    plot_options(myEBSD,p,24.5)
end
%Block Plot
if p.Results.B == true
    Blocks = ceil(Vars/2);
    figure();
    plot(Martensite,mp,'micronbar',p.Results.mbar)
    hold on
    caxis([0,13]);
    colormap(myEBSD.Blocks.RGB);
    for ii = 1:13
        plot(Martensite(Blocks == ii-1),mp(Blocks == ii-1),'FaceColor',myEBSD.Blocks.RGB(ii,:));
    end
    plot_options(myEBSD,p,12.5)
end
%Packets Plot
if p.Results.P == true
    Packets = ceil(Vars/6);
    figure();
    plot(Martensite,mp,'micronbar',p.Results.mbar)
    hold on
    caxis([0,5]);
    colormap(myEBSD.Packets.RGB);
    for ii = 1:5
        plot(Martensite(Packets == ii-1),mp(Packets == ii-1),'FaceColor',myEBSD.Packets.RGB(ii,:));
    end
    plot_options(myEBSD,p,4.5)
end
end

function plot_options(myEBSD,p,c_max)
% ================ %
%  Optional Stuff
% ================ %
% Display colorbar
if p.Results.cbar == 1
    cm = colorbar;
    set(cm, 'Ticks', 0.5:c_max, 'TickLabels',[0:24]);
end
% Add Variant Border
if p.Results.V_outline == true
    plot(myEBSD.Variants.SubBlockBoundaries.boundary,'FaceColor','magenta','linewidth',0.5)
end
% Add Block Border
if p.Results.B_outline == true
    plot(myEBSD.Blocks.Boundaries.boundary,'FaceColor','cyan','linewidth',1.0)
end
% Add Packet Border
if p.Results.P_outline == true
    plot(myEBSD.Packets.Boundaries.boundary,'FaceColor','black','linewidth',1.5)
end
% Add Prior Austenite Border
if p.Results.A_outline == true
    plot(myEBSD.AusGrains.grains.boundary,'FaceColor','white','linewidth',2)
end
end
