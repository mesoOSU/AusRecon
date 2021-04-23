function PltSngGrnBlocks(myEBSD,Packets,grnId)
% Plots packet boundaries in a clear fashion along with the austenite 
% grain, martensite grains within the PAG, and weights if flagged. Note 
% that the martensite grains are plotted with a micronbar, but every other 
% choice is up to the user. The color scheme of the packet boundaries is:
    % Black  --> No Assigned Packet!
    Colors = {'black'};
%     cmap = colormap('hsv');
    
    % Extract individual packet structure
    Packet = Packets{grnId};
    
    % Austenite grain
    AusGrn = Packet.AusGrain;
    
    % Orientation color mapping key
    key = ipfHSVKey(AusGrn);
    key.inversePoleFigureDirection=zvector;
    
    % Packet orientations
    PackOrs = Packet.MartGrains.orientations;
    
    % Martensite grains
    MartGrns = Packet.MartGrains;
    mc = key.orientation2color(PackOrs);
    % Plot the martensite variants transformed from the PAG
    figure; plot(MartGrns,mc,'MicronBar','off')
    hold on
    plot(Packet.Grain.boundary);
    
    % Specific Block indices
    BlckInds = Packet.Block.Boundaries;
    for i = 1:13
        Blck{i} = BlckInds == i-1;
    end

    % Plot block boundaries    
    figure; plot(MartGrns,mc,'MicronBar','off')
    hold on
    plot(Packet.Grain.boundary);
    hold on
    
    % Adjust colormap values
    cmap = myEBSD.Blocks.RGB;
    caxis([0,13]);
    colormap(cmap);
    for ii = 1:13
        tmp = Blck{ii};
        if isempty(tmp)==0
            mp = key.orientation2color(PackOrs(tmp));
            plot(AusGrn(tmp),mp,'FaceColor',cmap(ii,:));
            hold on
        end
    end
    % Display colorbar
    cm = colorbar;
    set(cm, 'Ticks', 0.5:12.5, 'TickLabels',[0:12]);
    
    % Associated packet weights
    Wts = Packet.Block.Weights;
    % Plot weights  
    figure; plot(AusGrn,Wts,'MicronBar','off');
    hold on
    plot(Packet.Grain.boundary);

end

