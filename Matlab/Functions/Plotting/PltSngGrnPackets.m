function PltSngGrnPackets(myEBSD,Packets,grnId)
% Plots packet boundaries in a clear fashion along with the austenite 
% grain, martensite grains within the PAG, and weights if flagged. Note 
% that the martensite grains are plotted with a micronbar, but every other 
% choice is up to the user. 
    
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
    
    % Specific packet indices
    PckInds = Packet.Boundaries;
    for i = 1:5
        Pck{i} = PckInds == i-1;
    end

    % Plot packet boundaries
    figure; plot(MartGrns,mc,'MicronBar','off');
    hold on
    plot(Packet.Grain.boundary);
    hold on
    % Adjust colormap values
    cmap = myEBSD.Packets.RGB;
    caxis([0,5]);
    colormap(cmap);
    for ii = 1:5
        tmp = Pck{ii};
        if isempty(tmp)==0
            mp = key.orientation2color(PackOrs(tmp));
            plot(AusGrn(tmp),mp,'FaceColor',cmap(ii,:));
            hold on
        end
    end
    cm = colorbar;
    set(cm, 'Ticks', 0.5:4.5, 'TickLabels',[0:4]);
    
    % Associated packet weights
    Wts = Packet.Weights;
    % Plot weights  
    figure; plot(AusGrn,Wts,'MicronBar','off');
    hold on
    plot(Packet.Grain.boundary);

end

