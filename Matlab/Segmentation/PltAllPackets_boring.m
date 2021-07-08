function PltAllPackets_boring(varargin)
% Plots packet boundaries for the entire microstructure using the same
% coloring convection as described below.
    
    % Assign myEBSD structure
    myEBSD = varargin{1};
    
    % Remove micronbar if desired
    if length(varargin)==2
        mbar = 'MicronBar';
        mkey = 'off';
    else
        mbar = [];
        mkey = [];
    end
    
    % Austenite Reconstruction
    AusRecon = myEBSD.Recon.Ebsd;
    
    % Phase IDs
    TransID = myEBSD.Phase.ID{1};
    ReconID = myEBSD.Phase.ID{2};
    
    % Post-transformation martensite structure
    Martensite = myEBSD.Ebsd(myEBSD.Ebsd.phase==TransID);
    MartOrs = Martensite.orientations;
    
    % Remove indexing for initially indexed austenite so only the
    % martensite points will actually be plotted
    OrigEbsd = myEBSD.Ebsd;
    InitAus = find(OrigEbsd.phase==ReconID);
    
        
    % Orientation color mapping key
    key = ipfHSVKey(Martensite);
    key.inversePoleFigureDirection=zvector;
    mc = key.orientation2color(MartOrs);
    
    % Specific packet indices
    PckInds = myEBSD.Packets.IDs;
    PckInds(InitAus)=[];
    for i = 1:5
        Pck{i} = PckInds == i-1;
    end

    % Plot packets
    figure; plot(Martensite,mc,mbar,mkey);
    hold on
    
    % Adjust colormap values
    % cmap = myEBSD.Packets.RGB;
    caxis([0,5]);
    % colormap(cmap);
    for ii = 1:5
        tmp = Pck{ii};
        if isempty(tmp)==0
            mp = key.orientation2color(MartOrs(tmp));
            plot(Martensite(tmp),mp);%,'FaceColor',cmap(ii,:));
            hold on
        end
    end
    cm = colorbar;
    set(cm, 'Ticks', 0.5:4.5, 'TickLabels',[0:4]);
    
%     % Plot the associated packet weights
%     Wts = myEBSD.Packet.Wts;
%     Wts(InitAus)=[];
%     figure; plot(AusRecon,Wts,mbar,mkey);

end

