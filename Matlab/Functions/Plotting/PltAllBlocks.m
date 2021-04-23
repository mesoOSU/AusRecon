function PltAllBlocks(varargin)
% Plots packet boundaries for the entire microstructure using the same
% coloring convection as described below.
    % Black  --> No Assigned Packet!
    myEBSD = varargin{1};
    
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
    
    % Post-transformation martensite structure and homogenize it
    Martensite = myEBSD.Ebsd(find(myEBSD.Ebsd.phase==TransID));
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
    BlckInds = myEBSD.Blocks.IDs;
    BlckInds(InitAus)=[];
    for i = 1:13
        Blck{i} = BlckInds == i-1;
    end

    % Plot packet boundaries
    figure; plot(Martensite,mc,mbar,mkey);
    hold on
    
    % Adjust colormap values
    cmap = myEBSD.Blocks.RGB;
    caxis([0,13]);
    colormap(cmap);
    for ii = 1:13
        tmp = Blck{ii};
        if isempty(tmp)==0
            mp = key.orientation2color(MartOrs(tmp));
            plot(Martensite(tmp),mp,'FaceColor',cmap(ii,:));
            hold on
        end
    end
    % Display colorbar
    cm = colorbar;
    set(cm, 'Ticks', 0.5:12.5, 'TickLabels',[0:12]);
    
    % Overlay packet boundaries ontop of the block-ID map
    PckBounds = myEBSD.Packets.Boundaries;
    
    hold on
    plot(PckBounds.boundary)
    
%     % Plot the associated packet weights
%     Wts = myEBSD.Packet.Wts;
%     Wts(InitAus)=[];
%     figure; plot(AusRecon,Wts,mbar,mkey);

end

