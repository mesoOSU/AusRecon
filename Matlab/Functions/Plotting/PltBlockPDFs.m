function PltPacketPDFs(myEBSD,Packets,grnId,Dirs)
    % Plot pole figures based on user-input pole axis and orientation data.
    % Coloring used for variants is consistent with other block segment
    % plotting factors
        
    Packet = Packets{grnId};
    Marker = 10;
    
    OR = myEBSD.OR;

    % Prior austenite grain and corresponding variants
    PAG = Packet.AusGrain;
    PAG_Vars = Packet.MartGrains;

    % Austenite/Martensite orientation(s) and CS
    AusOr = PAG(1).orientations;
    CS_A = PAG.CS;
    AusOr = rotation('euler',AusOr.phi1,AusOr.Phi,AusOr.phi2,CS_A);
    MartGrnOrs = PAG_Vars.orientations;
    CS_M = PAG_Vars.CS;

    % Orientation color mapping key
    key = ipfHSVKey(AusOr);
    key.inversePoleFigureDirection=zvector;

    len = length(MartGrnOrs);
    if len > 2000
        permnums = randperm(len);
        nums = permnums(1:2000);
    else
        nums = [1:len];
    end
    MartOrs = MartGrnOrs(nums);
    Ors = inv(AusOr) * MartOrs;

    % Pole Figure Directions
    if isempty(Dirs)
         PFdir = Miller({0,0,1},CS_M);
    else
        if length(Dirs(:)) == 3
            PFdir = Miller({Dirs(1) Dirs(2) Dirs(3)},CS_M);
        elseif length(Dirs(:))==6
            PFdir = Miller({Dirs(1) Dirs(3) Dirs(5)},{Dirs(2) Dirs(4) Dirs(6)},CS_M);
        elseif length(Dirs(:))==9
            PFdir = ...
                Miller({Dirs(1) Dirs(4) Dirs(7)},{Dirs(2) Dirs(5) Dirs(8)},{Dirs(3) Dirs(6) Dirs(9)},CS_M);
        else
            warning('Max of 3 pole directions. Defaulting to {001} direction')
            PFdir = Miller({0,0,1},CS_M); 
        end
    end

    BlckInds = Packet.Block.Boundaries;
    PDF_Inds = BlckInds(nums);

    for i = 1:13
        Blck{i} = PDF_Inds == i-1;
    end
    
    figure;
    % Adjust colormap values
    cmap = myEBSD.Blocks.RGB;
    caxis([0,13]);
    colormap(cmap);
    for ii = 1:13
        tmp = Blck{ii};
        if isempty(tmp)==0
            plotPDF(Ors(tmp),PFdir,'antipodal','points','all','marker','.','MarkerColor',cmap(ii,:),'MarkerSize',Marker) 
            hold on
        end
    end

% Now plot theoretical variants
    % Compute variants and corresponding groupoid from euler angles
    [V,flag] = YardleyVariants(OR);
    [G,CT] = GroupoidVariantsCubic(V);
    
    % Since result is a passive rotation, take inverse of variant
    % transformation indices to convert to active (Consistent with MTEX)
    % NOTE TOREN: Cryspy is all passive, so skip this transpose step when 
    % converting reconstruction code for cryspy
    for j = 1:24
        Vt{j,1} = transpose(V{j});
    end
    
    % Now transform the transposed variants to euler angles
    for jj = 1:24
        VTeul(jj) = orientation('matrix',Vt{jj},CS_M);
    end
    
    figure;
    % Reset colormap for this figure
    colormap(cmap);
    caxis([0,12])
    
    for k = 1:12
        c1 = (k-1)*2+1;
        c2 = k*2;
        plotPDF(VTeul(c1:c2),PFdir,'antipodal','points','all','marker','.','MarkerColor',cmap(k+1,:),'MarkerSize',10) 
        hold on
    end        
end


