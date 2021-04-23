function [myEBSD] = CharMartBoundaries(myEBSD) 
    % Characterizes the packet, block, and sub-block boundaries determined
    % by the reconstruction algorithm and saves them in the return
    % structure
    
    % Transformation microstructure
    AusRecon = myEBSD.Recon.Ebsd;
    
    % Reconstruction phase ID
    ReconID = myEBSD.Phase.ID{2};
    
    % Remove indexing for initially indexed austenite so only the
    % martensite points will actually be plotted
    OrigEbsd = myEBSD.Ebsd;
    InitAus = find(OrigEbsd.phase==ReconID);
    
    % ID each specific martensite boundary throughout the microstructure.
    % First iteration designates packet boundaries; second designates block
    % boundaries; third characterizes the sub-block boundaries.
    for ii = 1:3
        if ii == 1
            Bounds = myEBSD.Packets.IDs;
            len = 5;
        elseif ii == 2
            Bounds = myEBSD.Blocks.IDs;
            len = 13;
        else
            Bounds = myEBSD.Variants.IDs;
            len = 121;
        end
        
        % Remove the Ids of those corresponding to non-martensite in the
        % post-transformation microstructure (such as retained austenite)
        Bounds(InitAus)=[];
        for jj = 1:len
            CharBnd{jj} = find(Bounds == jj-1);
        end

        % Choose 5 random orientations
        RandOrs = orientation.rand(len,1,AusRecon('Austenite').CS,myEBSD.SS);
        tmpEb = AusRecon;

        % Assign a specific orientation to a packet ID
        for kk = 1:len
            tmp = CharBnd{kk};
            if isempty(tmp)==0
                tmpEb(tmp).orientations = RandOrs(kk);            
            end
        end
        
        % Compute grain structure
        [tmpGrains] = calcGrains(tmpEb,'angle',2*degree);
        
        if ii == 1
            myEBSD.Packets.Boundaries = tmpGrains;
        elseif ii == 2
            myEBSD.Blocks.Boundaries = tmpGrains;
        else
            myEBSD.Variants.SubBlockBoundaries = tmpGrains;
        end
    end
end

