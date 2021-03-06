function PltSubBlocks(varargin)
    % Plot variants within blocks in a systematic manner and see if we can
    % determine a rough estimate of the sub-block boundaries
    
    % Var ID
    myEBSD = varargin{1};
    
    if length(varargin)==2
        mbar = 'MicronBar';
        mkey = 'off';
    else
        mbar = [];
        mkey = [];
    end
    % Id the variants and the reconstruction phase
    Vars = myEBSD.Variants.IDs;
    ReconID = myEBSD.Phase.ID{2};
    
    % Remove non-martensite (retained austenite) in the post-transformation
    % microstructure
    InitAus = find(myEBSD.Ebsd.phase==ReconID);
    Vars(InitAus) = [];
    
    % Call the martensite and reconstruction microstructures
    % Phase IDs
    TransID = myEBSD.Phase.ID{1};
    ReconID = myEBSD.Phase.ID{2};
    
    % Post-transformation martensite structure and homogenize it
    Martensite = myEBSD.Ebsd(find(myEBSD.Ebsd.phase==TransID));
    RecEbsd = myEBSD.Recon.Ebsd;
    
    % Find all unique variants and reassign based on a 1 to 24 variant
    % labeling system
    
    ChngVars = Vars;
    for jj = 1:121
        if jj == 1
            SubBlck_Full{jj} = find(Vars==0);
        else
            SubBlck_Full{jj} = find(Vars == jj-1);
            altVar = (ceil((jj-1)/24)-1)*24;
            ChngVars(SubBlck_Full{jj}) = ChngVars(SubBlck_Full{jj})-altVar;
        end
    end 
    
    for kk = 1:25
        if kk == 1
            SubBlck{kk} = SubBlck_Full{kk};
        else
            SubBlck{kk} = find(ChngVars == kk-1);
        end
    end
    
%     % Classify 120 unique orientations
%     RandOrs = orientation.rand(120,1,RecEbsd.CS,myEBSD.SS);
%     tmpEb = RecEbsd;
%     
%     % Assign orientation to individual variants based on label
%     for kk = 1:120
%         tmp = SubBlck{kk};
%         if isempty(tmp)==0
%             tmpEb(tmp).orientations = RandOrs(kk);            
%         end
%     end
    
    % Extract sub-block RGB values
    SubBlckRGB = myEBSD.Variants.RGB;
    
    % Call orientation mapping factors
    key = ipfHSVKey(Martensite);
    key.inversePoleFigureDirection=zvector;
    Mor = Martensite.orientations;
    mp = key.orientation2color(Mor);
    
    figure; 
    plot(Martensite,mp,mbar,mkey);
    hold on
    % Colormap
    cmap = SubBlckRGB;
    caxis([0,25]);
    colormap(cmap);
    
    % Assign FaceColor to each respective sub-block
    for ii = 1:25
        tmp = SubBlck{ii};
        if isempty(tmp)==0
            plot(Martensite(tmp),mp(tmp),'FaceColor',cmap(ii,:));
            hold on
        end
    end 
%     cm = colorbar;
%     set(cm, 'Ticks', 0.5:24.5, 'TickLabels',[0:24]);
    
    hold on
    
%     % Overlay the block boundaries
%     BlckBounds = myEBSD.Blocks.Boundaries;
%     plot(BlckBounds.boundary,'FaceColor','black')
end



