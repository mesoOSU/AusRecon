function [Grains,Twin,Parent,myEBSD] = GrainSize_MixedSpace(varargin)
%%                      Function Description

% This function takes the output from the reconstruction functions and
% determines the grain boundary structure and corresponding bookkeeping for
% all of the parent phase grains. 

% The function utilizes the MTEX grain calculation function, which uses 
% voronoi tesselation with misorientation tolerances to define the grain 
% boundaries. Some modifications are then implemented for proper
% bookkeeping that allow for ensuing functions to work off of the results
% presented here.

% The function follows ASTM standards by combining any annealling twins 
% found in the reconstructed parent grain
% with the parent grain such that the parent-twin system becomes a single
% entity whose grain boundaries are only defined by the parent grain.

% Input Variables: myEBSD structure; Twin grain structure; Parent grain
%                  structure; VARIABLE: Grain size pixel tolerance

% Output Variables: Grain structure containing pertinent information;
%                   modified myEBSD, Twin, and Parent structures

%%                      Function Implementation

    % Set the constant varargin values
    myEBSD = varargin{1};
    Twin = varargin{2};
    Parent = varargin{3};
    
    % If a pixel tolerance is defined by the user, respect it. If not, the
    % default is ASTM standard 100 pixel grain size tolerance
    if length(varargin) == 4
        PixTol = varargin{4};
    else
        PixTol = 100;
    end
    
    % Preallocate global variables
    Ebsd = myEBSD.Recon.FullEbsd;       % Reconstructed (only) parent phase EBSD
    ChldID = myEBSD.Phase.ID{1};        % Child phase ID
    ParID = myEBSD.Phase.ID{2};         % Parent phase ID
    CS_T = myEBSD.CS{2};                % Transformed Phase
    CS_R = myEBSD.CS{3};                % Reconstructed Phase
    
    % Account for all points phases
    PhaseIDs = myEBSD.Recon.PhaseIDs;
    
    % If there exists unindexed points, assign a phase value of -1 to
    % ensure we flag them correctly
    PhaseIDs(Ebsd.isIndexed==0)=-1;

    % If assuming each twin creates a new grain or no twins exist in
    % reconstructed microstructure and save the parent-twin orientations
    TP_Ors = [];
    TP_Count = 0;
    if isempty(Twin.AllIndices)
        Merged_Ors = orientation('Euler',[0,0,0]);
        
    % Otherwise, merge the parent and twins into single grain
    else
        Combined = Twin.Merged;
        TwinPar = Twin.Parent.Or;
        for i = 1:length(Combined)
            Ebsd(Combined{i}).orientations = TwinPar{i};
            Merged_Ors(i) = TwinPar{i};
            TP_Ors = vertcat(TP_Ors,TwinPar{i});
        end
        TP_Count = i;
        Merged_Ors = unique(Merged_Ors);
    end

    % Identify the Sigma3 parent and twin systems
    Twn = vector3d([1,1,1;-1,-1,1;-1,1,1;1,-1,1]');
    TwnG_eul = orientation('axis',Twn,'angle',60*degree,CS_R);

    % Find original austenite points and any points that were misindexed 
    % and delete from the grain structure for cleaner images
    OrigAus = find(myEBSD.Ebsd.phase==ParID);
    UnIndPts = find(Ebsd.isIndexed==0);

    % Find the ebsd indexed points that were not assigned an austenite
    % orientation
    EbisInd = Ebsd(Ebsd.isIndexed);
    TransPts = find(PhaseIDs==ChldID);

    % Now assign a non-misorientation (euler=[0,0,0]) to all retained
    % transformation phase indices for proper bookeeping.
    if length(unique(EbisInd.phase)) > 1
        % Cluster martensite grains as non-oriented martensite points
        o0 = orientation('Euler',[0,0,0],CS_T);
    else
        % Due to MTEX's multiple phase limitations, assign the
        % untransformed martensite points a non-orientation while retaining
        % the transformation phase crystal symmetry
        o0 = orientation('Euler',[0,0,0],CS_R);
    end
    Ebsd(TransPts).orientations = o0;

    % Delete originally indexed parent phase and any unindexed points
    DelPts = [OrigAus;UnIndPts];
    Ebsd(DelPts)=[];
    
    % Due to the high texture with titanium, ensure a lower misorientation
    % tolerance than steel
    if strcmp(myEBSD.Material,'Titanium')
        MisDeg = 2 * degree;
    else
        MisDeg = 5 * degree;
    end

    % Calculate grain boundary structure using MTEX function
    [grains,Ebsd.grainId] = calcGrains(Ebsd,'angle',MisDeg);

    % Find any grains that have less than user-defined number of pixels. 
    % These grains are essentially ``noise'', and are removed from the 
    % grain boundaries assigned to the corresponding parent microstructure[
    if PixTol > 0
        [Ebsd,grains,SmllGrns] = RemSmGrns(Ebsd,grains,PixTol);
    end

    % Homogenize grain phases to bypass MTEX restrictions
    HomGrns = grains;
    if length(unique(HomGrns.phase))>1
        HomGrns.phase=ParID;
        TransGrns = find(HomGrns.meanOrientation==o0);
    else
        % All transformation phase indices have been assigned to the parent
        % phase, thus no faux' transformation grains exist and we assign it
        % an empty vector.
        TransGrns=[];
    end
    
    % Try to rectify any grains assigned through poor regularization
    % parameterization
    % NOTE: I'm pretty sure I fixed this issue in the Recon_Mixed_SymOP
    % function. It should enter the function, realize it's not needed, and
    % then exit. Problem is I'm not entirely positive this is true, so a
    % breakpoint is inserted in here if the function continues. I didn't
    % debug because I couldn't catch it, so if a problem is encountered
    % here, please contact me at:
    %                   brust.15@buckeyemail.osu.edu
    [Parent,Twin,TP_Count] =...
        GrnMerging(myEBSD,Ebsd,HomGrns,TransGrns,Parent,Twin,TP_Count);

    % Add grain size, grain area, and grain id to returned Grain function
    Grn_sz = grains.grainSize;
    Grn_area = area(grains);
    Grn_ID = grains.id;

    % Number of grains
    num_grains = length(Grn_ID);
    twn_count = 0;
    untrans_count = 0;
    grn_count = 1;
    EBSDs = cell(num_grains,1);

    % Save either full parents or combined twins 
    for j = 1:num_grains
        tmpId = Grn_ID(j);
        EbId = find(Ebsd.grainId==tmpId)';
        EBSDs{grn_count} = EbId;
        grn_count = grn_count+1; 
    end

    % Create bounding box around EBSD data set
    BBInds = find(Ebsd.x == min(Ebsd.x));
    BBInds = vertcat(BBInds,find(Ebsd.x == max(Ebsd.x)));
    BBInds = vertcat(BBInds,find(Ebsd.y == min(Ebsd.y)));
    BBInds = vertcat(BBInds,find(Ebsd.y == max(Ebsd.y)));
    BBInds = unique(BBInds);
    % Now extract the unique grain IDs bordering the IPF map and delete
    % from the area calculations
    EbsdBB = unique(Ebsd(BBInds).grainId);
    for jj = 1:length(SmllGrns)
        if any(EbsdBB == SmllGrns(jj))
            EbsdBB(find(EbsdBB == SmllGrns(jj))) = [];
        end
    end

    % Ebsd step size
    stepsz = Ebsd.scanUnit;
    if stepsz == 'um'
        delta = 1e-6;
        scale = 1e6;
    elseif stepsz == 'mm'
        delta = 1e-3;
        scale = 1;
    else
        warning('Ebsd step size not recorded')
    end

    % Multiply MTEX measured grain area by step size squared
    ASTMGrnArea = Grn_area * delta^2;

    % If hex-grid, modify to follow ASTM E2627 standards
     if length(Ebsd.unitCell) == 6
         ASTMGrnArea = ASTMGrnArea * sqrt(3)/2;
     end

    % Remove edge-bounded grains from any grain size calculations
    InternalGrnIds = Grn_ID;
    InternalGrnIds(EbsdBB) = [];
    ExternalGrnIds = Grn_ID(EbsdBB);

     % Compute ASTM number (What ASTM is this?)
    %      ASTM_num = 16.97737-3.321928*log10(mean(ASTMGrnArea))*1.0;    
     % Compute ASTM E2627 grain size number
    ASTM_num = -3.3223*log10(mean(ASTMGrnArea(InternalGrnIds)*scale))-2.995;
    Sz = size(grains(InternalGrnIds),1);
    freqGS = 1/Sz * sum(ASTMGrnArea(InternalGrnIds));
    
    StDev = sqrt(1/(Sz-1) * sum((ASTMGrnArea(InternalGrnIds) - freqGS).^2));
    StError = StDev / sqrt(Sz);
    CI95 = tinv([0.025,0.975],Sz-1);
    yCi95 = bsxfun(@times, StError, CI95(:));

    % Now measure grain sizes using the lineal intercept method
    [myEBSD,GrnInt] = BoundaryIntersect(myEBSD,grains);

    % Fill our grain structure with pertinent variables related to both the
    % assigned grain structure and the ASTM E2627 data
    Grains.grains = grains;
    Grains.Area = Grn_area;
    Grains.Size = Grn_sz;
    Grains.Indices = EBSDs;
    Grains.ASTM.E2627.G = ASTM_num;
    Grains.ASTM.E2627.ExcludedGrns = ExternalGrnIds;
    Grains.ASTM.E2627.IncludedGrns = InternalGrnIds;
    Grains.ASTM.E2627.Ci95 = yCi95;

    % Now set up the grain id so we can determine whether an indexed grain
    % is a pure parent or if it 
    GrnIds = zeros(num_grains,7)-1;
    GrnIds(:,1) = Grn_ID;
    GrnIds(:,2) = ChldID;

    % Fill in untransformed ebsd dataset 
    if untrans_count > 0
        Grains.Untransformed.Ebsd = UnTrans;
    end

    % Separate indexing for merged parents and twins in same system
    if twn_count>0
        Grains.Ebsd.Merged.Recon = Merge_grn;
        Grains.Ebsd.Merged.Trans = Merge_trans;
        Grains.Ebsd.Merged.ID = Merge_Id;
    end
    
    % Establish the maximum phase number for the number of phases in the
    % data set
    PhIds = [ChldID; ParID];
    mxPh = max(PhIds);

    % Only perform the array characterization if twins exists in the
    % microstructure
    if TP_Count > 0

        % How many twin types (rotations) exist
        TwnSz = size(Twin.Single.ID);
        TwnType = TwnSz(2); 

        % Identify which specific twin (based on rotation) is present in 
        % the merged parent-twin cases
        for jj = 1:length(grains)
            if jj == 68
                keyboard
            end
            
            % If the grain is oriented normal to the reference frame, it
            % is an unassigned martensite region and we ignore it
            if HomGrns.meanOrientation(jj) ~= o0
                % Loop through the twin-parent orientations to find the
                % match
                twnflg = 0;
                for kk = 1:TP_Count
                    % If match is found, ID the twin rotations (1 -> 4)
                    if HomGrns.meanOrientation(jj) == TP_Ors(kk)
                        
                        tmp = [];
                        for mm = 1:TwnType
                            tmp = vertcat(tmp,Twin.Single.ID{kk,mm});
                        end
                        GrnIds(jj,tmp+2) = tmp;
                        GrnIds(jj,2) = mxPh+1;
                        GrnIds(jj,7) = kk;
                        twnflg = 1;
                    end
                end
                if twnflg == 0
                    GrnIds(jj,7) = 0;
                    GrnIds(jj,2) = ParID;
                end
            end
        end
    else
        if length(unique(Ebsd.phase))==1
            ParGrns = find(grains.meanOrientation~=o0);
            GrnIds(ParGrns,2) = ParID;
        end
    end

    % Send the parent-twin id back with the Grain structure
    Grains.grainId = GrnIds;
    Grains.EbsdId = Ebsd.grainId;
    Grains.ASTM.E112_13 = GrnInt;
    % Add grain Ids to the Full, reconstructed Ebsd data set 
    myEBSD.Recon.MergedEbsd = Ebsd;
    
end

