function [Grains,Twin,Parent] = GrainSize(myEBSD,Twin,Parent,Merge)

    Ebsd = myEBSD.Recon.FullEbsd;
    FullEbsd = myEBSD.Recon.FullEbsd;
    PhaseIDs = myEBSD.Recon.PhaseIDs;
    TransID = myEBSD.Phase.ID{1};
    ReconID = myEBSD.Phase.ID{2};
    CS_T = myEBSD.CS{2};                        % Transformed Phase
    CS_R = myEBSD.CS{3};                        % Reconstructed Phase
    TransEbsd = myEBSD.TransEbsd;
    ReconPts = myEBSD.Recon.TransformedPoints;
    TransPts = find(PhaseIDs==TransID);
    
    % Cluster martensite grains as non-oriented martensite points
    o0 = orientation('Euler',[0,0,0],CS_T);
    Ebsd(TransPts).orientations = o0;
    % If assuming each twin creates a new grain or no twins exist in
    % reconstructed microstructure and save the parent-twin orientations
    TP_Ors = [];
    if isempty(Twin.AllIndices) || Merge == 0
        Merged_Ors = orientation('Euler',[0,0,0]);
        
    % If ignoring twin boundaries, merge twins and parent grains together
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
    
    % Find original austenite points and delete from the grain structure
    % for cleaner images
    OrigAus = find(myEBSD.Ebsd.phase==ReconID);
    Ebsd(OrigAus)=[];
    Ebsd('Martensite') = [];
    [grains,Ebsd.grainId] = calcGrains(Ebsd,'angle',2*degree);
    grains(grains.grainSize<25)=[];
    % Find the unique orientations within the grain structure and make sure
    % we don't have disconnected parent-twin systems (twin grain actually
    % lies outside of parent grain [poor regularization], thus making it
    % it's own unique grain instead)
    [unqOrs,UnqOrInds] = unique(grains.meanOrientation);
    GrnFlg = zeros(length(grains),1);
    GrnFlg(UnqOrInds) = 1;
    BadPTGrns = find(GrnFlg==0);
    
    % et the grain ids for the full ebsd dataset, with previously
    % indexed austenite being labeled with a 0
    FullEbGrnIds = zeros(length(FullEbsd),1);
    FullEbGrnIds(Ebsd.id) = Ebsd.grainId;
    myEBSD.Recon.FullEbsd.grainId = FullEbGrnIds;
    FEb = myEBSD.Recon.FullEbsd;
    
    % Search through any disconnected parent-twin systems and reassign
    % these as individual parent grains instead of split parent-twin
    % systems
    szTwn = size(Twin.Single.Or);
    Tflg = 0;
    for kk = 1:length(BadPTGrns)
        RepGrns = find(grains.meanOrientation==grains(BadPTGrns(kk)).meanOrientation);
        findPT = 0;
        while findPT < length(Twin.Merged)
            findPT = findPT + 1;
            if grains.meanOrientation(RepGrns(1))==Twin.Parent.Or{findPT}
                LP = length(Parent.Indices);
                
                % Determine whether there still exists a twin within the
                % parent and, if not, reassign both as new parent grains
                % then eliminate that twin index
                if length(RepGrns)>2
                    % Remove faux twin and make it it's own grain, but 
                    % maintain the parent-twin system
%                     for mm = 1:szTwn(2)
%                         if isempty(Twin.Single.Or{findPT})==0
%                             TOr = Twin.Single.Or{findPT};
%                             if TOr == 
                else
                    POr = Twin.Parent.Or{findPT};
                    TOr = Twin.Single.Or{findPT};
                    IndP = find(FEb.grainId==RepGrns(1));
                    IndT = find(FEb.grainId==RepGrns(2));
                    % Set Parent Grain
                    if FEb(IndP(1)).orientations==TOr
                        tInd = IndT; 
                        IndT = IndP;
                        IndP = tInd;
                    end
                    % Now add the ''grain'' indices and orientations to the
                    % parent structure
                        Parent.Indices{LP+1} = IndP;
                        Parent.Or{LP+1} = POr;
                        Parent.Indices{LP+2} = IndT;
                        Parent.Or{LP+2} = TOr;
                        Inds = [IndP;IndT];
                        PInds = Parent.AllIndices;
                        PInds = sort(vertcat(PInds,Inds));
                        Parent.AllIndices = PInds;
                        
                        % Remove the elements of the twin structure that we
                        % just reassigned to be parent grains for the
                        % proper bookkeeping
                        Twin.Merged(findPT) = [];
                        Twin.Parent.Or(findPT) = [];
                        Twin.Single.Or(findPT,:) = [];
                        Twin.Single.ID(findPT,:) = [];
                        TwnTmpInds = zeros(length(FEb),1);
                        TwnTmpInds(Twin.AllIndices) = 1;
                        TwnTmpInds(Inds) = 0;
                        Twin.AllIndices = find(TwnTmpInds==1);
                        TP_Ors(findPT) = [];
                        
                        % Reduce the parent-twin counter
                        TP_Count = TP_Count-1;     
                end
            end
        end
    end
    % Add grain size, grain area, and grain id to returned Grain function
    Grn_sz = grains.grainSize;
    Grn_area = area(grains);
    Grn_ID = grains.id;
    
    % Number of grains
    num_grains = length(Grn_ID);
    par_count = 1;
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
    
    % Compute ASTM number
    ASTM_num = 16.97737-3.321928*log10(mean(Grn_area))*1.0;
    
    % Fill our grain structure with pertinent variables
    Grains.grains = grains;
    Grains.Area = Grn_area;
    Grains.Size = Grn_sz;
    Grains.Indices = EBSDs;
    Grains.ASTMNum = ASTM_num;
    
    % Homogenize grain phases for mtex purposes
    HomGrns = grains;
    HomGrns.phase=ReconID;
    
    % Now set up the grain id so we can determine whether an indexed grain
    % is a pure parent or if it 
    GrnIds = zeros(num_grains,7)-1;
    GrnIds(:,1) = Grn_ID;
    GrnIds(:,2) = TransID;
    
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
  
    PhIds = [TransID; ReconID];
    mxPh = max(PhIds);
    
    % Only perform the array characterization if twins exist in the
    % microstructure
    if TP_Count > 0
    
        % How many twin types (rotations) exist
        TwnSz = size(Twin.Single.ID);
        TwnType = TwnSz(2); 

        % Identify which specific twin (based on rotation) is present in the 
        % merged parent-twin cases
        for jj = 1:length(grains)
            % If the grain is orientated normal to the reference frame, it
            % is an unassigned martensite region and we ignore
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
                    GrnIds(jj,2) = ReconID;
                end
            end
        end
    end

    % Send the parent-twin id back with the Grain structure
    Grains.grainId = GrnIds;
    Grains.EbsdId = Ebsd.grainId;
                
    
    
end

