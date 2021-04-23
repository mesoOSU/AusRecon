function [Grains] = GrainSize(myEBSD,Twin,Merge)

    Ebsd = myEBSD.Recon.Ebsd;
    PhaseIDs = myEBSD.Recon.PhaseIDs;
    TransID = myEBSD.Phase.ID{1};
    ReconID = myEBSD.Phase.ID{2};
    CS_T = myEBSD.CS{2};                        % Transformed Phase
    CS_R = myEBSD.CS{3};                        % Reconstructed Phase
    EbsdOrig = Ebsd;
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
    [grains,Ebsd.grainId] = calcGrains(Ebsd,'angle',2*degree);
%     grains(grains.grainSize<25)=[];
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
%        EbId = find(Ebsd.grainId==tmpId)';
% Quickfix: should return single id
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

