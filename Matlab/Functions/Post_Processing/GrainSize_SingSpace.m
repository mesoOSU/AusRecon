function [Grains] = GrainSize(myEBSD,Twin,Merge)

    Ebsd = myEBSD.Recon.Ebsd;
    recon_pts = myEBSD.Recon.TransformedPoints;
    Ebsd = Ebsd(recon_pts);
    EbsdOrig = Ebsd;
    TransEbsd = myEBSD.TransEbsd(recon_pts);

    % If assuming each twin creates a new grain or no twins exist in
    % reconstructed microstructure
    if isempty(Twin) || Merge == 0
        Merged_Ors = orientation('Euler',[0,0,0]);
        
    % If ignoring twin boundaries, merge twins and parent grains together
    else
        Combined = Twin.Merged;
        TwinPar = Twin.Parent.Or;
        for i = 1:length(Combined)
            Ebsd(Combined{i}).orientations = TwinPar{i};
            Merged_Ors(i) = TwinPar{i};
        end

        Merged_Ors = unique(Merged_Ors);
    end
    
    grains = calcGrains(Ebsd,'angle',2*degree);
    grains(grains.grainSize<25) = [];
    Grn_sz = grains.grainSize;
    Grn_area = area(grains);
    Grn_ID = grains.id;

    % Number of grains
    num_grains = length(Grn_ID);
    par_count = 0;
    twn_count = 0;
    % Save either full parents or combined twins 

    for j = 1:num_grains
        tmpId = Grn_ID(j);
        Grn_Or = grains(j).meanOrientation;
        EbId = find(Ebsd.orientations==Grn_Or);
        Merge_flag = any(Merged_Ors==Grn_Or);
        if Merge_flag == 0
            par_count = par_count+1;
            Par_grn{par_count,1} = Ebsd(EbId);
            Par_trans{par_count,1} = TransEbsd(EbId);
            Par_Id(par_count,1) = tmpId;
        else
            twn_count = twn_count+1;
            Merge_grn{twn_count,1} = EbsdOrig(EbId);
            Merge_trans{twn_count,1} = TransEbsd(EbId);
            Merge_Id(twn_count,1) = tmpId;
        end   
    end       
        
    Grains.grains = grains;
    Grains.Area = Grn_area;
    Grains.Size = Grn_sz;
    Grains.Ebsd.Parent.Recon = Par_grn;
    Grains.Ebsd.Parent.Trans = Par_trans;
    Grains.Ebsd.Parent.ID = Par_Id;
    
    % If there are twins
        if isempty(Twin)==0
        %     % Merged Parent-Twin grain average size
        %     Merged_Grainsize = mean(grains(Merged_ID).grainSize);
            % Don't think this is calculated correctly
            ASTM_num = 16.97737-3.321928*log10(mean(Grn_area))*1.0;
%             ASTM_num=[];

            Grains.Ebsd.Merged.Recon = Merge_grn;
            Grains.Ebsd.Merged.Trans = Merge_trans;
            Grains.Ebsd.Merged.ID = Merge_Id;
            Grains.ASTMNum = ASTM_num;
        end
end

