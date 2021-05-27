function [Parent,Twin,TP_Count] = GrnMerging(myEBSD,Ebsd,HomGrns,TransGrns,Parent,Twin,TP_Count)
% The purpose of this function was to rectify any "discontinuous" grains
% that may have been assigned as separate grains in Recon_Mixed_SymOP but
% should actually be a continuation of a single grain. However, a different
% approach was taken that should absolve this issue, although I'm not
% entirely sure it solved the entire issue. Please contact me at:
% brust.15@buckeyemail.osu.edu if you run into any issues.

    % Find the unique orientations within the grain structure and make sure
    % we don't have disconnected parent-twin systems (twin grain actually
    % lies outside of parent grain [poor regularization], thus making it
    % it's own unique grain instead)
    [unqOrs,UnqOrInds] = unique(HomGrns.meanOrientation);
    GrnFlg = zeros(length(HomGrns),1);
    GrnFlg(UnqOrInds) = 1;
    GrnFlg(TransGrns) = 1;
    BadPTGrns = find(GrnFlg==0);

    % Get the grain ids for the full ebsd dataset, with previously
    % indexed austenite being labeled with a 0
    FullEbGrnIds = zeros(length(myEBSD.Recon.FullEbsd),1);
    FullEbGrnIds(Ebsd.id) = Ebsd.grainId;
    myEBSD.Recon.FullEbsd.grainId = FullEbGrnIds;
    FEb = myEBSD.Recon.FullEbsd;
    
    % Search through any disconnected parent-twin systems and reassign
    % these as individual parent grains instead of split parent-twin
    % systems
    if isempty(Twin.AllIndices)==0
        for kk = 1:length(BadPTGrns)
            RepGrns = find(HomGrns.meanOrientation==HomGrns(BadPTGrns(kk)).meanOrientation);
            findPT = 0;
            while findPT < length(Twin.Merged)
                findPT = findPT + 1;
                if HomGrns.meanOrientation(RepGrns(1))==Twin.Parent.Or{findPT}
                    if isempty(Parent.AllIndices) == 0
                        LP = length(Parent.Indices);
                    else
                        LP = 0;
                    end

                    % Determine whether there still exists a twin within the
                    % parent and, if not, reassign both as new parent grains
                    % then eliminate that twin index
                    if length(RepGrns)<2
                        G1Or = Twin.Parent.Or{findPT};
                        G2Or = Twin.Single.Or{findPT};
                        IndG1 = find(FEb.grainId==RepGrns(1));
                        IndG2 = find(FEb.grainId==RepGrns(2));
                        FullGInds = [IndG1;IndG2];
                        GIndsStrct{1} = IndG1;
                        GIndsStrct{2} = IndG2;

                        % Depending on the size of the variable orientation
                        % portions, assign the parent grain to the larger
                        % section (1) and the twin grain to the smaller
                        % section (2)
                        if FEb(IndG1(1)).orientations==G2Or
                            G2Ind = IndG2; 
                            IndG2 = IndG1;
                            IndG1 = G2Ind;
                            Rep2 = RepGrns(2);
                            RepGrns(2) = RepGrns(1);
                            RepGrns(1) = Rep2;
                        end

                        % Remove the elements of the twin structure that we
                        % just reassigned to be parent grains for the
                        % proper bookkeeping
                        Twin.Merged(findPT) = [];
                        Twin.Parent.Or(findPT) = [];
                        Twin.Single.Or(findPT,:) = [];
                        Twin.Single.ID(findPT,:) = [];
                        TwnTmpInds = zeros(length(FEb),1);
                        TwnTmpInds(Twin.AllIndices) = 1;
                        TwnTmpInds(FullGInds) = 0;
                        Twin.AllIndices = find(TwnTmpInds==1);
                        TP_Ors(findPT) = [];

                        % Reduce the parent-twin counter
                        TP_Count = TP_Count-1;  

                        % Ensure that, if any twins exist within a grain
                        % we're newly indexing, it's not just noise and
                        % actually consists of a tangible grain (i.e. > 25
                        % indexed pixels)
                        for jk = 1:2
                            tmpGInds = GIndsStrct{jk};
                            Gid = RepGrns(jk);
                            tmpGrn_Ors = unique(FEb(tmpGInds).orientations);
                            LenGOrs = length(tmpGrn_Ors);
                            if LenGOrs > 1
                                % Determine the sizes of the parent and twinned
                                % portions
                                GOrInds = [];
                                GOrsLen = 0;
                                for ij = 1:LenGOrs
                                    GOrInds{ij} =...
                                        find(FEb(tmpGInds).orientations==tmpGrn_Ors(ij));
                                    GOrsLen(ij) = length(GOrInds{ij});
                                end
                                % Flag these regions of varying orientations to
                                % ensure that they meet the size requirement (>25
                                % pixels)
                                GOrFlg = GOrsLen < 25;

                                % Identify the largest indexing (parent or
                                % grain) orientation
                                [~,maxG] = max(GOrsLen);
                                maxOr = tmpGrn_Ors(maxG);

                                % Remove any noise that may exist to ensure a
                                % correct indexing of the parent grain or
                                % parent-twin system
                                for ij = 1:LenGOrs
                                    if GOrFlg(ij)
                                        tmpG =...
                                            find(FEb(tmpGInds).orientations==tmpGrn_Ors(ij));
                                        FEb(tmpGInds(tmpG)).orientations = maxOr;
                                        LenGOrs = LenGOrs-1;
                                    end
                                end
                                % Remove this ``twin'' altogether
                                tmpGrn_Ors(GOrFlg)=[];
                                GOrInds(GOrFlg)=[];
                                GOrsLen(GOrFlg)=[];
                                % Reidentify the largest grain
                                [~,maxG] = max(GOrsLen);
                            end

                            % Search through first set of indices to see if
                            % there exists a twinned grain within this region
                            % or if it's a true parent grain now
                            if LenGOrs == 1
                                Parent.Indices{LP+1} = tmpGInds;
                                Parent.Or{LP+1} = tmpGrn_Ors;
                                PInds = Parent.AllIndices;
                                PInds = sort(vertcat(PInds,tmpGInds));
                                Parent.AllIndices = PInds;
                                LP = LP+1;
                            else
                                % If two or more orientations exist, reclassify
                                % this portion as a twinned grain instead and fill
                                % in the Twin structure accordingly

                                % Add index to Twin structure
                                TP_Count = TP_Count+1;
                                % Add the merged indices indicating the parent-twin
                                % system
                                Twin.Merged{TP_Count} = tmpGInds;
                                % Based on whichever orientation index has more
                                % indexed points, classify the parent orientation
                                % as the largest ``portion'' then delete from the
                                % G1_Ors array
                                Twin.Parent.Or{TP_Count} = maxOr;
                                HomGrns(Gid).meanOrientation = maxOr;
                                TP_Ors(TP_Count) = maxOr;

                                tmpGrn_Ors(maxG) = [];
                                LenGOrs = LenGOrs-1;

                                % Now rotate this maximum orientation by the number
                                % of twins existing in the system to identify which
                                % twin(s) are present
                                GPar_Or = rotation('euler',maxOr.phi1,maxOr.Phi,maxOr.phi2,CS_R);
                                TwnGOrs = transpose(GPar_Or * TwnG_eul);

                                % Now identify the existing twins accordingly and
                                % add to the Twin Structure
                                TG_Count = 1;
                                for ij = 1:LenGOrs
                                    TG_Flg = TwnGOrs==tmpGrn_Ors(ij);
                                    TG_ID = find(TG_Flg);
                                    Twin.Single.Or{TP_Count,TG_Count} = TwnGOrs(TG_ID);
                                    Twin.Single.ID{TP_Count,TG_Count} = TG_ID;
                                    TG_Count = TG_Count+1;
                                end
                                TInds = Twin.AllIndices;
                                TInds = sort(vertcat(TInds,tmpGInds));
                                Twin.AllIndices = TInds;
                            end
                        end

                        % Since we found the faulty Parent-Twin system, set
                        % our findPT counter to high number to exit loop
                        findPT = 1e6;
    %                         
                    end % Repeated grains loop
                end % Loop through parent-twin orientations
            end % While looping through twins condition
        end % If bad parent-twin grain systems are present (kk)
    end % If twins exist in our reconstructed microstructure
    
end

