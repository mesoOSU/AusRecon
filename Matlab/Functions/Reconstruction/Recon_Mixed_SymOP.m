function [myEBSD,FullParent,FullTwin] = Recon_Mixed_SymOP(myEBSD,num_iters,OP,IP)
%%
    % Initialization of relevant parameters independent of quadrant setup
    CS_T = myEBSD.CS{2};                    % Transformed Phase
    CS_R = myEBSD.CS{3};                    % Reconstructed Phase
    SS = myEBSD.SS;                         % Sample Symmetry
    T2R = myEBSD.TransMats{1};              % Transformation to Reconstruction phase
    R2T = myEBSD.TransMats{2};              % Transformation from Reconstruction Phase
    Trans_MODF = myEBSD.MODF;               % Misorientation Distribution Function
    TransID = myEBSD.Phase.ID{1};           % Phase ID for Transformation phase
    ReconID = myEBSD.Phase.ID{2};           % Phase ID for Reconstruction phase 
    psi = myEBSD.noise.psi;                 % Psi value 
    n_phases = myEBSD.Phase.n_phases;       % Number of phases
    FullEbsd = myEBSD.Ebsd;                 % Full martensite microstructure
    FullWts = zeros(length(FullEbsd),1);    % Zeroed vector of weighted values
    ci = myEBSD.ci;                         % Confidence interval for each orientation
    
    % Vector over entire dataset exhibiting the phase IDs
    FullPhIds = FullEbsd.phase; 

    % Number of quadrants
    nQds = length(myEBSD.TransEbsd);
    
    TotCount = 0;
    if myEBSD.Quad.Flag
        ConsNoCuts = 10;
    else
        num_iters = round(num_iters/nQds);
        ConsNoCuts = 7;
    end
    
    % Variable to be used later
     UnqOrs = [];
   
    for ii = 1:nQds
        % Initialization of parameters dependent on quadrant setup
        Trans_ebsd = myEBSD.TransEbsd{ii};              % Transformation-phase ebsd
        TransEbsd_ci = ci(Trans_ebsd.id);               % Confidence interval for specific ebsd
        IP_wts = myEBSD.Weights.inplane{ii};            % In-plane weights
        c2g_wts = myEBSD.Weights.outplane{ii};          % Initial c2g weights
        adjpts = myEBSD.Adj_pts{ii};                    % Adjacent Points Array
        Rec_ebsd = Trans_ebsd;                          % Reconstructed ebsd dataset

        % Establish the in-plane weights
        IP_wts = (IP_wts+IP(1)).*IP(2);

        % Add source node
        DiGraph=digraph;

        % Call function to set up the initial graph
        [DiGraph,Endnode,Sinknode] = graph_setup(DiGraph,adjpts,IP_wts,c2g_wts,1);

        % Index original EBSD dataset and the active dataset, which will delete
        % indices of transformed points
        EbsdInd = [1:length(Rec_ebsd)]';
        ActiveTransInd = EbsdInd;

        % Find untransformed points
        PhaseIDs = ones(length(Rec_ebsd),1).*TransID;
        PostTransID = find(PhaseIDs==TransID);
        Parent = [];
        Twin = [];

        % Initialize the active ebsd dataset and post-transformation 
        % orientation list
        ActTrans_ebsd = Trans_ebsd(PostTransID);
        ActTrans_Ors = Trans_ebsd(PostTransID).orientations;

        % Define active directed graph that changes 
        ActDiGraph = DiGraph;
        ActEndnode = Endnode;
        ActSinknode = Sinknode;

        % Parent, Twin, and Initial No-cut counters
        ParCount = 0;
        TwnCount = 0;
        NoCutCount=1;
        
        % Start the iteration counter
        iters = 0;

        while iters < num_iters

            % Increase iteration count
            iters = iters+1;
            if isempty(UnqOrs)
                SkipMart = 0;
                % Choose point from within untransformed microstructure based on the
                % confidence interval
                ActTrans_ci = TransEbsd_ci(PostTransID).*100;
                ci_id = [1:length(ActTrans_ci)]';
                id_pt = randsample(ci_id,1,true,ActTrans_ci);
                Orien = ActTrans_Ors(id_pt);

                % Remove OP edges from previous weights
                ActDiGraph=rmedge(ActDiGraph,(1:length(ActTrans_Ors))+ActEndnode,ActSinknode);

                % Compute the misorientations and corresponding OP weights based on the 
                % chosen point
                Misos = inv(Orien).*ActTrans_Ors;
%Austin
                %g2t_wts = (eval(Trans_MODF,Misos)+OP(1)).*OP(2);
                g2t_wts = (alt_eval(Trans_MODF,Misos)+OP(1)).*OP(2);

                % Eliminate any possible negative weights
                g2t_wts(find(g2t_wts<0))=0;

                % Add edges corresponding to the guess
                ActDiGraph=addedge(ActDiGraph,(1:length(ActTrans_Ors))+ActEndnode,ActSinknode,g2t_wts);

                % Perform graph cut 
                [mf_c,gf,cs,ct]=maxflow(ActDiGraph,1,ActSinknode);

                % Copy the cut for convenience
                ctcopy = ct;
                ctcopy(end)=[];

                % Index and length of index for first cut
                Guess_cutEBSD = ActiveTransInd(ctcopy-ActEndnode);
                Guess_cutActive = ctcopy-ActEndnode;
            else
                Guess_cutEBSD = ActiveTransInd;
                Guess_cutActive = [1:length(ActiveTransInd)]';
                SkipMart = 1;
            end
            len_preCT  = length(Guess_cutEBSD);

            % Preallocate arrays for twinning features
            CutInds = [];
            TwnLocInds= [];

            % If the first cut is empty, add to the NoCut Counter and continue
            % on. After a set amount of consecutive non-cuts (try 5 for now), 
            % this will force the while loop to stop and exit
            if len_preCT < 50

            NoCutCount = NoCutCount+1;
            % If the cut isn't empty, perform a second cut to refine the grain
            % boundaries
            else

    %             if len_preCT < 500
    %                 aus_guesses = symmetrise(ActTrans_Ors(Guess_cutActive),CS_R)*T2R;
    %                 Guess_Or = global_pole_figure_estimation(aus_guesses,CS_R,SS,1);
    %                 Guess_Or = Guess_Or(1);
    %             else
                    if SkipMart
                        Guess_Or = UnqOrs(1);
                    else
                        % Ignore orientations coming from (very) poorly indexed
                        % points in our guess computation
                        tmp_ci = find(ActTrans_ci(Guess_cutActive) < .05);
                        ActOrs = ActTrans_Ors(Guess_cutActive);
                        ActOrs(tmp_ci) = [];
                        if isempty(ActOrs)
                            ActOrs = ActTrans_Ors(Guess_cutActive);
                        end
                        % Compute odf and determine max recon orientation
                        local_odf=calcODF(symmetrise(ActOrs)*T2R,'kernel',psi);
                        [~,Guess_Or] = max(local_odf);
                    end
    %             end

                % Now create second graph based solely on the first cut
                New_Digraph = digraph;

                % Establish new EBSD dataset based on the initial guess cut
                TransCT_ebsd = ActTrans_ebsd(Guess_cutActive);

                % Call function to compute adjacency array and corresponding
                % misorientations for the neighboring points
                [adjpts,mori,~] = adjpt_moris(TransCT_ebsd,1);

                % Evaluated misorientation distribution function values (serving as
                % in-plane weights)
%Austin                 modf_vals=eval(myEBSD.MODF,mori);
                modf_vals=alt_eval(myEBSD.MODF,mori);
                modf_vals(modf_vals<0)=0;
                IPwts = (modf_vals+IP(1)).*IP(2);

                % Redefine the parent orientation as a rotation
                Par_Or = rotation('euler',Guess_Or.phi1,Guess_Or.Phi,Guess_Or.phi2,CS_R);

                % Check for annealing twins based on parent orientation if 
                % reconstructing austenite
                if strcmp(myEBSD.Material,'Steel')
                    Twns = vector3d([1,1,1;-1,-1,1;-1,1,1;1,-1,1]');
                    Twn_eul = orientation('axis',Twns,'angle',60*degree,CS_R);
                    Twn_guesses = transpose(Par_Or * Twn_eul);
                    Sym_Ors = [Guess_Or; Twn_guesses];
                elseif strcmp(myEBSD.Material,'Titanium')
                    % Don't search for annealing twins
                    Sym_Ors = Guess_Or;
                end   

                % Compute new g2t_wts that encompasses both the parent and
                % twin orientations
                From_Or_guess=symmetrise(Sym_Ors)*R2T;
                From_guess_ODF=calcODF(From_Or_guess,'kernel',psi);
                OPWts=(eval(From_guess_ODF,ActTrans_Ors(Guess_cutActive))+OP(1)).*OP(2);

                % c2g weights become the inverse of the g2t weights in order to
                % ensure a precise cut
                c2gWts = 4./OPWts;

                % Call function to set up the initial graph
                [New_Digraph,endnode2,sinknode2] = graph_setup(New_Digraph,adjpts,IPwts,c2gWts,2);

                % Add edges corresponding to the guess reconstruction orientation    
                New_Digraph=rmedge(New_Digraph,(1:len_preCT)+endnode2,sinknode2);
                New_Digraph=addedge(New_Digraph,(1:len_preCT)+endnode2,sinknode2,OPWts);

                % Perform graph cut 
                [mf_c2,gf2,cs2,ct2]=maxflow(New_Digraph,1,sinknode2);

                % Copy the cut for convenience
                ctcopy2 = ct2;
                ctcopy2(end)=[];

                % Index of newly cut points
                NewCT_ind = ctcopy2-endnode2;
                if length(NewCT_ind) < 25
                    NoCutCount = NoCutCount+1;
                else

                    % Now re-index for the entire ebsd scope and collect the uncut
                    % points as well
                    Curr_cutEBSD = Guess_cutEBSD(NewCT_ind);
                    Curr_cutActive = Guess_cutActive(NewCT_ind);
                    SymInds = Curr_cutEBSD;
                    SymIndsFlg = zeros(length(SymInds),1);
                    c2g_wts(Curr_cutEBSD) = OPWts(NewCT_ind);

                    % Update the transformed points indexing
                    PhaseIDs(Curr_cutEBSD) = ReconID;
                    len_Cpts = length(Curr_cutEBSD);

                    % Add cut to reconstructed portion of ebsd
                    Rec_ebsd(Curr_cutEBSD).orientations = Guess_Or;

                    % If titanium, don't search through for twins. The entire
                    % system will be the parent grain, as bcc doesn't typically
                    % produce annealing twins.
                    if length(Sym_Ors) == 1

                        ParInds = Curr_cutEBSD;
                        CutInds = [];

                    % If steel, search through the entire parent-twin system to
                    % locate fcc annealing twins that could have formed during
                    % the formation of austenite.
                    elseif length(Sym_Ors) == 5

                        % Now create second graph based solely on the first cut
                        PT_DiGraph = digraph;

                        % Create weights related to likelihood of parent
                        % orientation
                        ParOr_guess=symmetrise(Guess_Or)*R2T;
                        From_Par_ODF=calcODF(ParOr_guess,'kernel',psi);
                        ParWts=(eval(From_Par_ODF,ActTrans_Ors(Curr_cutActive))+OP(1)).*OP(2);

                        % Establish new EBSD dataset based on parent-twin system generated
                        % above
                        PT_ebsd = ActTrans_ebsd(Curr_cutActive);

                        % Call function to compute adjacency array and corresponding
                        % misorientations for the neighboring points
                        [adjpts,mori,~] = adjpt_moris(PT_ebsd,1);

                        % Evaluated misorientation distribution function values (serving as
                        % in-plane weights)
%Austin                        modf_vals=eval(myEBSD.MODF,mori);
                        modf_vals=alt_eval(myEBSD.MODF,mori);
                        modf_vals(modf_vals<0)=0;
                        PT_IPwts = (modf_vals+IP(1)).*IP(2);

                        % Call function to set up the initial graph
                        [PT_DiGraph,endnodePT,sinknodePT] = graph_setup(PT_DiGraph,adjpts,PT_IPwts,ParWts,2);

                        % Twin counter and twin indexing arrays
                        twn_count = 1;  % Corresponds to the four possible twins for each parent
                        twn_accpt = 1;  % Corresponds to the accepted (cut) twins
                        LenTwnInd = 0;

        %                 % Loop through our four twin systems to fill in possible twin
        %                 % orientations within parent-twin system
                        while twn_count < 5

                            % Temporary twin, odf, and corresponding OP weights
                            tmp_twn = symmetrise(Twn_guesses(twn_count))*R2T;    
                            tmpODF = calcODF(tmp_twn,'kernel',psi);
                            tmpWts = ((eval(tmpODF,ActTrans_Ors(Curr_cutActive))+OP(1)).*OP(2));
                            % Add edges corresponding to the guess reconstruction orientation    
                            PT_DiGraph=rmedge(PT_DiGraph,(1:len_Cpts)+endnodePT,sinknodePT);
                            PT_DiGraph=addedge(PT_DiGraph,(1:len_Cpts)+endnodePT,sinknodePT,tmpWts);

                            % Perform graph cut 
                            [mf_cPT,gfPT,csPT,ctPT]=maxflow(PT_DiGraph,1,sinknodePT);

                            % Copy the cut for convenience
                            ctcopyPT = ctPT;
                            ctcopyPT(end)=[];

                            % Index of newly cut points
                            PT_CTind = ctcopyPT-endnodePT;

                            % Now re-index for the entire ebsd scope and collect the uncut
                            % points as well
                            PT_cut = Curr_cutEBSD(PT_CTind);

                            % If the cut seems reasonable (is larger than 500 points),
                            % allow it
                            ParWts(PT_CTind) = tmpWts(PT_CTind);
                            SymIndsFlg(PT_CTind) = twn_count;
                            TwnEBSDInds{twn_count,1} = PT_cut;
                            TwnLocInds{twn_count,1} = PT_CTind;
                            % Continue on to next twin
                            twn_count = twn_count+1;
                        end

                        % Delete the twin indices, leaving only the parent indices
                        ParInds = SymInds((SymIndsFlg==0));                   
                        c2g_wts(Curr_cutEBSD) = ParWts;
                        clear TwnOrs TwnID
                        twn_accpt = 0;                
                        % Loop through twins and index correspondingly
                        for twn = 1:4
                            EBSDInds = TwnEBSDInds{twn};
                            LocInds = TwnLocInds{twn};
                            if length(EBSDInds) > 25
                                twn_accpt = twn_accpt+1;
                                Rec_ebsd.orientations(EBSDInds) = Twn_guesses(twn);
                                CutInds{twn_accpt,1} = EBSDInds;
                                TwnOrs(twn_accpt) = Twn_guesses(twn);
                                TwnID(twn_accpt) = twn;
                            end
                        end
                    end

                    % Determine whether twins exist within parent grain and index
                    % responsibly
                    if isempty(CutInds)
                        [Parent,ParCount] = ParAssign(Parent,ParInds,ParCount,Trans_ebsd,Guess_Or);
                        ParFlg = 1;
                        TwnFlg = 0;
                    else
                        [Twin,TwnCount] = TwnAssign(Twin,TwnCount,Trans_ebsd,Curr_cutEBSD,Guess_Or,TwnOrs,TwnID);
                        TwnFlg = 1;
                        ParFlg = 0;
                    end

                    % Indicate array of indices of post-transformation points that have
                    % transformed to the pre-transformation phase and those that have
                    % yet to be transformed
                    PostTransID = find(PhaseIDs==TransID);
                    PreTransID = find(PhaseIDs==ReconID);

                    % Eliminate the cut nodes in the active graph 
                    ActiveTransInd(Curr_cutActive) = [];

                    % Reset counter to 1, indicating a successful cut has been
                    % made
                    NoCutCount = 1;
                    
                    % If the unique orientation list contains something, call it
                    % here and delete the orientation used. If this becomes empty,
                    % we will not revert back to the previous method where we
                    % entertain martensite space again
                    if SkipMart
                        
                        % ID the quadrant again and the corresponding
                        % indices for this new cut that needs to be merged
                        quadId = UnqOrsId(1,1);
                        MergeInds = Rec_ebsd(Curr_cutEBSD).id;
                        
                        % Location of the overlapping grain within the
                        % previous quadrants structure
                        QId = UnqOrsId(1,4);
                        
                        % If a twin exists, things get weird. Identify the
                        % new orientations and save in an array, then flag
                        % this newly assigned grain as irrelavant when
                        % stitching everything back together at the end
                        if TwnFlg
                            lenT = length(TwnOrs);
                            for qq = 1:lenT
                                QuadOrs(qq+1) = TwnOrs(qq);
                            end
                            Twin.Flg{TwnCount} = 1;
                            RecTwn{quadId}.Flg{QId} = 0;
                        elseif (ParFlg && UnqOrsId(1,2) == 1)
                            Parent.Flg{ParCount} = 1;
                            RecPar{quadId}.Flg{QId} = 0;
                        end
                        
                        % Now Determine whether the connecting grain was
                        % from a parent or twin
                        if UnqOrsId(1,2) == 1 % Flag for parent grain
                            NewInds = RecPar{quadId}.Indices{QId};
                            NewInds = vertcat(NewInds,MergeInds);
                            % If only a parent, we only change the indices
                            % recorded in the previous quadrant
                            if ParFlg
                                RecPar{quadId}.Indices{QId} = NewInds;
                                Parent.Indices{ParCount} = NewInds;
                            % If a twin exists, we need to reindex
                            % everything and add in the new information
                            elseif TwnFlg
                                
                                RecPar{quadId}(QId) = [];
                                ParCount = ParCount-1;
                                Ltwn = length(RecTwn{quadId}.Merged); 
                                RecTwn{quadId}.Merged{Ltwn+1} = NewInds;
                                RecTwn{quadId}.Parent.Or{Ltwn+1} = Guess_Or;
                                for qq = 1:lenT
                                    RecTwn{quadId}.Single.Or{Ltwn+1}(qq) = QuadOrs(qq+1);
                                    RecTwn{quadId}.Single.ID{Ltwn+1}(qq) = TwnID(qq);
                                end
                                % Now reclassify the twin structure
                                MrgEb = FullEbsd(NewInds);
                                MrgInds = 1:length(MrgEb);
                                [Twin,TwnCount] =...
                                    TwnAssign(Twin,TwnCount,MrgEb,MrgInds,Guess_Or,TwnOrs,TwnID);
                            end
                        % Flag indicating the cut portion connected to a
                        % previous quadrant was a twin
                        elseif UnqOrsId(1,2) == 2 
                            % Merge parent-twin structure across quadrants 
                            NewInds = RecTwn{quadId}.Merged{QId};
                            NewInds = vertcat(NewInds,MergeInds);
                            RecTwn{quadId}.Merged{QId} = NewInds;
                            
                            % Maximum number of twins that have been
                            % identified within a single PAG (max == 4)
                            szT = size(RecTwn{quadId}.Single.Or);
                            
                            % If indicated as a parent, delete the local
                            % quadrant parent structure and reclassify as a
                            % Twin structure
                            if ParFlg
                                ParCount = ParCount-1;
                                for qq = 1:szT(2)
                                    
                                    TmpOr = RecTwn{quadId}.Single.Or{QId,qq};
                                    TmpID = RecTwn{quadId}.Single.ID{QId,qq};
                                    if isempty(TmpOr)==0
                                        TwnOrs(qq) = TmpOr;
                                        TwnID(qq) = TmpID;
                                    end
                                end
                            else
                                TwnCount = TwnCount-1;
                            end
                            % Reclassify the local twin quadrant structure
                            MrgEb = FullEbsd(NewInds);
                            MrgInds = 1:length(MrgEb);
                            [Twin,TwnCount] =...
                                TwnAssign(Twin,TwnCount,MrgEb,MrgInds,Guess_Or,TwnOrs,TwnID);
                            
                            % If the orientation passed from connected
                            % quadrant came from the parent
                            if UnqOrsId(1,3) == 1
                                % If this one contains a twin, need to
                                % reindex our structure to accommodate this
                                if TwnFlg
                                    tcnt = 2;
                                    for qq = 1:lenT
                                        TOr = QuadOrs(qq+1);
                                        TId = TwnID(qq);
                                        twnchk = [];
                                        % Check previous quadrant structure
                                        % for all twins
                                        for tt = 1:szT(2)
                                            RecTwnOr = RecTwn{quadId}.Single.Or{QId,tt};
                                            if isempty(RecTwnOr)==0
                                                twnchk(tt) = TOr == RecTwnOr;
                                            end
                                        end
                                        % Only add to structure if the new
                                        % orientation is unique
                                        if all(twnchk==0)
                                            RecTwn{quadId}.Single.Or{QId,tcnt} = TOr;
                                            RecTwn{quadId}.Single.ID{QId,tcnt} = TId;
                                            tcnt=tcnt+1;
                                        end
                                    end
                                end
                            % Orientation passed over from connecting
                            % quadrant came from the twinned grain
%                             else 
%                                 % If a new twin was flagged, we need to add
%                                 % it to the local twin quadrant structure
%                                 % only if it's not the parent
%                                 if TwnFlg
%                                     % Find out how many twins exist in
%                                     % previous structure
%                                     go = 1;
%                                     while go <= szT(2)
%                                         twnchk = RecTwn{quadId}.Single.Or{QId,qcount};
%                                         if isempty(twnchk)
%                                             qcount = tt;
%                                             go = szT(2)+1;
%                                         end
%                                     end
%                                     % If the twin orientation was sent
%                                     % over, we need to make sure a flagged
%                                     % 'twin' isn't just the rest of the
%                                     % parent grain
%                                     for qq = 1:lenT
%                                         TOr = QuadOrs(qq+1);
%                                         TId = TwnID(qcount);
%                                         if TOr ~= RecTwn{quadId}.Parent.Or{QId,qq}
%                                             RecTwn{quadId}.Single.Or{QId,qcount+1} = TOr;
%                                             RecTwn{quadId}.Single.ID{QId,qcount+1} = TId;
%                                             qcount = qcount+1;
%                                         end
%                                     end
%                                     % If there wasn't a twin, reclassify
%                                     % the 'Parent' structure as a 'Twin'
%                                     % instead for proper bookkeeping.
%                                     if lenT == 0
%                                         ParCount = ParCount-1;
%                                         MrgEb = FullEbsd(NewInds);
%                                         MrgInds = 1:length(MrgEb);
%                                         [Twin,TwnCount] = TwnAssign(Twin,TwnCount,MrgEb,MrgInds,Guess_Or,TwnOrs,TwnID);
%                                     end
%                                 end 
                                
                            end % Determine type of twin         
                        end % Grain type within SkipMart loop
                        UnqOrs(1) = [];
                        UnqOrsId(1,:) = [];
                    end % SkipMart loop
                    
                end % Austenite space cut
            end % Martensite space cut
            
            % Total iteration count across quadrants
            TotCount = TotCount+1;

            % Display iterations to show progress of code
            display(sprintf('---------------Iteration #: %d---------------', TotCount))
            
            % Determine whether to continue on or end the current run
            if length(PostTransID) < 25
                EndFlag = 1;
            else
                % Truncate Transformation orientations and EBSD dataset
                ActTrans_ebsd = Trans_ebsd(PostTransID);
                ActTrans_Ors = Trans_ebsd(PostTransID).orientations;
                % Compute updated adjacent point
                [Act_adjpts,Act_mori,EndFlag] = adjpt_moris(ActTrans_ebsd,1);
            end

            % EndFlag represents either empty ebsd data set (fully transformed)
            % or sets that have no adjacent points (noise)
            if EndFlag == 1 || NoCutCount==ConsNoCuts
                iters = num_iters;
            else
                % Active inplane and first set of out of plane weights
%Austin                ActMODF_vals=eval(myEBSD.MODF,Act_mori);
                ActMODF_vals=alt_eval(myEBSD.MODF,Act_mori);
                ActMODF_vals(ActMODF_vals<0)=0;
                ActIPwts = (ActMODF_vals+IP(1)).*IP(2);
                Actc2gWts = c2g_wts(PostTransID);

                % Reset the active graph
                ActDiGraph = digraph;
                % Recreate an active graph that neglects all nodes that have 
                % already been cut
                [ActDiGraph,ActEndnode,ActSinknode] = graph_setup(ActDiGraph,Act_adjpts,ActIPwts,Actc2gWts,2);
            end          
        end
        % Store the reconstructed ebsd, weights, parent and twin
        % structures, and phase ids for the previous quadrant
        RecEbsd{ii} = Rec_ebsd;
        RecWts{ii} = c2g_wts;
        RecPar{ii} = Parent;
        RecTwn{ii} = Twin;
        RecPhIds{ii} = PhaseIDs;
        tp_count = 1;
        T1Ors = orientation('euler',[0,0,0],myEBSD.CS{3});
        T1Inds = 0;
        
        if ii ~= nQds
            
            % Now, for the next quadrant (if more than one exists), let's pull 
            % out the austenite orientations associated with a previous,
            % adjacent edge (if they exist) and use as our first guesses in the 
            % next quadrant.
            [rq,~] = find(myEBSD.Quad.IDs(:,1:(end-1))==ii+1);
            Qedgs = myEBSD.Quad.IDs(rq,:);
            Qedgs(Qedgs(:,2)>ii+1,:)=[];
            szQ = size(Qedgs);
            UnqOrs = [];
            UnqOrsId = [];
            Ptmpvec = [];
            Ttmpvec = [];

            % Fill in arrays that will be used in the next quadrant
            for E = 1:szQ(1)
                
                % Assign the ID of the adjacent quadrant and Ebsd
                AdjQd = Qedgs(E,1);
                tmpEb = RecEbsd{AdjQd};
                % Now extract the corresponding edges
                AdjEdg = myEBSD.Quad.Edges{AdjQd,Qedgs(E,end)};
                % Find any points in edge that were assigned transformation
                % phase and extract these points
                QdRecPh = find(RecPhIds{AdjQd}(AdjEdg)==ReconID);
                
                % If the orientation is found, save the quadrant it was 
                % found in (1st index), whether it was a parent  or twin 
                % (2nd), if twin, whether orientation belongs to a parent 
                % or twin (3rd), the location within the structure (4th) 
                % and finally the position within the 'UnqOrs' vector (5th)
                if isempty(QdRecPh)==0
                    EdgOrs = tmpEb(AdjEdg(QdRecPh)).orientations;
                    UnqOrs = vertcat(UnqOrs,unique(EdgOrs));
%                     UnqOrs = unique(UnqOrs);
                    
                    % A twinned grain exists from previous quadrant. Find 
                    % out whether bordering orientation is parent or twin
                    if isempty(RecTwn{AdjQd})==0
                        for jj = 1:length(RecTwn{AdjQd}.Merged)
                            T1tmp = find(UnqOrs==RecTwn{AdjQd}.Parent.Or{jj});
                            if isempty(T1tmp)==0
                                TmpTParOr = RecTwn{AdjQd}.Parent.Or{jj};
                                if any(T1Ors==TmpTParOr)==0
                                    T1Ors(tp_count) = RecTwn{AdjQd}.Parent.Or{jj};
                                    T1Inds(tp_count,1) = AdjQd;
                                    T1Inds(tp_count,2) = jj;
                                    tp_count = tp_count+1;
                                end
                            end
                            
                            % Determine if a twin (or multiple twins from
                            % same parent) overlap from previous quadrant 
                            szT = size(RecTwn{AdjQd}.Single.Or);
                            twincheck = 0;
                            tcount = 0;
                            for tt = 1:szT(2)
                                chktwnOr = RecTwn{AdjQd}.Single.Or{jj,tt};
                                if isempty(chktwnOr)==0
                                    tmp = find(UnqOrs==chktwnOr);
                                    if isempty(tmp)==0
                                        tcount = tcount+1;
                                        T2tmp(tcount) = tmp;
                                    end
                                end
                            end
                            
                            % Flag an overlapping twin
                            if tcount > 0
                                twincheck = 1;
                                % If multiple twins exist within same
                                % grain, remove superfluous orientations
                                % from UnqOrs matrix
                                if length(T2tmp) > 1
                                    for tt = 2:length(T2tmp)
                                        UnqOrs(T2tmp(tt)) = [];
                                    end
                                    T2tmp(2:end) = [];
                                end
                            end
                            
                            % A twin orientation resides at the 
                            % quadrant border
                            if twincheck==1   
                            % If twin orientation is overlapping, not the 
                            % parent, choose parent orientation instead
                                if isempty(T1tmp)
                                    UnqOrs(T2tmp) = RecTwn{AdjQd}.Parent.Or{jj};
                                    T1Ors(tp_count) = RecTwn{AdjQd}.Parent.Or{jj};
                                    T1Inds(tp_count,1) = AdjQd;
                                    T1Inds(tp_count,2) = jj;
                                    tp_count = tp_count+1;
                                    T1tmp = T2tmp;
                                % If both parent and twin(s) are on border, 
                                % delete the twin and reindex
                                else
                                    UnqOrs(T2tmp) = [];
                                end
                            end
                            clear T2tmp T1tmp
                        end
                    end
                    
                    % If untwinned parent grain exists within the previous
                    % quadrant, check to see if it overlaps into next one
                    if isempty(RecPar{AdjQd}) == 0
                        for jj = 1:length(RecPar{AdjQd}.Or)
                            Ptmp = find(UnqOrs==RecPar{AdjQd}.Or{jj});
                            if isempty(Ptmp)==0
                                if length(Ptmp)>1
                                    UnqOrs(Ptmp(2:end)) = [];
                                    Ptmp(2:end)=[];
                                end
                                Ptmpvec = vertcat(Ptmpvec,[AdjQd,1,0,jj,Ptmp]);
                            end
                        end
                    end
                end
                
            end
            
            % Finalize indexing for twinned-grain  orientations 
            % extending across the quadrant border
            if T1Ors(1) ~= orientation('euler',[0,0,0],myEBSD.CS{3});
                for kk = 1:length(T1Ors)
                    T1tmp = find(UnqOrs==T1Ors(kk));
                    Ttmpvec = vertcat(Ttmpvec,[T1Inds(kk,1),2,1,T1Inds(kk,2),T1tmp]);
                end
            end
            
            % Sort everything in the UnqOrsId vector to pass on to the
            % next quadrant reconstruction. This vector holds all the
            % information needed to restructure the finalized 
            % Parent/Twin structures in case a merge needs to happen
            UnqOrsId = vertcat(Ptmpvec,Ttmpvec);
            
            % Sort UnqOrsId so it matches up with the ordering of UnqOrs. This
            % will allow us to systematically delete variables for easier
            % bookkeeping
            if isempty(E)==0
                [~,srtUnq] = sort(UnqOrsId(:,end));
                tmpUnq = UnqOrsId(srtUnq,:);
                UnqOrsId = tmpUnq;
            end
        end
    end
    
    %%
    % Now put humpty-dumpty back together again
%Austin addin - Get rid of unneeded stuff to help with error checking
    clear Twin twn twn_accpt twn_count Twn_eul Twn_guesses TwnCount
    clear Twns TwnEBSDInds TwnFlg TwnLocInds UnqOrs
    clear TransEBSD_ci TransCT_ebsd Trans_MODF tp_count tmpWts tmpODF tmp_ci tmp_twn
    clear Act_adjpts Act_mori Actc2gWts ActDiGraph ActEndnode ActIPwts ActiveTransInd ActMODF_vals ActOrs
    clear ActSinknode ActTrans_ci ActTrans_ebsd ActTrans_Ors
    clear c2g_wts c2gWts adjpts ci ci_id ConsNoCuts cs cs2 csPT
    clear ct ctw ctcopy ctcopy2 ctcopyPT ctPT ct2 Curr_cutActive Curr_cutEBSD
    clear Curr_cutActive Curr_cutEBSD CutInds DiGraph EbsdInd EBSDInds EndFlag
    clear Endnode endnode2 endnodePT From_guess_ODF From_Or_guess From_Par_ODF
    clear mf_c mf_c2 mf_cPT Misos modf_vals mori Sinknode sinknode2 sinknodePT
    clear SkipMart psi Rec_ebsd num_iters New_Digraph NewCT_ind ParFlg ParInds
    clear PT_CTind PT_cut PT_DiGraph PT_PIwts PT_ebsd PT_IPwts ParWts Par_Or
    clear OPWts OP TransEbsd_ci Trans_ebsd TotCount T2R T1Inds T1Ors SymInds SymIndsFlg
    clear R2T PreTransID PostTransID ParOr_guess Parent ParCount Orien NoCutCount
    clear gfPT IP_wts IPwts IP ii iters len_Cpts len_preCT LenTwnInd local_odf LocInds
    
    ParCnt = 1;
    TwnCnt = 1;
    TwnOrLst = orientation('euler',[0,0,0],CS_R);
    FullParent = [];
    FullTwin = [];
    for jj = 1:nQds
       % Add reconstructed quadrants to entire microstructure and properly
       % arrange the weights and phase ids
       FullEbsd(RecEbsd{jj}.id) = RecEbsd{jj};
       FullWts(RecEbsd{jj}.id) = RecWts{jj};
       FullPhIds(RecEbsd{jj}.id) = RecPhIds{jj};
       % Establish the parent and twin structures for each quadrant
       tmpPar = RecPar{jj};
       tmpTwn = RecTwn{jj};
       
       % Loop through parents first and, if flagged, add to the structure
       % containing all of the parent grains
       if isempty(RecPar{jj})==0
           for kk = 1:length(RecPar{jj}.Flg)
                if RecPar{jj}.Flg{kk}==1
                    FullParent.Indices{ParCnt,1} = RecPar{jj}.Indices{kk};
                    FullParent.Or{ParCnt,1} = RecPar{jj}.Or{kk};
                    ParCnt = ParCnt+1;
                end
           end
       end
       % Loop through twins to do the same
       if isempty(RecTwn{jj})==0
           for ll = 1:length(RecTwn{jj}.Flg)
               if RecTwn{jj}.Flg{ll}==1
                   FullTwnFlg = 1;
                   szF = size(RecTwn{jj}.Single.Or);
                   FinTParOr = RecTwn{jj}.Parent.Or{ll};
                   FinTParInds = RecTwn{jj}.Merged{ll};
                   % Finalize one last flag to make sure no twins are being
                   % repeated across quadrants and added to the structure
                   % containing all parent-twin grains
                   if any(TwnOrLst == FinTParOr)
                       RedunTwn = find(TwnOrLst == FinTParOr);
                       Tnew = length(FinTParInds);
                       Told = length(FullTwin.Merged{RedunTwn,1});
                       % Assume structures with the most indices are updated,
                       % and keep these as the true parent-twin structure for
                       % the reconstruction
                       if Tnew > Told
                           FullTwin.Merged(RedunTwn) = [];
                           FullTwin.Parent.Or(RedunTwn) = [];
                           FullTwin.Single.Or(RedunTwn) = [];
                           FullTwin.Single.ID(RedunTwn) = [];
                           TwnCnt = TwnCnt-1;
                       else
                           FullTwnFlg = 0;
                       end
                   end
                   if FullTwnFlg
                       TwnOrLst(TwnCnt) = FinTParOr;
                       FullTwin.Merged{TwnCnt,1} = FinTParInds;
                       FullTwin.Parent.Or{TwnCnt,1} = FinTParOr;
                       for mm = 1:szF(2)
                           FullTwin.Single.Or{TwnCnt,mm} = RecTwn{jj}.Single.Or{ll,mm};
                           FullTwin.Single.ID{TwnCnt,mm} = RecTwn{jj}.Single.ID{ll,mm};
                       end
                       TwnCnt = TwnCnt+1;
                   end
               end
           end
       end
    end
    
    % Assign phase ids for entire microstructure
    FullEbsd.phase = FullPhIds;

    % If we have untransformed points and our initial ebsd dataset had more
    % than one phase, assign points
    if length(unique(PhaseIDs)) > 1 && n_phases > 1
       find_untrans = find(FullPhIds==TransID);
       FullEbsd(find_untrans).CS = CS_T;
    end

    % Sort all Parent and Twin Indices into a single array for indexing
    % purposes
    ParTot_Inds = [];
    TwinTot_Inds = [];

    if isempty(FullParent)
       ParTot_Inds = [];
    else
       for mm = 1:length(FullParent.Indices)
           ParTot_Inds = vertcat(ParTot_Inds,FullParent.Indices{mm});
       end
    end

    if isempty(FullTwin)
       TwinTot_Inds = [];
    else
       for nn = 1:length(FullTwin.Merged)
           TwinTot_Inds = vertcat(TwinTot_Inds,FullTwin.Merged{nn});
       end
    end

    % Assign the reconstructed points, which only amount to those first
    % identified as the post-transformation phase being assigned to the
    % pre-transformation phase
    Recon_pts = find(FullPhIds==ReconID);
    ReconIds = FullPhIds;
    ReconWts = FullWts;
    ReconEbsd = FullEbsd;
    PreTransInd = find(myEBSD.Ebsd.phase==ReconID);
    % Austin Edit
    % Recon_pts(PreTransInd) = [];
    ReconIds(PreTransInd) = [];
    ReconWts(PreTransInd) = [];
    ReconEbsd(PreTransInd) = [];

    % Add variables to myEBSD structure
    myEBSD.Recon.FullEbsd = FullEbsd;
    myEBSD.Recon.Ebsd = ReconEbsd;
    myEBSD.Recon.TransformedPoints = Recon_pts;
    myEBSD.Recon.PhaseIDs = ReconIds;
    myEBSD.Recon.Likelihood = ReconWts;

    % Add variables to Parent and Twin structures
    FullParent.AllIndices = ParTot_Inds;
    FullTwin.AllIndices = TwinTot_Inds;
    
end
