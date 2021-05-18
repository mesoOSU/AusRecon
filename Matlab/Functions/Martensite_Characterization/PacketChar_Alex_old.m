function [Packet] = PacketChar(myEBSD,Grains,gId,Packet,Twin)
    %% Intialization
    
%     GrnId = find(Grains.Ebsd.Parent.ID==GrainID);
    Aus_Grain   = myEBSD.Recon.Ebsd(Grains.Indices{gId});
    MartRecon   = myEBSD.Ebsd;
    CS_T        = myEBSD.CS{2};
    CS_R        = myEBSD.CS{3};
    TransID     = myEBSD.Phase.ID{1};               
    ReconID     = myEBSD.Phase.ID{2}; 
    SS          = myEBSD.SS;
    ksi         = myEBSD.OR;
    halfwidth   = myEBSD.noise.halfwidth;
    psi         = myEBSD.noise.psi;
    T2R         = myEBSD.TransMats{1};
    R2T         = myEBSD.TransMats{2};
    %% Start Characterization of Martensite Packet Boundaries
    
    if Grains.grainId(gId,2) ~= TransID
        
        % Interface parameters (constant for now, see if these do change
        % much with variable steel and ferrous alloy data sets)
        IP = 5;
        OP = 2e-1;

        % Make sure martensite grain doesn't have residual austenite
        AusIDs            = find(MartRecon.phase==ReconID);
        MartRecon(AusIDs) = [];
        Mart_Grain        = MartRecon(Grains.Indices{gId});
        
        count = 0;           
        twnnum = unique(Grains.grainId(gId,3:6));
        
        twnnum(1) = [];
        V = [];
        G = [];
        tmpInds = [];
        
        while count <= length(twnnum)
            count = count+1;
            % Assign the appropriate austenite orientation (parent or twin)
            if Grains.grainId(gId,2) == ReconID
                % Assign austenite orientation and rotate it to Transformation crystal
                % frame
                AusOr = Aus_Grain.orientations(1);
            else
                % Twin Id
                tId = Grains.grainId(gId,end);
                if count == 1
                    AusOr = Twin.Parent.Or{tId};
                else
                    AusOr = Twin.Single.Or{tId,count-1};
                end
            end
            
            % Extract the parent or twin (if it exists) indices
            tmpInds{count} = find(Aus_Grain.orientations == AusOr);
            
            % Make the chosen austenite orientation a rotation
            Aus_or(count) = rotation('Euler',AusOr.phi1,AusOr.Phi,AusOr.phi2,CS_T);

            % Compute variants and corresponding groupoid from euler angles
            [Vtmp,~] = YardleyVariants(ksi);
            [Gtmp,~] = GroupoidVariantsCubic(Vtmp);
            
            % Now concatenate the variant and groupoid (misorientation)
            % vectors to account for the existance of any twins
            V = vertcat(V,Vtmp);
            G = vertcat(G,Gtmp);
            
        end
        
        % Set up structure of weighted values
        WTs = zeros(length(Aus_Grain),5);
        for nn = 1:count
            
            % Since result is a passive rotation, take transpose of variant
            % transformation indices to convert to active (Consistent with MTEX)
            % NOTE TOREN: Cryspy is all passive, so skip this transpose step when 
            % converting reconstruction code for cryspy. Then convert the
            % transposed variants to euler angles.
            
            Vst = 24*(nn-1)+1;
            Vnd = 24*nn;
            for i = Vst:Vnd
                Vt{i,1} = transpose(V{i});
                VTeul(i) = orientation('matrix',Vt{i},CS_T);
            end

            % Convert group axis to actual axis and then to Euler angles
            Ax = vector3d(G(:,2),G(:,3),G(:,4));
            G_eul = orientation('axis',Ax,'angle',G(:,1),CS_T);
            % Kappa value
            kappa=psi.kappa;

            martOrs = Mart_Grain(tmpInds{nn}).orientations;
            len = length(tmpInds{nn});
            % Initialize the graph
            grn_graph = digraph;
            
            % Indexing start and finish for pertinent misorientations
            Gst = 16*(nn-1)+1;
            Gnd = Gst+3;

            % Create misorientation distribution functions for block and packet
            % boundaries, respectively
            mori_block = symmetrise(G_eul(Gst:Gnd),CS_T);
            mori_pack = symmetrise(G_eul(Gnd+1:nn*16),CS_T);
            martblock_modf=calcODF(mori_block,'kernel',psi);
            martpack_modf = calcODF(mori_pack,'kernel',psi);

            % Create adjacency array based on specific interactions between points
            % within the input grain
            [grn_adj,grn_mori,~] = adjpt_moris(Mart_Grain(tmpInds{nn}),1);

            % Assign a sample symmetry to the grain
            grn_mori.SS = SS;

            % Evaluate misorientations with respect to block-specific values
            martblock_modfvals = eval(martblock_modf,grn_mori);
            martblock_modfvals(martblock_modfvals<0)=0;
            % Now set up our graph
            grainIP_wts = martblock_modfvals.*IP;

            Pack_Bounds = zeros(len,1);
            Pack_Pts = [];

            for iter = 1:4
                % ID for group of 6 variants being analyzed and the rest of them
                tmp1 = (iter-1)*6+1;
                tmp2 = (iter)*6;

                % Pull variants out
                Varsg2t = VTeul(tmp1:tmp2);
                Varsc2g = VTeul;
                Varsc2g(tmp1:tmp2) = [];

                % Rotate variant orientations by respective austenite orientation
                Var_rotg2t = Aus_or(nn) * Varsg2t;
                Var_rotc2g = Aus_or(nn) * Varsc2g;

                % Compute ODF for both sets of variants
                Varg2tODF = calcODF(symmetrise(Var_rotg2t),'kernel',psi);
                Varc2gODF = calcODF(Var_rotc2g,'kernel',psi);

                % Compute the c2g weights (first set of OP weights)
                c2g_wts = eval(Varc2gODF,martOrs);
                c2g_wts(find(c2g_wts<0))=0;

                % Compute the g2t weights (second set of OP weights)
        %         g2t_wts = (eval(Varg2tODF,martOrs)+OP(1)).*OP(2); 
                g2t_wts = eval(Varg2tODF,martOrs);
                g2t_wts(find(g2t_wts<0))=0;

                % Only if we're on the first iteration should graph be
                % constructed
                if iter==1
                    % Setup graph within the single austenite grain boundaries
                    [grn_graph,grn_endnode,grn_sinknode] = graph_setup(grn_graph,grn_adj,grainIP_wts,c2g_wts,2);
                end

                % Add g2t weights to graph
                grn_graph=rmedge(grn_graph,(1:len)+grn_endnode,grn_sinknode);
                grn_graph=addedge(grn_graph,(1:len)+grn_endnode,grn_sinknode,g2t_wts);

                % Commence graph cutting
                [mf_c,gf,cs,ct]=maxflow(grn_graph,1,grn_sinknode);
                ctcopy = ct;
                ctcopy(end)=[];
                grn_ind = ctcopy-grn_endnode;

                % Assign cut out block boundaries, indices and corresponding
                % likelihoods
                Pack_Weights{iter} = g2t_wts;

                % Remove necessary c2t weights and add again
                grn_graph=rmedge(grn_graph,(1:len)+1,(1:len)+grn_endnode);
                grn_graph=addedge(grn_graph,(1:len)+1,(1:len)+grn_endnode,c2g_wts);

            end
            
            % Vector of the four corresponding packet boundary weights with a zero
            % column
            Weights = zeros(5,len);
            for i = 2:5
                Weights(i,:) = Pack_Weights{i-1};
            end

            % Now fill in the graph for indices of parents and twins (if
            % they exist)
            WTs(tmpInds{nn},:) = Weights';
        end
%     end
        %% ``Denoising'' Section (Don't worry about)
        
            % Identify max weights
            [max_wts,maxid_wts] = max(WTs');
            zerwts = find(max_wts==0);

            % Assign a packet boundary to each point based on the maximized weights
            if isempty(zerwts)
                maxid_wts = maxid_wts.*10;
            else
                maxid_wts = (maxid_wts-1).*10;
            end

            % Set denoising weights to the same ascribed above in initial cut
            DeNoiseWTs = WTs;
            
            clear grn_adj
            
            % Create adjacency array based on specific interactions between points
            % within the input grain
            [grn_adj,grn_mori,~] = adjpt_moris(Mart_Grain,1);
                       
            % For all zeroed values (noise and otherwise), find adjacent points and
            % reassign new weighted value as the mean between its surrounding
            % points weights
            
            for jj = 1:length(zerwts)
                tmpid = zerwts(jj);
                [r,~]=find(grn_adj==tmpid);
                unqadj = unique(grn_adj(r,:));
                unqadj(unqadj==tmpid)=[];               
                DeNoiseWTs(tmpid,:) = mean(WTs(unqadj,:));
            end

            % Take the max of the denoised weights to use as our inplane weights
            [maxnew_wts,maxnew_wtsID] = max(DeNoiseWTs');

            % Create new graph structure
            DeNoise_GrnGraph = digraph;
            len = length(Aus_Grain);

            % Identify the indexing for the maxed weights, and our new inplane
            % weights become the inverted difference between assigned boundary
            % value for each adjacent pair
            maxDN_new_wts1 = maxnew_wtsID(grn_adj(:,1));
            maxDN_new_wts2 = maxnew_wtsID(grn_adj(:,2));
            DN_grnIP_wts = IP./abs(maxDN_new_wts1-maxDN_new_wts2);
            DN_grnIP_wts(DN_grnIP_wts==Inf) = 1e2;

            % First set of OP weights are scaled to 10 
            DNc2g_wts = ones(length(DeNoiseWTs),1).*10;

            % Create graph
            [DeNoise_GrnGraph,grn_endnode,grn_sinknode] =...
                graph_setup(DeNoise_GrnGraph,grn_adj,DN_grnIP_wts,DNc2g_wts,2);
            
            Pack_Bounds = zeros(len,1);

            % Perform denoising cuts
            for ii = 1:4
                % DNg2t_wts are just scaled original weights
                DNg2t_wts = DeNoiseWTs(:,ii+1).*OP;

                % Add g2t weights to graph
                DeNoise_GrnGraph = ...
                    rmedge(DeNoise_GrnGraph,(1:len)+grn_endnode,grn_sinknode);
                DeNoise_GrnGraph = ...
                    addedge(DeNoise_GrnGraph,(1:len)+grn_endnode,grn_sinknode,DNg2t_wts);

                % Commence graph cutting
                [mf_c,gf,cs,ct]=maxflow(DeNoise_GrnGraph,1,grn_sinknode);
                ctcopy = ct;
                ctcopy(end)=[];

                % Indices for cut
                DNgrn_ind = ctcopy-grn_endnode;

                % Assign  boundary values and weights corresponding to cuts
                Pack_Bounds(DNgrn_ind) = ii;
                Pack_Pts{ii} = DNgrn_ind;
                DN_PackWts{ii} = DNg2t_wts;
            end
            
            FinalWts = zeros(length(WTs),1);
            for i = 1:length(WTs)
                FinalWts(i) = WTs(i,Pack_Bounds(i)+1);
            end
    else
            Mart_Grain = [];
            Pack_Bounds = zeros(length(Aus_Grain),1);
            Pack_Pts = [];
            FinalWts = zeros(length(Aus_Grain),1);
            maxid_wts = zeros(length(Aus_Grain),1);
    end
    % Fill in an EBSD data set with the packet boundaries

        % Assign pertinent variables to our ``Packet'' structure
        Packet{gId}.MartGrains  = Mart_Grain;
        Packet{gId}.AusGrain    = Aus_Grain;
        Packet{gId}.Boundaries  = Pack_Bounds;
        Packet{gId}.BoundaryIDs = Pack_Pts;
        Packet{gId}.Weights     = FinalWts;
        Packet{gId}.MaxBounds   = maxid_wts';
        Packet{gId}.Grain       = Grains.grains(gId);
        
        [Packet] = BlockChar(myEBSD,Grains,gId,Packet,Twin);
    
end

