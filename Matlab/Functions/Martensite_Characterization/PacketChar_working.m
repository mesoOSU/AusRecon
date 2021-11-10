function [Packet] = PacketChar(myEBSD,Grains,gId,Packet,Twin)
%%                      Function Description

% This function performs martensite packet segmentations on each indexed
% prior austenite grain using a similar graph cutting approach. Due to
% inherent noise in the resolution, a de-noising step is also implemented
% post-transformation to adequately cluster noisier indices into
% corresponding packets. 

% At the conclusion of the function, the determined packets are sent to 
% another function for block (and sub-block) characterization.

% Input Variables: myESBD structure; grain boundary structure; grain id
%                  array; Packet structure; Twin structure

% Output Variables: Structure containing packet information

%%                      Function Implementation

    % Intialization variables
    RecEb = myEBSD.Recon.Ebsd(myEBSD.Recon.Ebsd.isIndexed);
    AusGrn   = RecEb(Grains.Indices{gId});
    MartRecon   = myEBSD.Ebsd(myEBSD.Ebsd.isIndexed);
    CS_T        = myEBSD.CS{2};
    TransID     = myEBSD.Phase.ID{1};               
    ReconID     = myEBSD.Phase.ID{2}; 
    SS          = myEBSD.SS;
    ksi         = myEBSD.OR;
    psi         = myEBSD.noise.psi;
    
    % Begin characterization of martensite packet boundaries if the
    % submitted grain region was assigned an austenite phase (not residual,
    % unassigned martensite)
    if Grains.grainId(gId,2) ~= TransID

        % Regularization parameters
        % NOTE: These are constant for now, but you can see if these do
        % alter the results much. However, I found they were pretty 
        % consistent across different steel compositions)
        IP = 5;
        OP = 2e-1;

        % Make sure martensite variants within austenite grain don't have
        % any points marked as residual austenite
        AusIDs                = find(MartRecon.phase==ReconID);
        MartRecon(AusIDs)     = [];
        MartVars              = MartRecon(Grains.Indices{gId});
        
        % Counter and unique twin variant number (can be from -1; 1 to 4 if a 
        % twin was ID'd with in the austenite grain)
        % NOTE: The -1 index only doesn't hold iff ALL 4 twins have been 
        % identified within the hgihlighted austenite grain
        TwnCount = 0;           
        twnnum = unique(Grains.grainId(gId,3:6));
        
        % Delete -1 index if it exists (redundant)
        if any(twnnum == -1)
            twnnum(1) = [];
        end
        
        % Establish Variant (V), Groupoid (G), and temporary indexing
        % structures
        V = [];
        tmpInds = [];
        
        % Loop through all indexed twins to identify the correct variant
        % orientation 
        while TwnCount <= length(twnnum)
            TwnCount = TwnCount+1;
            % Assign the appropriate austenite orientation (parent or twin)
            if Grains.grainId(gId,2) == ReconID
                % Assign the parent orientation if no twins
                tmpOr = Grains.grains(gId).meanOrientation;
            else
                % Find the correct twin orientation if a twin exists
                tId = Grains.grainId(gId,end);
                if TwnCount == 1
                    % Assign parent orientation of parent-twin system
                    tmpOr = Twin.Parent.Or{tId};
                else
                    % Assign twin orientations 
                    tmpOr = Twin.Single.Or{tId,TwnCount-1};
                end
            end
            
            % Extract either the parent or twin EBSD in
            tmpInds{TwnCount} = find(AusGrn.orientations == tmpOr);
            
            % Make the chosen austenite orientation a rotation with
            % martensite crystal symmetry
            AusOr(TwnCount) = rotation('Euler',tmpOr.phi1,tmpOr.Phi,tmpOr.phi2,CS_T);

            % Compute variants and corresponding groupoid from euler angles       
            % NOTE: For more thorough understanding of groupoids, refer to:
            % Brust, A.F. et al, "Analysis of Misorientation Relationships
            % Between Austenite Parents and Twins" (2018)
            [Vtmp,~] = YardleyVariants(ksi);
            
            % Now concatenate the variant and groupoid (misorientation)
            % vectors to account for the existance of any twins
            V = vertcat(V,Vtmp);
            
        end

        % Set up structure of weighted values and loop through all possible
        % packets to establish likelihoods for all
        % NOTE: Graph cutting is not performed at this time, only
        % likelihoods are establish for each set of packet-specific
        % variants
        WTs = zeros(length(AusGrn),5);
        for nn = 1:TwnCount
            
            % Since result is a passive rotation, take transpose of variant
            % transformation and convert to active rotation (Consistent
            % with MTEX) and denote as euler angles
            % NOTE: Cryspy uses passive notation, so skip this transpose 
            % step when converting reconstruction code to Python. 
            Vst = 24*(nn-1)+1;
            Vnd = 24*nn;
            for i = Vst:Vnd
                Vt{i,1} = transpose(V{i});
                VTeul(i) = orientation('matrix',Vt{i},CS_T);
            end
            if gId == 5
                keyboard
            end
            % Set all martensite variant orientations within the chosen 
            % austenite grain boundaries into an array
            martOrs = MartVars(tmpInds{nn}).orientations;
            len = length(tmpInds{nn});

            % Indexing start and finish for relevant misorientations
            Gst = 16*(nn-1)+1;
            Gnd = Gst+3;

            % Pre-allocate the packet boundaries and packet indices
            Pack_Bounds = zeros(len,1);
            Pack_Pts = [];

            % Iterate through all 4 sets of packets related to either the
            % parent or twin orientation for segmentation
            for iter = 1:4
                % ID for group of 6 variants being analyzed and the rest of them
                tmp1 = (iter-1)*6+1;
                tmp2 = (iter)*6;

                % Extract the variant orientations specific to the packet
                % in question
                Varsg2t = VTeul(tmp1:tmp2);

                % Rotate variant orientations by austenite orientation
                % NOTE: These are the theoretical variants that should have
                % occured upon transformation of the chosen PAG grain
                TheoVar_g2t = AusOr(nn) * Varsg2t;

                % Compute ODF for both sets of variants
                Varg2tODF = calcODF(symmetrise(TheoVar_g2t),'kernel',psi);
                % Compute the g2t weights (second set of OP weights)
                g2t_wts = eval(Varg2tODF,martOrs);
                g2t_wts(find(g2t_wts<0))=0;

                % Assign the packet-specific weights for all indices
                Pack_Weights{iter} = g2t_wts;
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
        %% ``Denoising'' Section of Weighted Values
       
            % Identify the maximum weight w/r/t the packets for each
            % indexed node within the PAG/T
            [max_wts,maxid_wts] = max(WTs');
            zerwts = find(max_wts==0);

            % Assign a packet boundary to each point based on the maximized weights
            if isempty(zerwts)
                maxid_wts = maxid_wts.*10;
            else
                maxid_wts = (maxid_wts-1).*10;
            end

            % Set denoising weights to above packet-specific likelihoods
            DeNoiseWTs = WTs;
            
            % Create adjacency array based on specific interactions between 
            % neighboring nodes
            [grn_adj,grn_mori,Flg] = adjpt_moris(MartVars,1);
                       
            % If the maximum weight is 0 for all 4 packets, it's likely 
            % noise. Find the adjacent points with weights > 0 and reassign
            % a newly weighted value based on the mean weights of
            % neighboring points
            for jj = 1:length(zerwts)
                tmpid = zerwts(jj);
                [r,~]=find(grn_adj==tmpid);
                unqadj = unique(grn_adj(r,:));
                unqadj(unqadj==tmpid)=[];               
                DeNoiseWTs(tmpid,:) = mean(WTs(unqadj,:));
            end

            % Take the max of the denoised weights to use as IP weights
            [maxnew_wts,maxnew_wtsID] = max(DeNoiseWTs');
            
            % Pre-allocate variables
            len = length(AusGrn);
            Pack_Bounds = zeros(len,1);
            
            % Check to make sure the assigned packet is large enough to
            % establish an graph
            if Flg
                Pack_Pts = zeros(length(AusGrn),1);
            else
                % Index the  maxed weights and assign new IP weights as the
                % inverted difference between assigned boundary value for
                % each adjacent pair
                maxDN_new_wts1 = maxnew_wtsID(grn_adj(:,1));
                maxDN_new_wts2 = maxnew_wtsID(grn_adj(:,2));
                DN_grnIP_wts = IP./abs(maxDN_new_wts1-maxDN_new_wts2);
                DN_grnIP_wts(DN_grnIP_wts==Inf) = 1e2;

                % First set of OP weights are scaled to 10 
                DNc2g_wts = ones(length(DeNoiseWTs),1).*10;

                % Create graph consisting of nodes within the PAG
                [DeNoise_Graph,Enode,Snode] =...
                    graph_setup(grn_adj,DN_grnIP_wts,DNc2g_wts,2);

                % Perform the denoising cuts, which do utilize the graph
                % cutting technique to cluster variants into the most
                % likely packets
                for ii = 1:4
                    % Scaled original weights
                    DNg2t_wts = DeNoiseWTs(:,ii+1).*OP;
                    
                    % Ensure real positive weighted values
                    if any(isnan(DNg2t_wts))
                        DNg2t_wts(isnan(DNg2t_wts))=0;
                        DNg2t_wts(DNg2t_wts < 0) = 0;
                    end

                    % Add g2t weights to graph
                    DeNoise_Graph = ...
                        rmedge(DeNoise_Graph,(1:len)+Enode,Snode);
                    DeNoise_Graph = ...
                        addedge(DeNoise_Graph,(1:len)+Enode,Snode,DNg2t_wts);

                    % Initiate graph cutting
                    [mf_c,gf,cs,ct]=maxflow(DeNoise_Graph,1,Snode);
                    ct(end)=[];

                    % Record cut indices
                    DNgrn_ind = ct-Enode;

                    % Assign packet IDs and corresponding weights
                    Pack_Bounds(DNgrn_ind) = ii;
                    Pack_Pts{ii} = DNgrn_ind;
                    DN_PackWts{ii} = DNg2t_wts;
                end
            end
            
            % Establish final weights to be output to RunRecon function
            FinalWts = zeros(length(WTs),1);
            for i = 1:length(WTs)
                FinalWts(i) = WTs(i,Pack_Bounds(i)+1);
            end
                
    else
        % If grain was ID'd as unassigned martensite, establish weighted
        % values for consistency
        MartVars = [];
        Pack_Bounds = zeros(length(AusGrn),1);
        Pack_Pts = [];
        FinalWts = zeros(length(AusGrn),1);
        maxid_wts = zeros(length(AusGrn),1);
    end

        % Assign pertinent variables to output ``Packet'' structure
        Packet{gId}.MartGrains  = MartVars;
        Packet{gId}.AusGrain    = AusGrn;
        Packet{gId}.Boundaries  = Pack_Bounds;
        Packet{gId}.BoundaryIDs = Pack_Pts;
        Packet{gId}.Weights     = FinalWts;
        Packet{gId}.MaxBounds   = maxid_wts';
        Packet{gId}.Grain       = Grains.grains(gId);
        
        % Now, based on packet designations, segment variants within
        % packets to blocks (and sub-blocks)
        [Packet] = BlockChar(myEBSD,Grains,gId,Packet,Twin);
    
end

