function [Packet] = BlockChar(myEBSD,Grains,gId,Packet)
% For each PAG (and underlying twin), identify the corresponding
% block boundaries for each respective packet

Aus_Grain   = myEBSD.Recon.Ebsd(Grains.Indices{gId});
MartRecon   = myEBSD.Ebsd;
CS_T        = myEBSD.CS{2};
TransID     = myEBSD.Phase.ID{1};
ReconID     = myEBSD.Phase.ID{2};
SS          = myEBSD.SS;
ksi         = myEBSD.OR;
psi         = myEBSD.noise.psi;

Ptmp = Packet{gId};
% Zero array for each pixel to be assigned a block and variant ID
BlockWts = zeros(length(Ptmp.AusGrain),1);
BlockBounds = zeros(length(Ptmp.AusGrain),1);
VarWts = zeros(length(Ptmp.AusGrain),1);
VarLst = zeros(length(Ptmp.AusGrain),1);

% If no austenite orientation was assigned to the distribution of
% martensite variants, don't assign block boundaries or a variant
% labeling to them.
if Grains.grainId(gId,2) == TransID
    Packet{gId}.Block.Boundaries = zeros(length(Ptmp.AusGrain),1);
    Packet{gId}.Block.Weights        = zeros(length(Ptmp.AusGrain),1);
    Packet{gId}.Variants.List     = zeros(length(Ptmp.AusGrain),1);
    Packet{gId}.Variants.Weights     = zeros(length(Ptmp.AusGrain),1);
else
    
    % Interface parameters (constant for now, see if these do change
    % much with variable steel and ferrous alloy data sets)
    IP = 5;
    OP = 2e-1;
    
    % Make sure martensite grain doesn't have residual austenite
    AusIDs            = find(MartRecon.phase==ReconID);
    MartRecon(AusIDs) = [];
    Mart_Grain        = Ptmp.MartGrains;
    
    [V,~] = YardleyVariants(ksi);
    [Gtmp,~] = GroupoidVariantsCubic(V);
    AusOr = Aus_Grain.orientations(1);
    % Extract the parent or twin (if it exists) indices
    tmpInds = find(Aus_Grain.orientations == AusOr);
    % Make the chosen austenite orientation a rotation
    Aus_or = rotation('Euler',AusOr.phi1,AusOr.Phi,AusOr.phi2,CS_T);
    % Compute variants and corresponding groupoid from euler angles
    % Assign non-misorientation since this will play a role down
    % the line
    G0 = orientation('euler',[0,0,0],CS_T);
    
    
    % Convert group axis to actual axis and then to Euler angles
    Ax = vector3d(Gtmp(:,2),Gtmp(:,3),Gtmp(:,4));
    GEul = orientation('axis',Ax,'angle',Gtmp(:,1),CS_T);
    % Kappa value
    kappa=psi.kappa;
    
    G_eul = [G0;GEul];
    
    % Compute MODF for only sub-block boundaries
    mori_subblck = symmetrise(G_eul(1:2),CS_T);
    mori_rest = symmetrise(G_eul(3:17),CS_T);
    mart_subblck_modf=calcODF(mori_subblck,'kernel',psi);
    martrest_modf = calcODF(mori_rest,'kernel',psi);
    
    % Now we loop through blocks and variants for assignments regarding
    % each parent and twin within the same system. ii denotes parent or
    % twin; Pck denotes specific packet; Blck denotes specific block
    % within respective packet; mm denotes variant within a block
    ii = 1;
    % Extract relevant variants related to a single PAG/PAT
    Vst = 24*(ii-1)+1;
    Vnd = 24*ii;
    for m = Vst:Vnd
        Vt{m,1} = transpose(V{m});
        VTeul(m) = orientation('matrix',Vt{m},CS_T);
    end
    
    % Parent-twin variant Id (ease of labeling later on)
    PTVar = 0;
    
    % Zero array for each block in the PAG
    PckBlckWts = [];
    PckBlckInds = [];
    VariantWts = [];
    VariantInds = [];
    % Extract only martensite indices that are related to a
    % specific parent OR twin, not the entire system (if a twin
    % exists within a parent)
    Mart_tmp = Mart_Grain(tmpInds);
    
    % Now extract packet-specific variants from the given variant
    % list to be used cycled through to find the block assignments
    for Pck = 1:4
        PckGraph = digraph;
        PckInds = find(Ptmp.Boundaries(tmpInds)==Pck);
        len = length(PckInds);
        
        % Create adjacency array based on specific interactions between points
        % within the input grain
        [PckAdj,PckMori,Flg] = adjpt_moris(Mart_tmp(PckInds),1);
        if Flg==0
            % Extract the martensite packet orientations
            PckOrs = Mart_tmp(PckInds).orientations;
            % Assign a sample symmetry to the grain
            PckMori.SS = SS;
            
            % Evaluate misorientations with respect to block-specific values
            mart_subblck_modfvals = eval(mart_subblck_modf,PckMori);
            mart_subblck_modfvals(mart_subblck_modfvals<0)=0;
            % Now set up our graph
            PckIP_wts = mart_subblck_modfvals.*IP;
            
            BlckWts = zeros(len,3);
            % Now (finally) loop through the specific blocks and
            % determine which martensite variants belong to which
            % block within the packet
            for Blck = 1:3
                % ID for group of 6 variants being analyzed within
                % a specific packet (jj)
                Vtmp1 = 6*(Pck-1)+((Blck-1)*2+1);
                Vtmp2 = 6*(Pck-1)+((Blck)*2);
                
                % Pull variants out
                Blck_g2t = VTeul(Vtmp1:Vtmp2);
                Blck_c2g = VTeul;
                Blck_c2g(Vtmp1:Vtmp2) = [];
                
                % Rotate variant orientations by respective austenite orientation
                Blck_rotg2t = Aus_or * Blck_g2t;
                
                % Compute ODF for both sets of variants
                Blck_g2tODF = calcODF(symmetrise(Blck_rotg2t),'kernel',psi);
                %                         Blck_c2gODF = calcODF(Blck_rotc2g,'kernel',psi);
                %
                %                         % Compute the c2g weights (first set of OP weights)
                %                         c2g_wts = eval(Blck_c2gODF,PckOrs);
                %                         c2g_wts(find(c2g_wts<0))=0;
                
                % Compute the g2t weights (second set of OP weights)
                g2t_wts = eval(Blck_g2tODF,PckOrs);
                g2t_wts(find(g2t_wts<0))=0;
                %
                %                         % Only if we're on the first iteration should graph be
                %                         % constructed
                %                         if Blck==1
                %                             % Setup graph within the single austenite grain boundaries
                %                             [PckGraph,PckEndnode,PckSinknode] = graph_setup(PckGraph,PckAdj,PckIP_wts,c2g_wts,2);
                %                         end
                
                %                         % Add g2t weights to graph
                %                         PckGraph=rmedge(PckGraph,(1:len)+PckEndnode,PckSinknode);
                %                         PckGraph=addedge(PckGraph,(1:len)+PckEndnode,PckSinknode,g2t_wts);
                
                %                         % Commence graph cutting
                %                         [mf_c,gf,cs,ct]=maxflow(PckGraph,1,PckSinknode);
                %                         ctcopy = ct;
                %                         ctcopy(end)=[];
                %                         BlckInd = ctcopy-PckEndnode;
                
                % Assign cut out block boundaries, indices and corresponding
                % likelihoods
                BlckWts(:,Blck) = g2t_wts;
                VWts = zeros(len,2);
                for Var = 1:2
                    
                    Vi = Blck_rotg2t(Var);
                    Vi_g2tODF = calcODF(symmetrise(Vi),'kernel',psi);
                    
                    % Compute the g2t weights (second set of OP weights)
                    g2ti_wts = eval(Vi_g2tODF,PckOrs);
                    g2ti_wts(find(g2t_wts<0))=0;
                    
                    VWts(:,Var) = g2ti_wts;
                end % Loop through variants
                VariantWts{Blck} = VWts;
                %                           % Remove necessary c2t weights and add again
                %                           PckGraph=rmedge(PckGraph,(1:len)+1,(1:len)+PckEndnode);
                %                           PckGraph=addedge(PckGraph,(1:len)+1,(1:len)+PckEndnode,c2g_wts);
            end % Loop through blocks
        else
            BlckWts = [];
        end
        % Conditional statement based on size of Packet
        PckBlckWts{Pck} = BlckWts;
        PckBlckInds{Pck} = tmpInds(PckInds);
        PckVarWts{Pck} = VariantWts;
    end % Packet Loop
    
    % Loop back through the packets and relate the block weights to the
    % most probabilistic variant orientation for each respective packet
    for jj = 1:4
        if isempty(PckBlckWts{jj})==0
            FinInds = PckBlckInds{jj};
            % Identify max weights
            [max_Blckwts,BlckIds] = max(PckBlckWts{jj}');
            BlockWts(FinInds) = max_Blckwts;
            BlockBounds(FinInds) = (jj-1)*3+BlckIds;
            
            % From the chosen block, loop through each index and find
            % the max variant weight within said block
            for kk = 1:length(FinInds)
                BlckChoice = (jj-1)*3+BlckIds(kk);
                VarChoice = PckVarWts{jj}{BlckIds(kk)};
                [max_Varwts,VarId] = max(VarChoice(kk,:));
                VarWts(FinInds(kk)) = max_Varwts;
                VarLst(FinInds(kk)) = PTVar+2*(BlckChoice-1)+VarId;
            end
        end
    end
 
    % Now add the block and pixels to the packet structure and bring
    % this baby back home
    Packet{gId}.Block.Boundaries = BlockBounds;
    Packet{gId}.Block.Weights    = BlockWts;
    Packet{gId}.Variants.List    = VarLst;
    Packet{gId}.Variants.Weights = VarWts;
    
end % Conditional statement based on martensite being assigned a g_aus

end

