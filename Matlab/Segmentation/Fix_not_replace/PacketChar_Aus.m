function [Packet] = PacketChar_Aus(myEBSD,Grains,gId,Packet,PA_Grains)
%% Intialization
Aus_Grain = myEBSD.Recon.FullEbsd(myEBSD.Recon.FullEbsd.grainId == gId);
Mart_Grain = myEBSD.Ebsd(myEBSD.Recon.FullEbsd.grainId == gId);
AusOr = find_PAG_orientation(Aus_Grain);
CS_T        = myEBSD.CS{2};
ksi         = myEBSD.OR;
psi         = myEBSD.noise.psi;
IP = 5;
OP = 2e-1;
[VTeul,G_eul] = Var_orientations_from_ksi(ksi,CS_T);
%% Possibly Obsolete creating of tmpInds?
% Theoretically, tmpInds is unnecessary. However, in the off-case where
% there is a stray voxel in the PAG, this will either fix everything, or
% break everything. Later on, need to test then when run on every image.
% Might also never have a problem though? new routine doesn't seem to ever
% mess up boundaries like the old one.
tmpInds = find(Aus_Grain.orientations == AusOr);
martOrs = Mart_Grain(tmpInds).orientations;
len = length(tmpInds);

%% Build the graph
% AUSTIN NOTE: Ask Steve what G_eul is. I think this section is weighting
% neighbor bonds by likelyhood they are packet boundaries, not block, but
% unsure.
% Initialize the graph
grn_graph = digraph;
% Create adjacency array or misorientations between neighboring voxels
[grn_adj,grn_mori,~] = adjpt_moris(Mart_Grain(tmpInds),1);
% Create MODFs for block boundaries
mori_block = symmetrise(G_eul(1:4),CS_T); 
martblock_modf=calcODF(mori_block,'kernel',psi);
% Evaluate misorientations with respect to block-specific MODF
grn_mori.SS = myEBSD.SS;
martblock_modfvals = eval(martblock_modf,grn_mori);
martblock_modfvals(martblock_modfvals<0)=0;
% Reweight the In-plane connections (maybe should be function option?)
grainIP_wts = martblock_modfvals.*IP;


for iter = 1:4
    % ID for group of 6 variants being analyzed and the rest of them
    tmp1 = (iter-1)*6+1;
    tmp2 = (iter)*6;
    
    % Pull variants out
    Varsg2t = VTeul(tmp1:tmp2);
    Varsc2g = VTeul;
    Varsc2g(tmp1:tmp2) = [];
    
    % Rotate variant orientations by respective austenite orientation
    Var_rotg2t = rotation('Euler',AusOr.phi1,AusOr.Phi,AusOr.phi2,CS_T) * Varsg2t;
    Var_rotc2g = rotation('Euler',AusOr.phi1,AusOr.Phi,AusOr.phi2,CS_T) * Varsc2g;
    
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

WT_Len = max(max(tmpInds),size(Mart_Grain,1));
WTs = zeros(WT_Len,size(Weights,1));
WTs(tmpInds,:) = Weights';
%Austin Addition
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
    %        unqadj(unqadj > size(WTs,1))=[]; %Austin addition
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
Pack_Pts = [];
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

% Fill in an EBSD data set with the packet boundaries

% Assign pertinent variables to our ``Packet'' structure
Packet{gId}.MartGrains  = Mart_Grain;
Packet{gId}.AusGrain    = Aus_Grain;
Packet{gId}.Boundaries  = Pack_Bounds;
Packet{gId}.BoundaryIDs = Pack_Pts;
Packet{gId}.Weights     = FinalWts;
Packet{gId}.MaxBounds   = maxid_wts';
%Packet{gId}.Grain       = Grains.grains(gId);

[Packet] = BlockChar_Aus(myEBSD,Grains,gId,Packet);

end

function [VT_eul,G_eul] = Var_orientations_from_ksi(ksi,CS_T)
% Since result is a passive rotation, take transpose of variant
% transformation indices to convert to active (Consistent with MTEX)
% NOTE TOREN: Cryspy is all passive, so skip this transpose step when
% converting reconstruction code for cryspy. Then convert the
% transposed variants to euler angles.
[V,~] = YardleyVariants(ksi);
[G,~] = GroupoidVariantsCubic(V);
for i = 1:24
    Vt{i,1} = transpose(V{i});
    VT_eul(i) = orientation('matrix',Vt{i},CS_T);
end
% Convert group axis to actual axis and then to Euler angles
Ax = vector3d(G(:,2),G(:,3),G(:,4));
G_eul = orientation('axis',Ax,'angle',G(:,1),CS_T);

end


function AusOr = find_PAG_orientation(Aus_Grain)
%Find most common orientation in the PAG.
% in a perfect world, there should only ever be one option, and this
% function should be pointless, but occasionally there are errors in how
% the grain data is passed around, and explicitly asking for the mode helps
% make this code more resiliant to errors
[PAG_oris,~,pos] =unique(Aus_Grain.orientations);
PAG_ori_counts = accumarray(pos,1);
AusOr = PAG_oris(PAG_ori_counts==max(PAG_ori_counts));
end
