function [EbAus,counter] = Auto_OR_Grn(myEBSD,counter)
% MTEX-driven, automated `reconstruction' algorithm used to segment a
% portion of a grain that can then be used for the automatic determination
% of the OR

% Initialization of requisite parameters
Ebsd = myEBSD.Ebsd;
CS_T = myEBSD.CS{2};
CS_R =  myEBSD.CS{3};
ksi  = [5.26,10.30,10.53]; % Begin with KS measurement
ksi = [2.8331    8.6910    8.9463]

% Compute variants and corresponding groupoid from euler angles
[V,flag] = YardleyVariants(ksi);
[G,CT] = GroupoidVariantsCubic(V);

% Since result is a passive rotation, take inverse of variant
% transformation indices to convert to active (Consistent with MTEX)
for i = 1:length(V)
    Vt{i,1} = transpose(V{i});
end

% Now transform the transposed variants to euler angles
for ii = 1:length(V)
    VTeul(ii) = orientation('matrix',Vt{ii},CS_T);
end

% Convert group axis to actual axis and then to Euler angles
Ax = vector3d(G(:,2),G(:,3),G(:,4));
G_eul = orientation('axis',Ax,'angle',G(:,1),CS_T);

%% Perform misorientation-based reconstruction

% Construct grain structure for martensitic variants
[grains,Ebsd.grainId] = calcGrains(Ebsd);

% Isolate the grain boundaries and adjust the phase to singular
gB = grains.boundary;
gB.phase=myEBSD.Phase.ID{1};

% Classify groupoid as a misorientation in MTEX terms
rot = orientation('Euler',[0,0,0],CS_T);
G_misos(1) = inv(rot) * rot;
G_misos(2:length(G_eul)+1) = inv(rot) * G_eul;
G_misos = transpose(G_misos);

% Set tolerance
if counter == 1 
    tol = 1*degree;
elseif counter == 2
    tol = 1.25;
elseif counter == 3
    tol = 0.75;
elseif counter > 3
    warning('Try Manual Segmentation of Potential PAG')
end

% Find the austenite grain boundaries based on the above tolerance
gB_aus = [];
for i = 1:length(G_misos)
    if i == 1
        gB_aus = gB(angle(gB.misorientation,G_misos(i)) < tol);
    else
        gB_aus = vertcat(gB_aus,gB(angle(gB.misorientation,G_misos(i)) < tol));
    end
end

% Merge grains together based on the chosen austenite grain boundaries
[mergedGrains,parentID]=merge(grains,gB_aus,'calcMeanorientation');
mergedGrains.CS = CS_R;

% Find the maximum `austenite' grain and extract the Ebsd coordinates
% related to the original dataset
[sz,maxGrnid] = max(mergedGrains.grainSize);
ParGrnId = find(parentID==maxGrnid);
EbId = [];
for i = 1:length(ParGrnId)
EbId = vertcat(EbId,find(Ebsd.grainId==ParGrnId(i)));
end

% Truncated ebsd dataset from `grain reconstruction'
EbGrn = Ebsd(EbId);
EbGrn.phase = myEBSD.Phase.ID{1};

% Compute misorientation distribution function based on KS-OR and
% corresponding groupoid
hw = 1.6*degree;
psi=deLaValeePoussinKernel('halfwidth',hw);
kappa=psi.kappa;
MartModf = calcODF(G_misos,'kernel',psi);

% Establish transformation matrices
[T2R,flag]=calc_T2R(ksi,CS_R,CS_T);
[R2T,flag]=calc_R2T(ksi,CS_R,CS_T);

% From our EBSD region, transform the martensite to austenite and find the
% modal austenite orientation to use as our guess
aus_odf=calcODF(symmetrise(EbGrn.orientations)*T2R,'kernel',psi);
[~,AusOr] = max(aus_odf);
%%
% Add source node
DiGraph=digraph;

% If this is the first iteration, initialize variables that only need to be
% established once
if counter == 1
    
    % Call function to compute adjacency array and corresponding
    % misorientations for the neighbgoring points
    [adjpts,mori,~] = adjpt_moris(EbGrn,1);

    % In-plane and out-of-plane weights
    modf_vals=eval(MartModf,mori);
    modf_vals(modf_vals<0)=0;
    IP = [3,5];
    OP = [2,5e-1];
    
else
    IP = IP + 3;
    OP = OP*2;
    OP1wts = OP1wts*2;
end
IPwts = (modf_vals+IP(1)).*IP(2);

% Compute ODF for parent-twin and establish the second set of out-of-plane 
% weights
From_OrGuess=symmetrise(AusOr)*R2T;
From_GuessODF=calcODF(From_OrGuess,'kernel',psi);
OP2wts=(eval(From_GuessODF,EbGrn.orientations)+OP(1)).*OP(2);

% Establish constant as first set of out-of-plane weights
OP1wts = mean(OP2wts)*ones(length(EbGrn),1);

% Call function to set up the initial graph
[DiGraph,endnode,sinknode] = graph_setup(DiGraph,adjpts,IPwts,OP1wts,2);

% Add edges corresponding to the guess reconstruction orientation    
DiGraph=rmedge(DiGraph,(1:length(EbGrn))+endnode,sinknode);
DiGraph=addedge(DiGraph,(1:length(EbGrn))+endnode,sinknode,OP2wts);

% Perform graph cut 
[mf_c,gf,cs,ct]=maxflow(DiGraph,1,sinknode);

% Truncate ct and finalize the cut
ct(end) = [];
Cut = ct - endnode;

end

