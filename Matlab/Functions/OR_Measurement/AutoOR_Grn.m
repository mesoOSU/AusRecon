function [myEBSD,EbAus,EbGrn,AusOr] = AutoOR_Grn(varargin)
% MTEX-driven, automated `reconstruction' algorithm used to segment a
% portion of a grain that can then be used for the automatic determination
% of the OR

myEBSD = varargin{1};
EbGrn = varargin{2};
counter = varargin{3};

ksiKS = [5.26,10.30,10.53];
ksiNW = [0,9.74,9.74];

% If the user inputs a specific starting OR, use as the base for the
% automatic determination of the OR. If not, default with the GT OR.
if length(varargin)==4
    ksi = varargin{4};
else
%     ksi  = [5.26,10.30,10.53]; % Begin with KS measurement
    ksi = [2.16,8.06,8.30];     % Begin with GT measurement
end

% Initialization of requisite parameters
CS_T = myEBSD.CS{2};
CS_R =  myEBSD.CS{3};
    
% Compute variants and corresponding groupoid from euler angles
[V,~] = YardleyVariants(ksi);
[G,~] = GroupoidVariantsCubic(V);

% Convert group axis to actual axis and then to Euler angles
Ax = vector3d(G(:,2),G(:,3),G(:,4));
G_eul = orientation('axis',Ax,'angle',G(:,1),CS_T);

% Classify groupoid as a misorientation in MTEX terms
rot = orientation('Euler',[0,0,0],CS_T);
G_misos(1) = inv(rot) * rot;
G_misos(2:length(G_eul)+1) = inv(rot) * G_eul;
G_misos = transpose(G_misos);
lenEb = length(EbGrn);

% Set tolerance
% if counter == 1 
if (all(ksi==ksiKS) || all(ksi==ksiNW)) 
    tol = 1.25*degree;
else
    tol = 0.85*degree;
end

% For the first iteration, perform misorientation-based reconstruction to
% capture a conglomerate of grains
if counter == 1 || lenEb < 5e2
    
    if counter > 1
        tmpEbsd = myEBSD.Ebsd;
        tol = ((tol/degree) + 0.25)*degree;
    else
        tmpEbsd = EbGrn;
    end

    % Construct grain structure for martensitic variants
    [grains,tmpEbsd.grainId] = calcGrains(tmpEbsd);

    % Isolate the grain boundaries and adjust the phase to singular
    gB = grains.boundary;
    gB.phase=myEBSD.Phase.ID{1};

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
        EbId = vertcat(EbId,find(tmpEbsd.grainId==ParGrnId(i)));
    end

    % Truncated ebsd dataset from `grain reconstruction'
    EbGrn = tmpEbsd(EbId);
    EbGrn.phase = myEBSD.Phase.ID{1};
end

% Compute misorientation distribution function based on KS-OR and
% corresponding groupoid
hw = 2*degree;
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

% If this is the first iteration, initialize variables that only need to be
% established once
% if counter == 1

% Call function to compute adjacency array and corresponding
% misorientations for the neighbgoring points
[adjpts,mori,~] = adjpt_moris(EbGrn,1);

% In-plane and out-of-plane weights
modf_vals=eval(MartModf,mori);
modf_vals(modf_vals<0)=0;

% If the input OR is KS or NW, adjust the parameterization for these. If
% it's GT or something similar, alter as well and include for twins.
if (all(ksi==ksiKS) || all(ksi==ksiNW)) 
    IP = [10,15];
    OP = [2,1.5e-1];
    TwnWts = zeros(length(EbGrn),1);
else
    IP = [7,8];
    OP = [2,1.5e-1];
    Twns = vector3d([1,1,1;-1,-1,1;-1,1,1;1,-1,1]');
    Twn_eul = orientation('axis',Twns,'angle',60*degree,CS_R);
    AusOrRot = rotation('euler',AusOr.phi1,AusOr.Phi,AusOr.phi2);
    TwnAus = AusOrRot * Twn_eul;
    TwnOrs = symmetrise(TwnAus)*R2T;
    TwnODF = calcODF(TwnOrs,'kernel',psi);
    TwnWts = (eval(TwnODF,EbGrn.orientations)+OP(1)).*OP(2);
end

% Compute ODF for parent and establish the second set of out-of-plane 
% weights
From_OrGuess=symmetrise(AusOr)*R2T;
From_GuessODF=calcODF(From_OrGuess,'kernel',psi);
OP2wts=(eval(From_GuessODF,EbGrn.orientations)+OP(1)).*OP(2);

% Establish constant as first set of out-of-plane weights
OP1wts = abs(1./OP2wts);
OP1wts(OP1wts==Inf)=1e2;
OP1wts = OP1wts + TwnWts;
    
IPwts = (modf_vals+IP(1)).*IP(2);

% Add source node
DiGraph=digraph;

% Establish constant as first set of out-of-plane weights
OP1wts = abs(4./OP2wts);
OP1wts(OP1wts==Inf)=1e2;

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

% Send back the austenite cut 
EbAus = EbGrn(Cut);
EbGrn(Cut) = [];

end

