clear all
close all

filename = "AF_003.ang";
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

% Extract filename information
%filename = "AF_008.ang"
%filename = "5894-500-coarsened.ang"
IP =  [4,6];
myEBSD = [];
% If an OR is provided, set proper variable and continue on. The flag
% sips the OR determination in lieu of what the user provided.
OR_Flag = 1;
%myEBSD.OR = [3.09,8.10,8.48];
%myEBSD.noise.halfwidth = 1.7*degree;
%OR_Flag = 0;


%OR_Flag = 0;
%myEBSD.OR = [3.09,8.10,8.48];
%myEBSD.OR = [2.8772,9.1181,9.3199];
%myEBSD.noise.halfwidth =  2.1484*degree;



% Define the material (Steel only for now, will incorporate Titanium next)
myEBSD.Material = 'Steel';
myEBSD.Phase.Name{1} = 'Martensite';
myEBSD.Phase.Name{2} = 'Austenite';
myEBSD.Celldims{1} = [2.87, 2.87, 2.87];
myEBSD.Celldims{2} = [3.65, 3.65, 3.65];
myEBSD.rec_space = 'Mixed';
[myEBSD] = import_EBSD_Alex(filename,myEBSD);
clear filetype pathname
myEBSD.Ebsd.y = flipud(myEBSD.Ebsd.y);

%% 3.) OR DETERMINATION

% Use truncation function to crop an ebsd region inside of an assumed prior
% austenite grain based on input min/max coordinates. From plot above, we
% can manually select the region, note the coordinates, and then input here.
% Note that we don't need to select the entire PAG, just a large enough
% region to capture a sufficient number of martensite variants.

tic
plt_mart = 1;   % Flag to plot pdf martensite orientations
plt_vars = 1;   % Flag to plot pdf austenite variants from generated ksi values
gen=0;          % Generate synthetic dataset (0 means use existing)
mart_Ors=2000;  % Number of martensite points to consider when plotting
plt_ksi = 0;    % Flag to plot spread of computed ksi angles

% Call autoOR_estimation if flagged
if OR_Flag == 1
    % Adds computed ksi values and the loglikelihood of generated ksi values
    % (ksi{1} and ksi{2}, respectively) to myEBSD structure
    [myEBSD] = AutoOR_estimation(myEBSD,plt_mart,plt_vars,mart_Ors,plt_ksi);
    %%% assignin('base','myEBSD',myEBSD)
    %         [myEBSD] = OR_estimation(myEBSD,plt_mart,plt_vars,gen,mart_Ors,plt_ksi);
end
%% 4a.) Parameters for Reconstruction
% Note: Only need to run once in order to capture the following:
% transformation MODF, initial in-plane and out-of-plane weights (both are
% time consuming), adjacent points array, post-transformation phase ebsd
% dataset, and grid regions used for reconstruction parameters

% Returns an MODF for transformation microstructure and psi value to noise
% faction
if isfield(myEBSD,'MODF')==0
    [myEBSD] = calcMODF(myEBSD);
    %%% assignin('base','myEBSD',myEBSD)
end

% If dataset is deemed large enough, splits the dataset into 4
% quadrants to increase computational efficiency
if isfield(myEBSD,'TransMats')==0
    [myEBSD] = DataQuadrants(myEBSD);
    
    % Initialization parameters for reconstruction. Adds values to structure
    % for the transformation-phase ebsd (TransEbsd), adjacent point array
    % (Adj_pts), weights for both in-plane (Weights{1}) and out-of-plane
    % (first set; Weights{2}) and transformation values (orientations) allowing
    % for smooth transition from transformation phase to reconstruction phase
    [myEBSD] = initialize_recon(myEBSD);
    
    % Grid dimensions for transformation phase space that is used to either
    % produce austenite guess orientations (Austenite space; make large enough
    % to capture variant space) or randomly select martensite orientations to
    % choose from (martensite space: make smaller)
    
    if strcmp(myEBSD.rec_space,'Mixed')
        %             myEBSD.Or_guess = [];
        %             myEBSD.Like_Inds = [];
    else
        grid_x = 20;    % Spanned units in x direction
        grid_y = 18;    % Spanned units in y direction
        % grid_x = 32;
        % grid_y = 32;
        plot_grid = 1;  % Flag to plot grid over microstructure
        %
        % Produces a grid based on user-defined dimensions, resulting in an array
        % consisting of either post-transformation  or pre-transformation guess
        % orientations (Or_guesses) and corresponding smaller grid indices used in
        % "Reconstruction" to designate low-likelihood regions
        myEBSD = GridStructure(myEBSD,grid_x,grid_y,plot_grid);
    end
    %%% assignin('base','myEBSD',myEBSD)
end

%% 4b.) Perform Austenite Reconstruction and Plot

tic
% Regularization and scaling parameters for in-plane and out-of-plane weights
%     IP_reg = 4;         % Constant
IP_reg = IP(1);       % AFRL
%     IP_scale = 6;       % Constant
IP_scale = IP(2);     % AFRL
IP_params=[IP_reg,IP_scale];
% OP_reg = 2;       % Collin
OP_reg = 2;
%     OP_scale = 2e-1;  % Constant
OP_scale = 1.75e-1; % Diehl;
OP_params = [OP_reg,OP_scale];
% Number of iterations. Leave empty for default
iters = [];

% Perform Reconstruction in user-designated space. Returns added items to
% myEBSD structure corresponding to the total likelihood (Weights{3}),
% reconstructed ebsd (Recon), and the maximum likelihood reconstruction and
% corresponding weights (Max{1},{2}, respectively)
[myEBSD,Parent,Twin] = Call_Recon(myEBSD,IP_params,OP_params,iters);
recon_time = (toc/60)

%% 5.) Twin and Parent Grain Identification

if strcmp(myEBSD.rec_space,'Mixed')
else
    %     Angular tolerance for twins (in degrees)
    ang_tol = 2;
    
    %     Adds Parent and Twin structures to our workspace, which include Clean and
    %     Twinned Grains, both structures that contain the necessary indexing
    %     information and ID of each specific entity (or pairing)
    [myEBSD,Parent,Twin] = TwinParentID(myEBSD,ang_tol);
end

% CRAZY AUSTIN MOVE- Just going to skip this whole section and do my own grain definition
% 6.) Compute Pre-Transformation Grain Size and Plot Distribution

% Ignore Twin boundaries by setting to 1
Mergetwins = 1;

% Returns Grains structure which includes grain size, corresponding id,
% and ASTM information if requested based on the reconstruction space
if strcmp(myEBSD.rec_space,'Mixed')
    [Grains,Twin,Parent,myEBSD] = GrainSize_MixedSpace(myEBSD,Twin,Parent);
else
    [Grains] = GrainSize_SingSpace(myEBSD,Twin,Mergetwins);
end

% Save the 'Grain' structure to the EBSD workspace for the user
%%% assignin('base','AusGrains',Grains)

% 8.) Martensite Packet Boundary Characterization



for iiii = 1:4
    
    % Preallocate the necessary variables
    Packets         = [];
    PackRecon       = zeros(length(myEBSD.Ebsd),1);
    PackReconWts    = zeros(length(myEBSD.Ebsd),1);
    BlockReconWts   = zeros(length(myEBSD.Ebsd),1);
    BlockRecon      = zeros(length(myEBSD.Ebsd),1);
    VariantRecon    = zeros(length(myEBSD.Ebsd),1);
    VariantReconWts = zeros(length(myEBSD.Ebsd),1);
    
    
    
    %% Austin's fantastical Magical Segmentation Fun
    % Step 0: Erase evidence of last attempt
    Packets         = [];
    PackRecon       = zeros(length(myEBSD.Ebsd),1);
    PackReconWts    = zeros(length(myEBSD.Ebsd),1);
    BlockReconWts   = zeros(length(myEBSD.Ebsd),1);
    BlockRecon      = zeros(length(myEBSD.Ebsd),1);
    VariantRecon    = zeros(length(myEBSD.Ebsd),1);
    VariantReconWts = zeros(length(myEBSD.Ebsd),1);
    
    %Step 1: do a grain segmentation and assign threshold grain size
    % quick aside to future me, this does NOT delete the small grains, it just
    % rearranges the ordering to grains below the threshold get a new grain ID
    % that is AFTER the threshold "num_grains" and no longer have a polygon
    % related to them in the "PA_grains" object.
    [PA_Grains,myEBSD.Recon.FullEbsd.grainId] = calcGrains(myEBSD.Recon.FullEbsd,'angle',1*degree);
    [myEBSD.Recon.FullEbsd,PA_Grains,~] = RemSmGrns(myEBSD.Recon.FullEbsd,PA_Grains,100);
    num_grains = size(PA_Grains,1);
    %%
    
    %Step 2: loop over grains above threshold size
    
    figure()
    plot(myEBSD.Recon.FullEbsd("Austenite"),myEBSD.Recon.FullEbsd("Austenite").grainId)
    caxis([0 num_grains])
    
    for k = 1:length(PA_Grains)
        disp(['now working on Grain ' int2str(k) ' of ' int2str(length(PA_Grains))])
        try
            % Now determine packets for the characterized austenite grains (ignores
            % the unassigned martensite)
            [Packets] = PacketChar_Aus(myEBSD,k,Packets);
            % The actual austenite id number
            AusId = Packets{k}.AusGrain.id;
            % For each pixel in the EBSD microstructure, fill with the
            % corresponding packet, block, and variant assignments
            PackRecon(AusId)       = Packets{k}.Boundaries;
            PackReconWts(AusId)    = Packets{k}.Weights;
            BlockRecon(AusId)      = Packets{k}.Block.Boundaries;
            BlockReconWts(AusId)   = Packets{k}.Block.Weights;
            VariantRecon(AusId)    = Packets{k}.Variants.List;
            VariantReconWts(AusId) = Packets{k}.Variants.Weights;
        catch
            disp(['Error: on grain ' int2str(k) ' The code failed to properly segment the'])
            disp('packets. This is likely due to a twinning effect.')
        end
    end
    
    
    % Assign these values to the myEBSD structure and then send to the
    % workspace base
    myEBSD.Packets.IDs = PackRecon;
    myEBSD.Packets.Wts = PackReconWts;
    myEBSD.Packets.RGB = [0,0,0;
        1,0,0;
        0,1,0;
        0,0,1;
        1,1,0];
    myEBSD.Blocks.IDs = BlockRecon;
    myEBSD.Blocks.Wts = BlockReconWts;
    myEBSD.Blocks.RGB = [0,0,0;
        252,122,120;
        250,4,0;
        102,0,0;
        120,252,122;
        4,250,0;
        0,125,2;
        130,125,253;
        2,2,250;
        1,1,100;
        220,220,170;
        250,250,2;
        227,176,10]./255;
    myEBSD.Variants.IDs = VariantRecon;
    myEBSD.Variants.Wts = VariantReconWts;
    
    % Assign the sub-block RGB values as well
    SubBlockRGB = zeros(25,3);
    RGBcond = [-1*50/255; 50/255];
    counter = 2;
    for mm = 2:13
        tmpRGB = myEBSD.Blocks.RGB(mm,:);
        CondRGB = tmpRGB > RGBcond(2);
        BlckRGB = tmpRGB;
        for nn = 1:2
            tmpRGB(CondRGB) = tmpRGB(CondRGB) + RGBcond(nn)/2;
            tmpRGB(tmpRGB > 1) = 1;
            SubBlckRGB(counter,:) = tmpRGB;
            counter = counter+1;
            tmpRGB = BlckRGB;
        end
    end
    
    % Add RGB values to subblock designations
    myEBSD.Variants.RGB = SubBlckRGB;
    %%% assignin('base','AusGrnPackets',Packets);
    
    %
    myEBSD = CharMartBoundaries(myEBSD);
    %%% assignin('base','myEBSD',myEBSD);
    %
    % Now write pertinent data to the text file
    % genTxt(myEBSD,Twin,Grains,filename)
    myEBSD.TwinGrains = Twin;
    myEBSD.ParentGrains = Parent;
    myEBSD.AusGrains = Grains;
    myEBSD.AusGrnPackets = Packets;
    myEBSD.filename = filename;
    assignin('base','myEBSD',myEBSD);
    
    
    [path,name,~] = fileparts(filename);
%    mat_name = [path '/'  name '_1_' int2str(iiii) '_Recon.mat'];
    mat_name = [ name '_3_' int2str(iiii) '_Recon.mat'];
    mat_name = strjoin(mat_name, '');
    %save(mat_name, 'myEBSD', '-v7.3')
    save('dummy', 'myEBSD', '-v7.3')
    clear path name
    
    PltSubBlocks(myEBSD,1)
end

