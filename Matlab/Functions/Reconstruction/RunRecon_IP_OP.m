function RunRecon_IP_OP(myEBSD,filename,IP,OP,angular_tolerance)
% Home script with all requisite functions related to:
% 1.) Loading ebsd dataset
% 2.) Assistance in truncating ebsd dataset to specific region and plotting
%     of corresponding IPF map
% 3.) Assistance in truncating it to a region within a single
%     prior austenite grain and utilizing said region to evaluate the
%     steel's orientation relationship
% 4.) Reconstruction of the pre-transformation austenite microstructure for
%     either a portion of the microstructure or the entire one
% 5.) Identification of all computed parent-twin systems and standalone
%     parent austenite grains
% 6.) Plotting histogram for the prior austenite grain size and average
%     prior austenite grain size
% 7.) Adjustment of twin boundaries for most realistic representation
% 8.) Characterization of martensite packet boundaries within single
%     austenite grains

% Set up mtex plotting preference
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% 4b.) Perform Austenite Reconstruction and Plot

tic
% Number of iterations. Leave empty for default
iters = [];
% Perform Reconstruction in user-designated space. Returns added items to
% myEBSD structure corresponding to the total likelihood (Weights{3}),
% reconstructed ebsd (Recon), and the maximum likelihood reconstruction and
% corresponding weights (Max{1},{2}, respectively)
[myEBSD,Parent,Twin] = Call_Recon(myEBSD,IP,OP,iters);
recon_time = (toc/60)
AusRecon = myEBSD.Recon.Ebsd;
Likelihood = myEBSD.Recon.Likelihood;

%% 5.) Twin and Parent Grain Identification

if strcmp(myEBSD.rec_space,'Mixed')
else
    %     Adds Parent and Twin structures to our workspace, which include Clean and
    %     Twinned Grains, both structures that contain the necessary indexing
    %     information and ID of each specific entity (or pairing)
    [myEBSD,Parent,Twin] = TwinParentID(myEBSD,angular_tolerance);
end
%%
% 6.) Compute Pre-Transformation Grain Size and Plot Distribution

% Ignore Twin boundaries by setting to 1
Mergetwins = 1;

% Returns Grains structure which includes grain size, corresponding id,
% and ASTM information if requested based on the reconstruction space
[Grains,Twin,Parent,myEBSD] = GrainSize_MixedSpace(myEBSD,Twin,Parent);

% Save the 'Grain' structure to the EBSD workspace for the user
%%% assignin('base','AusGrains',Grains)

%% 8.) Martensite Packet Boundary Characterization

% Preallocate the necessary variables
Packets         = [];
PackRecon       = zeros(length(myEBSD.Ebsd),1);
PackReconWts    = zeros(length(myEBSD.Ebsd),1);
BlockReconWts   = zeros(length(myEBSD.Ebsd),1);
BlockRecon      = zeros(length(myEBSD.Ebsd),1);
VariantRecon    = zeros(length(myEBSD.Ebsd),1);
VariantReconWts = zeros(length(myEBSD.Ebsd),1);

% Loop through each PAG and assign a packet, block, and variant ID to
% each pixel within the PAG
for k = 1:length(Grains.grainId)
    disp(['now working on Grain ' int2str(k) ' of ' int2str(length(Grains.grainId))])
    try
        % Now determine packets for the characterized austenite grains (ignores
        % the unassigned martensite)
        [Packets] = PacketChar(myEBSD,Grains,k,Packets,Twin);
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

%%
myEBSD = CharMartBoundaries(myEBSD);
%%% assignin('base','myEBSD',myEBSD);
%%
% Now write pertinent data to the text file
% genTxt(myEBSD,Twin,Grains,filename)
myEBSD.AusRecon_Ebsd = AusRecon;
myEBSD.AusRecon_Likelihood = Likelihood;
myEBSD.TwinGrains = Twin;
myEBSD.ParentGrains = Parent;
myEBSD.AusGrains = Grains;
myEBSD.AusGrnPackets = Packets;
myEBSD.filename = filename;
assignin('base','myEBSD',myEBSD);
end
