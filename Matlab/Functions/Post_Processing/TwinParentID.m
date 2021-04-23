function [myEBSD,Parent,Twin] = TwinParentID(myEBSD,ang_tol,varargin)
    % Identification of twins within transformed austenite microstructure.
    % Output is a collection of EBSD datasets that include both parents and
    % parent-twin systems

    Recon_Ebsd = myEBSD.Recon.Ebsd;
    Trans_Ebsd = myEBSD.TransEbsd;
    CS_R = myEBSD.CS{3};
    Recon_pts = myEBSD.Recon.TransformedPoints;
    TransID = myEBSD.Phase.ID{1};
    ReconID = myEBSD.Phase.ID{2};

    if isempty(varargin)
        Parent = [];
        Twin = [];
    else
        Parent = varargin{1};
        Twin = varargin{2};
    end

    %% First Classification of Transformed Grains
    
    if strcmp(myEBSD.rec_space,'Mixed')
        PhaseIDs = myEBSD.Recon.PhaseIDs;
        TransPts = find(PhaseIDs==TransID);
        
        % Cluster martensite grains as non-oriented martensite points
        Recon_Ebsd(TransPts).phase = ReconID;
        ReconOrs = Parent.Or;
        ReconOrs_ind = Parent.Indices;
        for ii = 1:length(Twin.Merged)
            tmpOr{1,1} = Twin.Parent.Or{ii};
            tmpInd{1,1} = Twin.Parent.Indices{ii};
            tmpOr{2,1} = Twin.Single.Or{ii};
            tmpInd{2,1} = Twin.Single.Indices{ii};
            count = 3;
            for jj = 2:4
                if isempty(Twin.Single.Or{ii,jj})==0
                    tmpOr{count} = Twin.Single.Or{ii,jj};
                    tmpInd{count} = Twin.Single.Indices{ii,jj};
                    count = count+1;
                end
            end
            ReconOrs = [ReconOrs;tmpOr];
            ReconOrs_ind = [ReconOrs_ind;tmpInd];
            clear tmpOr tmpInd
        end
        
        len_ReconOrs = zeros(length(ReconOrs),1);
        for i = 1:length(ReconOrs)
            Recon_Ors(i) = ReconOrs{i};
            len_ReconOrs(i) = length(ReconOrs_ind{i});
        end
    else
        
        % Reassign ebsd datasets based on only transformed phase grains
        Recon_Ebsd = Recon_Ebsd(Recon_pts);
        Trans_Ebsd = Trans_Ebsd(Recon_pts);

        Recon_Ors = Recon_Ebsd.orientations; % Index the transformed orientations
        Recon_Ors = unique(Recon_Ors);                  % Classify the unique ones
        ReconOrs_ind = cell(length(Recon_Ors),1);
        len_ReconOrs = zeros(length(Recon_Ors),1);

        % Pull out indices for each indexed grain
        for i = 1:length(Recon_Ors)
            tmpind = find(Recon_Ebsd.orientations==Recon_Ors(i));
            ReconOrs_ind{i,1} = tmpind;
            len_ReconOrs(i) = length(tmpind);  
        end

        % Delete orientations of grains that are extremely small (likely noise, 
        % ergo won't contain a twin)
        deleteOrs = find(len_ReconOrs<25);
        len_ReconOrs(deleteOrs) = [];
        ReconOrs_ind(deleteOrs)=[];
        Recon_Ors(deleteOrs)=[];
    end

    % Now sort indexed austenite grains from largest to smallest
    [len_ReconOrs,Recon_hi2lo] = sort(len_ReconOrs,'descend');
    Recon_Ors = Recon_Ors(Recon_hi2lo);
    ReconOrs_ind = ReconOrs_ind(Recon_hi2lo);
    
    clear Parent Twin

    %% Capture All Twins Based on Misorientation Angle Information

    % Find min/max coordinates for aus_ebsd dataset
    min_x = min(Recon_Ebsd.x);
    max_x = max(Recon_Ebsd.x);
    min_y = min(Recon_Ebsd.y);
    max_y = max(Recon_Ebsd.y);

    % Determine symmetry
    PG = LaueName(Recon_Ebsd.CS);

    % If cubic m-3m
    if strcmp(PG,'m-3m')
        max_ang = 60;
    else
    % Add for hcp (Ti)
    end

    % Now we set up a loop that will effectively capture twins that are
    % adjacent to or nestled within a specific grain
    cub_twins =[max_ang-ang_tol,max_ang+ang_tol];
    twinIDs = [];
    twin_Ind = [];
%     Par_Ind = [];
%     Partwin_Or = [];
%     ParTwin = [];
%     ParTwinID = [];
%     twinpar_Ind = [];
    num_twins = 0;
    notwn_count = 1;
    twin_par = 1;
    counter = 0;
    for active_grns = 1:length(Recon_Ors)
        if any(active_grns==twinIDs)
            % Skip grains already identified as twins
        else
            % Create temporary ebsd dataset
            tmpebsd = Recon_Ebsd;
            % Active grain orientation
            gi = Recon_Ors(active_grns);

            % Local min/max values for each grain increased by some user-defined
            % coordinate tolerance
            gimin_x = min(tmpebsd(ReconOrs_ind{active_grns}).x); % Min-x
            if (gimin_x-ang_tol) > min_x
                gimin_x = gimin_x-ang_tol;
            else
                boundary = min_x-gimin_x;
                gimin_x = gimin_x-boundary;
            end
            gimax_x = max(tmpebsd(ReconOrs_ind{active_grns}).x)+ang_tol; % Max-x
            if (gimax_x+ang_tol) < max_x
                gimax_x = gimax_x+ang_tol;
            else
                boundary = max_x-gimax_x;
                gimax_x = gimax_x+boundary;
            end
            gimin_y = min(tmpebsd(ReconOrs_ind{active_grns}).y)+ang_tol; % Min-y
            if (gimin_y-ang_tol) > min_y
                gimin_y = gimin_y-ang_tol;
            else
                boundary = min_y-gimin_y;
                gimin_y = gimin_y-boundary;
            end
            gimax_y = max(tmpebsd(ReconOrs_ind{active_grns}).y)+ang_tol; % Max-y
            if (gimax_y+ang_tol) < max_y
                gimax_y = gimax_y+ang_tol;
            else
                boundary = max_x-gimax_x;
                gimax_y = gimax_y+boundary;
            end

            % Delete points outside of indexed grain plus tolerance
            tmpebsd(find(tmpebsd.x < gimin_x | tmpebsd.x > gimax_x | tmpebsd.y < gimin_y | tmpebsd.y > gimax_y))=[];
            % Unique orientations within tmpebsd
            Gi = unique(tmpebsd.orientations);

            % Misorientation angles within tmpebsd orientations
            mori_gi = inv(gi)*Gi; 
            mori_ang = (angle(mori_gi)./degree)';
            mori_ax= axis(mori_gi)';

            % Of all computed misorientations for the given parent orientation,
            % find the corresponding twin elements
            find_twins = find(mori_ang > cub_twins(1) & mori_ang < cub_twins(2));
            twin_plane = Miller(vector3d(1,1,1),CS_R);

            % Now make sure misorientation axes are close (within a user-defined
            % error) to the fcc twin plane {111}
            twin_axcomp = angle(mori_ax(find_twins),twin_plane)./degree;
            find_twins(twin_axcomp>ang_tol)=[];

            % Finally, compare twin orientations to see if we should combine
            % two separate twins as a single orientation
            len_twin = length(find_twins);
            j=1;
            while j < len_twin
                k = j+1;
                twin_angcomp = angle(Gi(find_twins(j)),Gi(find_twins(k)))./degree;
                if twin_angcomp < 1
                    tj = Gi(find_twins(j));
                    tk = Gi(find_twins(k));
                    tj_ind = find(Recon_Ebsd.orientations==tj);
                    tk_ind = find(Recon_Ebsd.orientations==tk);
                    lenj = length(tj_ind);
                    lenk = length(tk_ind);
                    if lenj > lenk
                        Recon_Ebsd(tk_ind).orientations=tj;
                        find_twins(k)=[];
                        twinIDs = vertcat(twinIDs,find(Recon_Ors==tk));
                    else
                        Recon_Ebsd(tj_ind).orientations=tk;
                        find_twins(j)=[];
                        twinIDs = vertcat(twinIDs,find(Recon_Ors==tj));
                    end
                    len_twin = length(find_twins);
                    clear len1 len2 tj tk tj_ind tk_ind lenj lenk
                else
                    j = j+1;
                end
            end

            % Twin orientations
            twin_Ors = Gi(find_twins);

            % Ensure that we can find twins within austenite ebsd or just 
            % discard as noise
            if j==1 && isempty(twin_Ors)==0
                lentwin = find(Recon_Ors==twin_Ors);
                if isempty(lentwin)==1
                    find_twins=[];
                end
            end
            clear j k twin_angcomp 

            % Find index of parent grain and corresponding twins, starting with
            % those from the parent grain
            counter = counter+1;
            par_ind = ReconOrs_ind{active_grns};

            if isempty(find_twins) 
                % No twins exist for parent grain
                Par_Ind{notwn_count,1} = par_ind;
                Parcount{notwn_count,1} = counter; 
                notwn_count = notwn_count+1;
                
            else
                pre_count = counter;
                ParTwin{twin_par,1} = par_ind;
                ParTwinID{twin_par,1} = counter;
                twin_par = twin_par+1;
                gi_twins = zeros(len_twin,1);
                
                % Inner loop to add indices to what we'll graph.
                for i = 1:len_twin
                    num_twins = num_twins+1;
                    
                    counter = counter+1;
                    % Increase counter on number of twins in dataset
                    twinInds = find(Recon_Ebsd.orientations==twin_Ors(i));
                    partwin_ind = vertcat(par_ind,twinInds);
                    
                    gi_twins(i,1) = find(Recon_Ors==twin_Ors(i));

                    % Save indices and orientations for each twin and parent
                    % within the same system
                    Partwin_Or{num_twins,1} = gi;
                    Twinonly_Or{num_twins,1} = twin_Ors(i);
                    twinpar_Ind{num_twins,1} = partwin_ind;  % Parent-Twin system
                    twinparind{num_twins,1} = par_ind;
                    twinparID{num_twins,1} = pre_count;      % Parent Counter
                    twin_Ind{num_twins,1} = twinInds;        % Twin Indices
                    twinID{num_twins,1} = counter;           % Twin Only

                    clear twin_ind partwin_ind
                end
                % Now find indexed twins and don't run through iteration with them
                % (eliminate redundant conclusions)
                twinIDs = vertcat(twinIDs,gi_twins);
            end
            clear tmpebsd gi gimin_x gimax_x gimin_y gimax_y boundary Gi...
                  morigi twin_ang tmp_ebsd gi_twins mori_gi mori_ang mori_ax...
                  find_twins twin_combos twin_count len_twin
        end
    end
    
    % Add relevant items to myEBSD structure
    myEBSD.Recon.ID = Recon_pts;
    Parent.Indices = Par_Ind;      % Single Parent Indices
%     Parent.Clean.ID = Parcount;          % Overall Grain ID
    
    % No twins exist
    if num_twins == 0
        Twin=[];
    else
%         Parent.Twinned.Indices = ParTwin;    % Twinned-Parent Indices
%         Parent.Twinned.ID = ParTwinID;       % Twinned-Parent ID
        Twin.Merged = twinpar_Ind;           % Parent-Twin system indices
        Twin.Parent.Indices = twinparind;    % Twinned-Parent Indices for parent
%         Twin.Parent.ID = twinparID;          % Twinned-Parent ID
        Twin.Parent.Or = Partwin_Or;         % Twinned-Parent Orientation
        Twin.Single.Indices = twin_Ind;      % Twin Indices
%         Twin.Single.ID = twinID;             % Twin ID
        Twin.Single.Or = Twinonly_Or;        % Twin Orientation
    end
end

