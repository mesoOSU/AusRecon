function [myEBSD] = Recon_SingSpace(myEBSD,IP,OP,numOrs)
    %% Define Accrued Variables
    CS_T = myEBSD.CS{2};                % Transformed Phase
    CS_R = myEBSD.CS{3};                % Reconstructed Phase
    SS = myEBSD.SS;                     % Sample Symmetry
    Phases = myEBSD.Phase.Name;         % Phase names
    n_phases = myEBSD.Phase.n_phases;   % Number of total phases
    Trans_ebsd = myEBSD.TransEbsd;      % Transformation-phase ebsd
    T2R = myEBSD.TransMats{1};          % Transformation to Reconstruction phase
    R2T = myEBSD.TransMats{2};          % Transformation from Reconstruction Phase
    adjpts = myEBSD.Adj_pts;            % Adjacent Points Array
    Or_guesses = myEBSD.Or_guess;       % Orientation Guesses
    IP_wts = myEBSD.Weights.inplane;         % In-plane weights
    c2g_wts = myEBSD.Weights.outplane;        % Initial c2g weights
    rec_space = myEBSD.rec_space;       % Reconstruction space
    Trans_MODF = myEBSD.MODF;           % Misorientation Distribution Function
    Like_ind = myEBSD.Like_Inds;        % Indices for Low-likelihood regions
    psi = myEBSD.noise.psi;              % Psi value 
    Trans_ID = myEBSD.Phase.ID{1};      % Phase ID for Transformation phase
    Recon_ID = myEBSD.Phase.ID{2};      % Phase ID for Reconstruction phase    
    %% Define Everything Else/Graph Setup
    % Array of transformed phase orientations
    Trans_Ors=Trans_ebsd.orientations;

    % Now create our reconstruction ebsd map that, prior to any graph
    % cutting, is equivalent to the post-transformation ebsd microstructure
    recon_ebsd = Trans_ebsd;
    recon_ebsd.CS = CS_R;
    % Set the variable for the maximum austenite microstructure
    max_recon_ebsd = recon_ebsd;
    
    % Permutate Or guess list
    rearrange = randperm(length(Or_guesses));
    Or_guesses = Or_guesses(rearrange);

    % Lengths of pertinent arrays
    len_guess = length(Or_guesses);
    len_like = length(Like_ind);

    IP_reg=IP(1);         % Regularization of IP weights
    IP_scale=IP(2);   % Scaling probabilities for IP weights
    g2t_reg = OP(1);  % 2.0e0 works well
    g2t_scale = OP(2);   % Scaling probabilities for guess to sink
    
    % Our default will run through all guess orientations plus some
    % fraction of the number of likelihood guesses
    if isempty(numOrs)
        numOrs = round(len_guess+(0.20*len_like));      % # of total orientations
    end

    % Preallocation for likelihood weighting
    like_sum = zeros(len_like,1);

    % Establish the in-plane weights
    IP_wts = (IP_wts+IP_reg)*IP_scale;

    % Add source node
    di_graph=digraph;

    % Call function to set up the initial graph
    [di_graph,endnode,sinknode] = graph_setup(di_graph,adjpts,IP_wts,c2g_wts,2);

    % Counter for the orientation changes w/r/t austenite variants that
    % maximize the likelihood of the microstructure guess
    OR_counter = 1;

    % Keep track of the maximum c2g_wts for each variant in an austenite guess
    % so that we can choose the best pick
    max_c2g_wts = 1;
    
    % Initial values so we don't have to run through everything again
    Ph_trans = ones(length(Trans_ebsd),1).*Trans_ID;
    
    % Counters and the like
    like_count = 1;
    min_transpnts = round(0.0001*length(Trans_ebsd));
    len_sort = 0;
    untrans_count = 1;
    
    c2g_wts_final = c2g_wts;

    %% Perform Graph Cut
    % Display Total Number of Iterations that code will run through
    display(sprintf('----Total # of Iterations to be Performed: %d----',numOrs))

    % Iterative loop to fill out the rest of the austenite microstructure
    for iter = 1:numOrs
        % Working in Pre-Transformation Phase Space
        if strcmp(rec_space,Phases{2})
            if iter>len_guess
                % Take into account points that have yet to be transformed
                if untrans_count <= len_sort && untrans_ID == 0
                    Guess_Or = postOr_guess(untrans_count);
                    untrans_count = untrans_count+1;
                else 
                % Compute guess from low-likelihood regions 
                loc_guess = like_ID(like_count);
                mart_lst = Trans_Ors(Like_ind{loc_guess,:});
                loc_odf=calcODF(symmetrise(mart_lst)*T2R,'kernel',psi);
                [mm,Guess_Or] = max(loc_odf);
                like_count = like_count+1;
                end
            else
                % Use pre-computed guess orientations
                Guess_Or = Or_guesses(iter);
            end
            % Generate guess OR ODF and add weights from guess to sink
            From_Or_guess=symmetrise(Guess_Or)*R2T;
            From_guess_ODF=calcODF(From_Or_guess,'kernel',psi);
            % Now set the g2t_wts
            g2t_wts=(eval(From_guess_ODF,Trans_Ors)+g2t_reg).*g2t_scale;
            
        else % Working purely in Post-Transformation phase space
            
            if iter>len_guess
                % Take into account points that have yet to be transformed
                if untrans_count <= len_sort && untrans_ID == 0
                    Guess_Or = postOr_guess(untrans_count);
                    untrans_count = untrans_count+1;
                else 
                % Compute guess from low-likelihood regions 
                ind = Like_ind{like_count,:};
                likeOrs = Trans_Ors(ind);
                like_ci = Trans_ebsd.ci(ind)*100;
                likeID = [1:length(like_ci)]';
                like_pt = datasample(likeID,1,1,'Weights',like_ci);
                Guess_Or = likeOrs(like_pt);
                like_count = like_count+1;
                end
            else
                % Use pre-computed guess orientations
                Guess_Or = Or_guesses(iter);
            end
            % Generate and add weights from guess to sink (g2t_weights)
            Misos = inv(Guess_Or).*Trans_Ors;
            % Now set the g2t_wts
            g2t_wts=(eval(Trans_MODF,Misos)+g2t_reg).*g2t_scale;
        end
        
        % Eliminate any possible negative weights
        g2t_wts(find(g2t_wts<0))=0;

        % Increase regularization parameters since untransformed regions are 
        % most likely to have low confidence indexes
        if isempty(numOrs)==0 && iter > len_guess
            tmpcount = 1;
            index_update = [];
            checkindex = isempty(index_update);
            while checkindex==1
                % Set up temporary scaling
                g2t_wts=g2t_wts*tmpcount;

                % Add edges corresponding to the guess
                di_graph=rmedge(di_graph,(1:length(Trans_Ors))+endnode,sinknode);
                di_graph=addedge(di_graph,(1:length(Trans_Ors))+endnode,sinknode,g2t_wts);

                % Perform graph cut 
                [mf_c,gf,cs,ct]=maxflow(di_graph,1,sinknode);

                % Copy the cut for convenience
                ctcopy = ct;
                ctcopy(end)=[];

                % Update the orientations that need updating
                index_update = ctcopy-endnode;
                checkindex=isempty(index_update);
                tmpcount=tmpcount+1;
                if tmpcount > 3
                    checkindex=0;
                end
                g2t_wts=g2t_wts/tmpcount;
            end
        else
            % Add edges corresponding to the guess
            di_graph=rmedge(di_graph,(1:length(Trans_Ors))+endnode,sinknode);
            di_graph=addedge(di_graph,(1:length(Trans_Ors))+endnode,sinknode,g2t_wts);

            % Perform graph cut 
            [mf_c,gf,cs,ct]=maxflow(di_graph,1,sinknode);

            % Copy the cut for convenience
            ctcopy = ct;
            ctcopy(end)=[];

            % Update the orientations that need updating
            index_update = ctcopy-endnode;
        end

        % Update current weights with the new weights based on each variant
        tmp_wts = c2g_wts;
        tmp_wts(index_update) = g2t_wts(index_update);
        sum_c2g_wts = sum(c2g_wts);
        sum_tmp_wts = sum(tmp_wts);

        % Test current and guess microstructure likelihoods
        r_prob = rand;
        like_frac = sum_tmp_wts./sum_c2g_wts;
        like_prob = like_frac > r_prob;

        % If microstructure is accepted, update current microstructure with the
        % guess microstructure
        if length(index_update)>min_transpnts && like_prob
            c2g_wts = tmp_wts;
            
            if strcmp(rec_space,Phases{1})
                % Now choose austenite orientation randomly based on the max ODF 
                % austenite orientation in that region
                tmp_odf=calcODF(symmetrise(Trans_Ors(index_update))*T2R,'kernel',psi);
                [mm,Guess_Or] = max(tmp_odf);
                % Generate guess OR ODF and add weights from guess to sink
                From_Or_guess=symmetrise(Guess_Or)*R2T;
                From_guess_ODF=calcODF(From_Or_guess,'kernel',psi);
                % Now set the g2t_wts
                tmpg2t_wts=(eval(From_guess_ODF,Trans_Ors)+g2t_reg).*g2t_scale;
                c2g_wts_final(index_update) = tmpg2t_wts(index_update);
            end
            
            % Now set the new austenite guess
            recon_ebsd.orientations(index_update)=Guess_Or;

            % Keep track of martensitic points that have yet to be transformed
            Ph_trans(index_update) = Recon_ID;
            
            % Remove and add new edge corresponding to the current to grain edges
            di_graph = rmedge(di_graph,(1:length(Trans_Ors))+1,(1:length(Trans_Ors))+endnode);
            di_graph=addedge(di_graph,(1:length(Trans_Ors))+1,(1:length(Trans_Ors))+endnode,c2g_wts);

            % Keep track of the austenite guess that maximizes the overall
            % likelihood of austenite microstructure
            if sum_tmp_wts > max_c2g_wts
                max_c2g_wts = c2g_wts;      
                max_recon_ebsd = recon_ebsd;
            end

            % Counter to track the chosen orientations that change
            OR_counter = OR_counter+1;
        end  

        % Display iterations to show progress of code
        display(sprintf('---------------Iteration #: %d---------------', iter))

        % If we're through all of the generated austenite guesses, sum the
        % regions of localized likelihoods
        if iter == len_guess
            % Find untransformed points (Ut), their length, and
            % Reconstructed phase points (Re)
            Ut_pts = find(Ph_trans==Trans_ID);
            num_untranspts = length(Ut_pts);
            Re_pts = find(Ph_trans==Recon_ID);

            if (num_untranspts < 100)
                    len_sort = 1;
                    untrans_ID = 1;
            else

                % Create duplicate ebsd dataset corresponding to only the
                % martensite points that haven't been transformed  yet
                untrans_ebsd = Trans_ebsd;
                if n_phases > 1
                    untrans_ebsd.phase = Ph_trans;
                    untrans_ebsd(Phases{1}).orientations = orientation('euler',[1,1,1],CS_T);
                    untrans_ebsd(Phases{2}).orientations = orientation('euler',[0,0,0],CS_R);
                else 
                    untrans_ebsd(Ut_pts).orientations = orientation('euler',[1,1,1],CS_T);
                    untrans_ebsd(Re_pts).orientations = orientation('euler',[0,0,0],CS_T);
                end
                % Calculate Grains
                [grains,untrans_ebsd.grainId,untrans_ebsd.mis2mean] = calcGrains(untrans_ebsd);
                % Assign only untransformed grains 
                if n_phases > 1
                    untrans_grains = grains(Phases{1});
                else
                    grn_Ut = find(grains.meanOrientation == orientation('euler',[1,1,1]));
                    untrans_grains = grains(grn_Ut);
                end

                % Sort grains, both indexed and nonindexed, by grain size and
                % then arrange in descending order
                untrans_grains(untrans_grains.grainSize<50)=[];
                untrans_grain_area = untrans_grains.area;
                [sort_area,sorted_id] = sort(untrans_grain_area,'descend');
                len_sort = length(sorted_id);
                untrans_ID = isempty(sorted_id);

            for k = 1:len_sort
                tmpId = find(untrans_ebsd.grainId==untrans_grains.id(sorted_id(k)));
                tmpOrs = Trans_Ors(tmpId);
                
                if strcmp(rec_space,Phases{2})
                    tmp_odf=calcODF(symmetrise(tmpOrs)*T2R,'kernel',psi);
                    [mm,tmp_Or] = max(tmp_odf);
                    clear tmp_odf mm
                else
                    loc_ci = Trans_ebsd.ci(tmpId)*100;
                    transID = [1:length(loc_ci)]';
                    trans_pt = datasample(transID,1,1,'Weights',loc_ci);
                    tmp_Or = tmpOrs(trans_pt);
                    clear loc_ci transID trans_pt
                end
                postOr_guess(k) = tmp_Or;
                clearvars tmpId tmpORs tmp_Or
            end 
            clearvars mart_grains grains untrans_ebsd untrans_grain_area sort_area sorted_id 
            end
        end

        % If we're through the untransformed regions, move on to regions of 
        % low-likelihood
        if (untrans_count > len_sort) || (untrans_ID == 1)
            % Save microstructure post martensite region utilization
    %         save ausrecon_postgrains.mat
            % For each of the smaller regions, compute the likelihood sum 
            for iter_like = 1:len_like
                like_sum(iter_like) = sum(c2g_wts(Like_ind{iter_like,:}));
            end
            % Sort these regions of summed likelihoods in order from least to
            % greatest to use as our next set of guesses
            [~,like_ID] = sort(like_sum);
        end
    end
    
    % Index points that transformed
    Recon_pts = find(Ph_trans==Recon_ID);
    
    % If we have untransformed points and our initial ebsd dataset had more
    % than one phase, assign points
    if length(unique(Ph_trans)) > 1 && n_phases > 1
        recon_ebsd.phase = Ph_trans;
        find_untrans = find(Ph_trans==Trans_ID);
        recon_ebsd(find_untrans).CS = CS_T;
    end
    
    if strcmp(rec_space,Phases{2})
        c2g_wts_final = c2g_wts;
    end
    
    % Add variables to myEBSD structure
    myEBSD.Recon.Ebsd = recon_ebsd;
    myEBSD.Recon.TransformedPoints = Recon_pts;
    myEBSD.Recon.MaxEbsd = max_recon_ebsd;
    myEBSD.Recon.MaxLikelihood = max_c2g_wts;
    myEBSD.Recon.Likelihood = c2g_wts_final;


end
