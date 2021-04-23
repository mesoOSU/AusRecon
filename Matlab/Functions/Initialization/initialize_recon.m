function [myEBSD] = initialize_recon(myEBSD)

    %% Define Crystal/Sample Symmetries and OR
    CS_T = myEBSD.CS{2};
    CS_R = myEBSD.CS{3};
    SS = myEBSD.SS;
    Phases = myEBSD.Phase.Name;
    OR = myEBSD.OR;
    psi = myEBSD.noise.psi;
    Ebsd = myEBSD.Ebsd;
    rec_space = myEBSD.rec_space;
    TransID = find(Ebsd.phase==myEBSD.Phase.ID{1});
    
    %% Calculate the Misorientation Data
    
    % Construct Transformation Matrices
    [T2R,flag]=calc_T2R(OR,CS_R,CS_T);
    [R2T,flag]=calc_R2T(OR,CS_R,CS_T);
    
    if myEBSD.Quad.Flag
        len = 1;
    else
        len = 4;
    end
    for i = 1:len
        Trans_ebsd = myEBSD.Quad.Ebsd{i};
        Trans_ebsd.CS = CS_T;

        nn=1;
        % Call function to compute adjacency array and corresponding
        % misorientations for the neighboring points
        [adjpts,mori,~] = adjpt_moris(Trans_ebsd,nn);

        % Evaluated misorientation distribution function values (serving as
        % in-plane weights)
        %Austin
%        trans_modf_vals=eval(myEBSD.MODF,mori);
        trans_modf_vals=alt_eval(myEBSD.MODF,mori);
        trans_modf_vals(trans_modf_vals<0)=0;

        % Set orientations for every point in the designated space equal to a 
        % specific, random austenite orientation. This will be our current image
        % that will continuously be updated based each iteration.
        if strcmp(rec_space,Phases{2})
            % Create uniform odf based on austenite crystal symmetry
            odfR = uniformODF(CS_R);
            recon_ebsd.orientations = calcOrientations(odfR,length(Trans_ebsd.orientations));
            % Calculate first set of weights corresponding to the random austenite
            % microstructure (TAKES LONG TIME!)
            Trans_from_Recon_guess=symmetrise(Trans_ebsd.orientations)*R2T;
            Trans_from_guess_ODF=calcODF(Trans_from_Recon_guess,'kernel',psi);
            % Evaluated 
            init_c2g_wts=eval(Trans_from_guess_ODF,Trans_ebsd.orientations);
        else
            init_c2g_wts = ones(length(Trans_ebsd),1);
        end
        
        % Add transformation ebsd datasets to relevant structure
        myEBSD.Adj_pts{i} = adjpts;
        myEBSD.Weights.inplane{i} = trans_modf_vals;
        myEBSD.Weights.outplane{i} = init_c2g_wts;
        myEBSD.TransEbsd{i} = Trans_ebsd;

    end
    % Add variables to myEBSD structure (adjacent points array, ebsd
    % corresponding to transformed and reconstructed points)
    myEBSD.TransMats{1} = T2R;
    myEBSD.TransMats{2} = R2T;
    
end