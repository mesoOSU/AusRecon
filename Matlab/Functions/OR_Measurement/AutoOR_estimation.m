function [myEBSD] = AutoOR_estimation(myEBSD,vis_m,vis_a,num_Ors,plt_ksi)
    tic;
    % Array to save data in 
    CS_T = myEBSD.CS{2};
    CS_R = myEBSD.CS{3};
    SS = myEBSD.SS;
    
    itermax = 10;
    iters = 0;
    
    while iters < itermax
        iters = iters + 1;
        % Extract region to perform OR estimation on
        if iters == 1
            FullEbsd = myEBSD.Ebsd;
        end
        [myEBSD,Ebsd,FullEbsd,austenite_initial] = AutoOR_Grn(myEBSD,FullEbsd,iters);

        TransID = find(Ebsd.phase == myEBSD.Phase.ID{1});
        ebsd = Ebsd(TransID);
        martensite=ebsd.orientations;
        %martensite(find(ebsd.ci<0))=[];
        martensite.CS=CS_T;
        size_martensite=length(martensite);

        %downsample ebsd dataset to number of orientations
        keep=randperm(length(martensite));
        if length(martensite) > num_Ors
            martensite=martensite(keep(1:num_Ors));
        else
        end
        %plot the martensite pole figure
        if (vis_m == 1)
            figure; plotPDF(martensite,Miller({0,0,1},martensite.CS),'antipodal','points','all','marker','.');
        end

        %% Generate initial guess for parameters to be estimated

        %initial guess for ksi based on paper
        ksi_initial=[5.26 10.3 10.53];      % KS
    %     ksi_initial=[0.00,9.74,9.74];       % NW
        % ksi_initial = ksi_1st_guess;

        % ksi_initial=[3.3,8.5,8.9];        % KS-Like
        %ksi_initial=[4.37,7.79,9.24];
        %ksi_initial=[2.3114    8.9269    9.1332];  % 1st Samp Exp Obs
        % ksi_initial = [2.986,8.228,8.584];      % Crop3 Ksi average for good fits
        halfwidth_in=2.5*degree;
        %initial guess on halfwidth
        %halfwidth_initial=2*degree;

        % initial austenite guess from global pole figure

        [T2R,flag]=calc_T2R(ksi_initial,CS_R,CS_T);

        if flag
            warning('non-physical initial ksi angles!');
        end

        %% set parameters for inference

        % Set mu and sigma for prior probability on ksi. Values based on rough
        % Estimates from Yardley Payton 2014 conference paper
        % ksi_prior_mu=[3,8,9];
        ksi_prior_mu = [5,9,10];
        % ksi_prior_sigma=[1.2,1.2,1.2];
        ksi_prior_sigma = [2,2,2];

        % Noise is modeled as a unimodal odf centered on the cube orientation.
        % Parameter to be estimated by Bayesian inference is the halfwidth of the
        % odf kernel. Halfwidth distribution is assumed folded Gaussian. This
        % approximates uniform from 0->~1 degree then decaying at larger noise values
        halfwidth_prior_mu=1;
        halfwidth_prior_sigma=2;

        %austenite - twp prior considered for austenite based on cases. If ksi
        %angles can be carefully measured then austenite prior will be a unimodal
        %odf about the global pole figure modal ODF. If initial ksi values are
        %approximate then austenite prior is uniform over the fundamental zone

        %leftover can ignore
        halfwidth_austenite_prior=3*degree;
        austenite_prior_odf=uniformODF(CS_R,SS);


        prior_pars=struct;
        prior_pars.ksi_prior_mu=ksi_prior_mu;
        prior_pars.ksi_prior_sigma=ksi_prior_sigma;
        prior_pars.halfwidth_prior_mu=halfwidth_prior_mu;
        prior_pars.halfwidth_prior_sigma=halfwidth_prior_sigma;
        prior_pars.austenite_prior_odf=austenite_prior_odf;
        prior_pars.CS_A=CS_R;
        prior_pars.CS_M=CS_T;
        prior_pars.SS=SS;

        %% MAP estimate of parameters by optimization

        options=optimset('fminsearch');
        options=optimset(options,'display','iter','algorithm','sqp');
        optimfunc= @(samples) -posterior_pdf_fminunc(samples,prior_pars,martensite);
        x0=[ksi_initial,halfwidth_in/degree];
        
        % Ensure the constraint that ksi_1 < ksi_2 and ksi_3. Since we
        % don't know for certain if ksi_2 SHOULD be < ksi_3, force this to
        % occur if the first constraint is met.
        init_guess=0;
        count = 0;
        while init_guess == 0
            count = count+1;
        % Optimization function outside of ML add-on
            [MAPpars,loglike,exitflag,output]=fminsearch(optimfunc,x0,options);
        
            if (MAPpars(1) < MAPpars(2) && MAPpars(3)) || count > 2
                init_guess = 1;
                % Enforce constraint if need be
                if(MAPpars(2) > MAPpars(3))
                    MAPpars(2) = (MAPpars(3)-1e-3);
                end
            end
        end
        ksi_1st_guess=MAPpars(1:3);
        halfwidth_in=MAPpars(4)*degree;

        %% Set MCMC sampler parameters

        %burnin is number of samples to disregard at the beginning maybe try 100
        %for paper
        burnin=0;
        num_samples=1;

        %initialize counters and placeholders
        MAP_step=0;
        MAP_likelihood=-1e6;


        %ksi values update parameters
        scale=0.15;
        %ksi_width=[0.05, 0, 0;
        %            0, 0.17, 0.15;
        %            0, 0.15, 0.16]*scale;

        % ksi_width=[0.05, 0, 0;
        %             0, 0.17, 0.0;
        %            0, 0.0, 0.17]*scale;
        %
        ksi1_width=.2*scale;
        ksi2_width=.08*scale;
        ksi3_width=.08*scale;
        halfwidth_width=0.02;
        austenite_width=1*degree;

        %leftover needs to be cleaned up
        austenite_current=austenite_initial;

        %% componentwise MCMC sampler

        ksi1_current=ksi_1st_guess(1);
        ksi2_current=ksi_1st_guess(2);
        ksi3_current=ksi_1st_guess(3);
        halfwidth_current=halfwidth_in;

        %austenite_posterior=austenite_current;
        ksi1_posterior=ksi1_current;
        ksi2_posterior=ksi2_current;
        ksi3_posterior=ksi3_current;
        halfwidth_posterior=halfwidth_current;


        M=martensite;
        num_mart_Ors=100;%;num_martensite_orientations;
        for ii=1:num_samples+burnin

            aaa=[0,0,0];
            mm=randperm(length(M));
            martensite=M(mm<=num_mart_Ors);
            halfwidth_proposal=halfwidth_current;

            %    halfwidth_proposal=(randn*halfwidth_width+halfwidth_current/degree)*degree;
            % ksi 3 sample in order ksi3, ksi2, ksi1, halfwidth

            p_current=martensite_posterior_log_likelihood(martensite,[ksi1_current,ksi2_current,ksi3_current],halfwidth_current,austenite_current,prior_pars);

            flag=1;
            flag_counter=0;
            skip=0;
            while flag
                flag_counter=flag_counter+1;
                ksi3_proposal=abs(randn*ksi3_width+ksi3_current);
                [Variants flag]=YardleyVariants([ksi1_current,ksi2_current,ksi3_proposal]);
                if flag_counter>=10
                    flag=0
                    ksi3_proposal=ksi3_currrent;
                    skip=1;
                end
            end

            [T2R,flag]=calc_T2R([ksi1_current,ksi2_current,ksi3_proposal],CS_R,CS_T);

            austenite=symmetrise(martensite)*T2R;

            [austenite_proposal] = global_pole_figure_estimation(austenite,CS_R,SS,1);

            for kk=1:length(austenite_proposal)
                temp(kk)=martensite_posterior_log_likelihood(martensite,[ksi1_current,ksi2_current,ksi3_proposal],halfwidth_current,austenite_proposal(kk),prior_pars);
            end
            id=find(temp==max(temp),1);
            p_proposal=temp(id);%martensite_posterior_log_likelihood(martensite,[ksi1_proposal,ksi2_proposal,ksi3_proposal],halfwidth_proposal,austenite_proposal,prior_pars);
            clear temp

            %stricly the log probability
            p_accept=p_proposal-p_current;

            if skip
                accept = false;
            else
            accept=log(rand)<p_accept;
            accept_log(ii)=accept;
            end

            if accept
                austenite_current=austenite_proposal(id);
                ksi3_current=ksi3_proposal;
                aaa(3)=1;
            end

            ksi3_posterior=[ksi3_posterior;ksi3_current];

            % ksi 2 sample in order ksi3, ksi2, ksi1, halfwidth

            p_current=martensite_posterior_log_likelihood(martensite,[ksi1_current,ksi2_current,ksi3_current],halfwidth_current,austenite_current,prior_pars);

            flag=1;
            flag_counter=0;
            skip=0;
            while flag
                flag_counter=flag_counter+1;
                ksi2_proposal=abs(randn*ksi2_width+ksi2_current);
                [Variants flag]=YardleyVariants([ksi1_current,ksi2_proposal,ksi3_current]);
                if flag_counter>=10
                    flag=0
                    ksi2_proposal=ksi2_currrent;
                    skip=1;
                end
            end

            [T2R,flag]=calc_T2R([ksi1_current,ksi2_proposal,ksi3_current],CS_R,CS_T);

            austenite=symmetrise(martensite)*T2R;

            [austenite_proposal] = global_pole_figure_estimation(austenite,CS_R,SS,1);

            for kk=1:length(austenite_proposal)
                temp(kk)=martensite_posterior_log_likelihood(martensite,[ksi1_current,ksi2_proposal,ksi3_current],halfwidth_current,austenite_proposal(kk),prior_pars);
            end
            id=find(temp==max(temp),1);
            p_proposal=temp(id);%martensite_posterior_log_likelihood(martensite,[ksi1_proposal,ksi2_proposal,ksi3_proposal],halfwidth_proposal,austenite_proposal,prior_pars);
            clear temp

            p_accept=p_proposal-p_current;

            if skip
                accept = false;
            else
            accept=log(rand)<p_accept;
            accept_log(ii)=accept;
            end

            if accept
                austenite_current=austenite_proposal(id);
                ksi2_current=ksi2_proposal;
                aaa(2)=1;
            end

            ksi2_posterior=[ksi2_posterior;ksi2_current];

            % ksi 1 sample in order ksi3, ksi2, ksi1, halfwidth

            p_current=martensite_posterior_log_likelihood(martensite,[ksi1_current,ksi2_current,ksi3_current],halfwidth_current,austenite_current,prior_pars);

            flag=1;
            flag_counter=0;
            skip=0;
            while flag
                flag_counter=flag_counter+1;
                ksi1_proposal=abs(randn*ksi1_width+ksi1_current);
                [Variants flag]=YardleyVariants([ksi1_proposal,ksi2_current,ksi3_current]);
                if flag_counter>=10
                    flag=0
                    ksi1_proposal=ksi1_currrent;
                    skip=1;
                end
            end

            [T2R,flag]=calc_T2R([ksi1_proposal,ksi2_current,ksi3_current],CS_R,CS_T);

            austenite=symmetrise(martensite)*T2R;

            [austenite_proposal] = global_pole_figure_estimation(austenite,CS_R,SS,1);

            for kk=1:length(austenite_proposal)
                temp(kk)=martensite_posterior_log_likelihood(martensite,[ksi1_proposal,ksi2_current,ksi3_current],halfwidth_current,austenite_proposal(kk),prior_pars);
            end
            id=find(temp==max(temp),1);
            p_proposal=temp(id);%martensite_posterior_log_likelihood(martensite,[ksi1_proposal,ksi2_proposal,ksi3_proposal],halfwidth_proposal,austenite_proposal,prior_pars);
            clear temp

            p_accept=p_proposal-p_current;

            if skip
                accept = false;
            else
            accept=log(rand)<p_accept;
            accept_log(ii)=accept;
            end

            if accept
                austenite_current=austenite_proposal(id);
                ksi1_current=ksi1_proposal;
                aaa(1)=1;
            end
            ksi1_posterior=[ksi1_posterior;ksi1_current];
        end
        
        % If the likelihood value is large enough, conclude with this OR.
        % If not, and we've gone through 10 iterations with no successful
        % OR, return error and make user input manually.
        if p_proposal > 4e2
            iters = 1e2;
        elseif iters == itermax
            error('Automatic Determination of OR Failed. Try Manual Segmentation of Potential PAG!')
        end
    end

    % Display current ksi values
    
    ksi_out=[ksi1_posterior(burnin+1:end),ksi2_posterior(burnin+1:end),ksi3_posterior(burnin+1:end)];
    %out2=[out,halfwidth_posterior(burnin+1:end)/degree];
    sampling_rates=[numel(find(diff(ksi1_posterior))),numel(find(diff(ksi2_posterior))),numel(find(diff(ksi3_posterior)))]/num_samples;
    ksi_curr = mean(ksi_out);

    % Plot the Discrete ODFs
    [R2T,flag]=calc_R2T(ksi_curr',CS_R,CS_T);
    curr_mart = symmetrise(austenite_current)*R2T;
    curr_mart_odf = calcODF(curr_mart,CS_T,SS,CS_R);
    % figure; plotPDF(curr_mart_odf,Miller({0,0,1},{1,1,0},{1,1,1},CS_M),'antipodal','points','all','marker','.');
    if vis_a == 1
        martensite=generate_simulated_data(austenite_proposal,[ksi_curr(1),ksi_curr(2),ksi_curr(3)],halfwidth_current,num_Ors,CS_T);
        figure; plotPDF(martensite,Miller({0,0,1},martensite.CS),'antipodal','points','all','marker','.');
    end
    
    myEBSD.OR = ksi_curr;
    myEBSD.ORLikelihood = p_proposal;
    myEBSD.ksi.ksi1 = ksi_out(:,1);
    myEBSD.ksi.ksi2 = ksi_out(:,2);
    myEBSD.ksi.ksi3 = ksi_out(:,3);
    myEBSD.noise.halfwidth = halfwidth_proposal;
    elapsed_time = toc/60
    if plt_ksi == 1
        for i = 1:3
            figure; histogram(ksi_out(:,i),20)
            if i == 1
                title('\xi_1')
            elseif i == 2
                title('\xi_2')
            else
                title('\xi_3')
            end
            xlim([min(ksi_out(:,i))-1,max(ksi_out(:,i)+1)])
        end
    end
end