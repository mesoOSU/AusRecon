function [myEBSD] = GridStructure(myEBSD,xdim,ydim,Grid)
    % Initialization Parameters
    rec_space = myEBSD.rec_space;
    Phases = myEBSD.Phase.Name;
    ebsd = myEBSD.Ebsd;
    trans_ebsd = myEBSD.TransEbsd;
    T2R = myEBSD.TransMats{1};
    psi = myEBSD.noise.psi;
    conf_wts = trans_ebsd.ci;
    
    %% Grid Setup
    % Minimum and maximum x,y coordinates
    xmin = min(ebsd.x);
    xmax = max(ebsd.x);
    ymin = min(ebsd.y);
    ymax = max(ebsd.y);
    
    if strcmp(rec_space,Phases{2})
        % Parameter that user can input for the box size and overlap
        x_ov = round(xdim/4);
        y_ov = round(ydim/4);
    else
        x_ov = 0;
        y_ov = 0;
    end
    
    xstart = xmin;
    ystart = ymax-ydim;
    
    % Show grid to user to see if they need to re-do
    if Grid == 1
       plotEBSD(myEBSD,ebsd,0);
    end
    
    % Set  my flag/counter value
    fin = 0;
    counter = 0;
    % Set up our austenite orientation guesses
    while fin == 0
        counter = counter+1;
        if ((xstart+xdim)>xmax) && (ystart<ymin)
            rr = [xstart ymin xmax-xstart ydim];
            fin = 1;
        elseif (xstart+xdim)>xmax
            rr = [xstart ystart xmax-xstart ydim];
            xstart = xmin-(xdim-x_ov);
            ystart = ystart-(ydim-y_ov);
        elseif ystart<ymin
            rr = [xstart ymin xdim ydim];
        else
            rr = [xstart ystart xdim ydim];
        end
        
        % Plot grid if desired
        if Grid==1
            rectangle('position',rr);
        end

        % Based on reconstruction space, either produce list of potential
        % pre or post-transformation orientations to use as either guesses
        % (austenite space) or misorientation comparisons (martensite
        % space), respectively
        cond=inpolygon(trans_ebsd,rr);
        tmp_ebsd = trans_ebsd(cond);
%         if isempty(tmp_ebsd)
%             keyboard
%         end
        trans_tmp = tmp_ebsd.orientations;
        % If too few points to choose from, skip this region altogether.
        % Else, determine which space we're working in and specify results
        % based on that
        if length(tmp_ebsd) > 200
            if strcmp(rec_space,Phases{2})
                loc_odf=calcODF(symmetrise(trans_tmp)*T2R,'kernel',psi);
                [mm,aus_locmax] = max(loc_odf);
                Or_guesses(counter,1) = aus_locmax;
                display(sprintf('---------------Martensite Region # %d Completed---------------', counter))
                clear loc_odf mm aus_locmax
            else
                loc_pts = find(cond);
                loc_ci = conf_wts(loc_pts).*100;
                transID = [1:length(loc_pts)]';
                trans_guess = datasample(transID,1,1,'Weights',loc_ci);
                Or_guesses(counter,1) = trans_tmp(trans_guess);
                clear loc_pts loc_ci transID trans_guess
            end
        else
            counter = counter-1;
        end
        
        % Shift our starting x-position
        xstart = xstart+(xdim-x_ov);
        clearvars tmp_ebsd cond rr trans_tmp
    end
    
    clear xstart ystart counter fin

    % Eliminate any repeat austenite orientations
    Or_guesses = unique(Or_guesses);
   
    %% Construct Smaller Grid for Reconstruction
    % Reset a lot of variables from above
    box_like = round(xdim/2);
    xstart = xmin;
    ystart = ymax-box_like;
    counter = 0;
    fin = 0;
    like_ind = [];

    % Now set up a smaller boxed grid that we'll use to pick out the spots
    % where poor likelihoods exist post
    while fin == 0

        counter = counter+1;

        if ((xstart+box_like)>xmax) && (ystart<ymin)
            rr = [xstart ymin (xmax-xstart) (yprev-ymin)];
            fin = 1;
        elseif (xstart+box_like)>xmax
            rr = [xstart ystart xmax-xstart box_like];
            yprev = ystart;
            xstart = xmin-box_like;
            ystart = ystart-box_like;
        elseif ystart<ymin
            rr = [xstart ymin box_like (yprev-ymin)];
        else
            rr = [xstart ystart box_like box_like];
        end
        cond=inpolygon(trans_ebsd,rr);
        like_ind{counter,1} = find(cond);
        clearvars cond rr

        xstart = xstart+box_like;
    end
    
    myEBSD.Or_guess = Or_guesses;
    myEBSD.Like_Inds = like_ind;
end
