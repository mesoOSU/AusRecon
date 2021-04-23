function plotEBSD(varargin)

    Ebsd = varargin{1};
    myEBSD = varargin{2};
    
    if length(varargin)==3
        mbar = 'MicronBar';
        mbar_flg = 'off';
    else
        mbar = [];
        mbar_flg = [];
    end

    % Number of phases and the phase name
    n_phases = length(unique(Ebsd.phase));
    Phases = myEBSD.Phase.Name;
    
    if n_phases == 1
        key = ipfHSVKey(Ebsd);
        key.inversePoleFigureDirection=zvector;
    else
        key = ipfHSVKey(Ebsd(Phases{1}));
        key.inversePoleFigureDirection=zvector;
    end
%     % Plot the key if desired
%     if Key == 1
%         plot(key,'MicronBar','off')
%     end
    
    % Determine if we need to partition Ebsd dataset
    partitions = ceil(length(Ebsd)/4e6);
    
    % Max cartesian coordinates
    max_X = max(Ebsd.x);
    max_Y = max(Ebsd.y);
    
    part_count = 1;
    while part_count <= partitions
    
        if max_X > max_Y
            if part_count == 1
                x_beg = min(Ebsd.x);
            else
                x_beg = x_end;
            end
            x_end = (part_count*max_X)/partitions;
            y_beg = min(Ebsd.y);
            y_end = max_Y;
        else
            if part_count == 1
                y_beg = min(Ebsd.y);
            else
                y_beg = y_end;
            end
            x_beg = min(Ebsd.x);
            x_end = max_X;
            y_end = (part_count*max_Y)/partitions;
        end
        % Coordinates for Ebsd portion we're plotting
        rr = [x_beg y_beg (x_end-x_beg) (y_end-y_beg)];
        cond=inpolygon(Ebsd,rr);
        tmpebsd = Ebsd(cond);
        
        % Now plot microstructures
            if n_phases == 1
                color = key.orientation2color(tmpebsd.orientations);
                figure; plot(tmpebsd,color,mbar,mbar_flg);
            else
                % Plot phase with most points first so as to avoid weird
                % IPF mappings
                for jj = 1:n_phases
                    len_ph(jj,1) = length(tmpebsd(Phases{jj}));
                    len_ph(jj,2) = jj;
                end
                % Sort them in descending order
                [~,sorted] = sort(len_ph(:,1),'descend');

                figure
                for i = 1:n_phases
                    hold on
                    ph_id = sorted(i);
                    color = key.orientation2color(tmpebsd(Phases{ph_id}).orientations);
                    plot(tmpebsd(Phases{ph_id}),color,mbar,mbar_flg);
                end
            end
        part_count = part_count+1;
    end      
end
