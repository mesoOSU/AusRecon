function [varargout] = BoundaryIntersect(varargin)

%%                      Function Description

% This function follows ASTM E112-13 standards (for grains) by combining 
% annealling twins with the parent s.t. the parent-twin system becomes a 
% single entity whose grain boundaries are defined by the parent grain.

% Single intercept lines are performed on any packet, block, and sub-block 
% boundaries for further analysis if the material is steel by randomly
% assigning a start and stop on the x_min and y_max edge, and then
% performing the same sort of standard segmentation analysis.

% FOR GRAINS ONLY!:

% The grain boundary lineal intercept method involves 12 total lineal 
% intercepts, with 8 vertical lines spanning the x-axis coordinates and 
% 4 horizontal lines spanning the y-axis, all lines being equally 
% distributed from each other.

% Input Variables: myESBD structure; grain boundary structure; if input, 
% the set of coordinates for a single user-defined line

% Output Variables: Modified myEBSD structure; modified grain structure
% containing all computed ASTM data

%%                      Function Implementation

    % Assign myEBSD and Grain structures (called Feature here)
    myEBSD = varargin{1};
    Feature = varargin{2};
    Bounds = Feature.boundary;
    
    % Length of packets, blocks, and variants (if they have already been
    % computed)
    if isfield(myEBSD,'Packets')
        if isfield(myEBSD.Packets,'Boundaries')
            LpB = length(myEBSD.Packets.Boundaries);
        else
            LpB = 0;
        end
        if isfield(myEBSD.Blocks,'Boundaries')
            LbB = length(myEBSD.Blocks.Boundaries);
        else
            LbB = 0;
        end
        if isfield(myEBSD.Variants,'SubBlockBoundaries')
            LsbB = length(myEBSD.Variants.SubBlockBoundaries);
        else
            LsbB = 0;
        end
    else
        LpB = 0;
        LbB = 0;
        LsbB = 0;
    end
    
    % Define parent microstructure ebsd data set and its length
    ParEb = myEBSD.Recon.Ebsd;
        
    % Reconstruction phase ID
    ReconID = myEBSD.Phase.ID{2};
    
    % Remove indexing for initially indexed austenite so only the
    % martensite points will actually be plotted
    OrigEbsd = myEBSD.Ebsd;
    InitAus = find(OrigEbsd.phase==ReconID);
    
    % If user wants to define their own line segment, use their
    % coordinates. Else, randomly confine a line to the minimum x value 
    % (left edge), a y value contained within this edge, and two randomly 
    % chosen x,y endpoints.
    if length(varargin)==3
        X1 = varargin{4}(1);
        X2 = varargin{4}(2);
        Y1 = varargin{4}(3);
        Y2 = varargin{4}(4);
    else
        X1 = min(ParEb.x);
        % Values greater than X, randomly choose one of these as X endpoint
        GrtX = find(ParEb.x > X1);
        X2 = ParEb(GrtX(randi(length(GrtX)))).x;
        % Find unique X values on the left edge and choose a Y value
        Xline = find(ParEb.x==X1);
        Y1 = ParEb(Xline(randi(length(Xline)))).y;
        % Now find a randomly chosen Y end point value greater than the
        % chosen Y start point
        Y2 = max(ParEb.y);
    end
    
    % Establish the line endpoints
    L1 = [X1,Y1];
    L2 = [X2,Y2];
    % Length of line
    lenLine = sqrt((L2(1)-L1(1))^2+(L2(2)-L1(2))^2);
    
    % Determine number of intersections and line-segment lengths
    [X,Y,SortedLineSegs] = Bounds.intersect(L1,L2);
    Xint = X(find(~isnan(X)));
    Yint = Y(find(~isnan(Y)));
    % Number of intersections
    P_i = length(SortedLineSegs);
    
    % Average grain size
    li = lenLine / P_i;
    % Surface Area
    Sv = 2/lenLine*1000;
    
    % Calculate ASTM Grain Size Number Per ASTM E112-13
    ASTM_Num = -6.643856*log10(lenLine/1000)-3.288;
    
    varargout{1} = [];
    % Add pertinent line segment information into relevant structures based
    % on whether the input boundaries are Packet, Block, Sub-block or Grain
    % boundaries.
    if length(Feature) == LpB
        myEBSD.Packets.LineIntercepts = [SortedXint;SortedYint];
        myEBSD.Packets.LineSegments = SortedLineSegs;
        myEBSD.Packets.LineLength = lenLine;
        myEBSD.Packets.EndPoints = [L1;L2]';
        
    elseif length(Feature) == LbB
        myEBSD.Blocks.LineIntercepts = [SortedXint;SortedYint];
        myEBSD.Blocks.LineSegments = SortedLineSegs;
        myEBSD.Blocks.LineLength = lenLine;
        myEBSD.Blocks.EndPoints = [L1;L2]';
        
    elseif length(Feature) == LsbB
        myEBSD.Variants.LineIntercepts = [SortedXint;SortedYint];
        myEBSD.Variants.LineSegments = SortedLineSegs;
        myEBSD.Variants.LineLength = lenLine;
        myEBSD.Variants.EndPoints = [L1;L2]';
        
    else % Do multiple lineal intercepts for ASTM grain size measurement
        
        % Extract triple point coordinates and set pixelated tolerance
        tpX = Bounds.triplePoints.x;
        tpY = Bounds.triplePoints.y;
        tpTol = 1;

        % Range of x and y values and lineal intercept line length
        Lx = range(ParEb.x);
        Ly = range(ParEb.y);
        TotLine = (Lx*4)+(Ly*8);

        % X-based values
        XPart = Lx / 10;
        XStrt = min(ParEb.x);
        XEnd = max(ParEb.x);   

        % Y-based values
        YPart = Ly / 6;
        YStrt = min(ParEb.y);
        YEnd = max(ParEb.y);

        % X-dependent line that runs across the x-axis coordinates at 
        % constant y values
        Xdep_X1 = ones(1,4)*XStrt;
        Xdep_X2 = ones(1,4)*XEnd;
        Xdep_Y1 = linspace(YStrt+YPart,YEnd-YPart,4);
        Xdep_Y2 = Xdep_Y1;

        % Y-dependent line that runs across the y-axis coordinates at 
        % constant x values
        Ydep_X1 = linspace(XStrt+XPart,XEnd-XPart,8);
        Ydep_X2 = Ydep_X1;
        Ydep_Y1 = ones(1,8)*YStrt;
        Ydep_Y2 = ones(1,8)*YEnd;

        % Number of lineal intercepts
        Sz = 12;

        for ii = 1:Sz

            if ii > 4
                X1 = Ydep_X1(ii-4);
                X2 = Ydep_X2(ii-4);
                Y1 = Ydep_Y1(ii-4);
                Y2 = Ydep_Y2(ii-4);
                lenLine = Ly;
            else
                X1 = Xdep_X1(ii);
                X2 = Xdep_X2(ii);
                Y1 = Xdep_Y1(ii);
                Y2 = Xdep_Y2(ii);
                lenLine = Lx;
            end

            % Establish the line endpoints
            L1 = [X1,Y1];
            L2 = [X2,Y2];
            
            % Length of line
            %         lenLine = sqrt((L2(1)-L1(1))^2+(L2(2)-L1(2))^2);

            % Determine number of intersections and line-segment lengths
            [X,Y,LineSegs] = Bounds.intersect(L1,L2);

            % Extract the intersections
            Ints = ~isnan(X);
            % Id grains of the intersections
            IntGrns = Bounds.grainId(Ints,:);

            % Find X and Y coordinates of intersections
            Xint = X(Ints);
            Yint = Y(Ints);

            % Extract coordinates of the unique grains at intersections
            [UnqIntGrns,UnqGrnIntId] = unique(IntGrns,'rows');
            UnqXint = Xint(UnqGrnIntId);
            UnqYint = Yint(UnqGrnIntId);
            %         UnqLineSeg = LineSegs(UnqGrnIntId);

            % Sort based on X intercepts unless x is constant, then sort by Y
            if X2 - X1 == 0
                [~,SortIntGrnId] = sort(UnqYint,'descend');
            else
                [~,SortIntGrnId] = sort(UnqXint);
            end

            % Extract the sorted, unique grain ids of intercepts and x and y 
            % coordinate intercepts. Although sorting is unneccessary, it 
            % allows us to follow the line from grain to grain either left to
            % right (if x isn't constant) or top to bottom
            SortedIntGrns = UnqIntGrns(SortIntGrnId,:);
            SortedXint = UnqXint(SortIntGrnId)';
            SortedYint = UnqYint(SortIntGrnId)';
            %         SortedLineSegs = UnqLineSeg(SortIntGrnId)';

            % Number of intersections
            P_i = length(SortedIntGrns)-1;
            
            % Check for triple points
            tpDsts = abs((L2(2)-L1(2)).*tpX-(L2(1)-L1(2)).*tpY+L2(1)*L1(2)-L2(2)*L1(1))/...
                sqrt((L2(2)-L1(2))^2+(L2(1)-L1(1))^2);
            
            % For square pixels, it's simple
            if length(ParEb.unitCell) == 4
                % Calculate the tolerance for triple point intersection
                % distances
                tolDst = 2*tpTol*abs(ParEb.unitCell(1,1));
            else
            % For hexagonal pixels, keep it simple and only use the
            % x-distance between pixel centroids
            XminHex = min(Ebsd.unitCell(:,1));
            XmaxHex = max(Ebsd.UnitCell(:,1));
            tolDst = tpTol*(XmaxHex-XminHex);
            end
            
            % Find the length of the number of triple points whose
            % tolerance is intersected by the lineal intercept line
            tpSegs = length(find(tpDsts < tolDst));
            
            % Adjust the number of intersections based on the number of
            % triple point intersections
            P_i = P_i + tpSegs*0.5;
              
            % Now check for any lineal intercepts that directly intersect a
            % triple point.
            % NOTE 1: The way MTEX searches for grain boundaries, triple
            % points don't necessarily (nor commonly) land on the vertices
            % of indexed pixels, but rather at points within the pixels.
            % NOTE 2: Because of this, if a triple point does happen to
            % land on the vertex of a pixel, it results in 3 detected 
            % grain boundary segmentations. Since ASTM standards call for 
            % 1.5 segmentations at triple points, we need to adjust P_i by 
            % subtracting 2 (since we already added the 0.5 in the previous
            % line)
            TPtChk = zeros(length(SortedXint),1);
            UnqIntGrn = unique(SortedIntGrns);
            TanChk = [];
            LenTanChk = zeros(length(SortedXint),1);
            
            % Loop through all of the intersected coordinates and find
            % those that happen to align with a triple point coordinate
            for kk = 1:length(SortedXint)
                TPtChk(kk) = any(tpX == SortedXint(kk) & tpY == SortedYint(kk));
            %             if kk <= length(UnqIntGrn)
            %                 TanChk{kk,1} = find(UnqIntGrn(kk) == SortedIntGrns);
            %             LenTanChk(kk) = length(find(UnqIntGrn(kk) == SortedIntGrns));
            end

%           % If any triple points exist, modify our intersection count
%           by reducing by 2xn, where n is the number of triple points
%           intersected
            if length(find(TPtChk)) > 0
                P_i = P_i - length(find(TPtChk == 1));
            end
            
            % If any tangents exist, modify our intersection count
            % NOTE: I haven't gotten this to work, so this does nothing.
            % Would need some directionality of the grain boundary line
            % segmentations that happen to run coincident with the lineal
            % intercept line.
            if any(LenTanChk > 2)
                TanPts = find(LenTanChk > 2);
                TripPts = find(TPtChk);
                TripPtGrns = SortedIntGrns(TripPts,:);
                TanPtGrns = SortedIntGrns(TanPts,:);

                for  jj = 1:length(TanPtGrns)
                    if any(TanPtGrns(jj,:) ~= TripPtGrns)
                        P_i = P_i -0.5;
                    end
                end
            end
            
            % Mu line segment length
            li(ii) = lenLine / P_i;
            % Number of intersections
            NSegs(ii) = P_i;
            % Surface Area
            Sv(ii) = 2/li(ii)*1000;
            
            % ASTM E112-13 grain size number
            ASTM_Nums(ii) = -6.643856*log10(li(ii)/1000)-3.288;

            % Save the intercepted grains for each line
            E112_13.GrnIntersects{ii} = SortedIntGrns;
            E112_13.XIntersects{ii} = SortedXint;
            E112_13.YIntersects{ii} = SortedYint;
        end
        
        % CI Functions
        CI95li = tinv([0.025,0.975],Sz-1);
        CI95ASTM = tinv([0.025,0.975],Sz-1);
        CI95Sv = tinv([0.025,0.975],Sz-1);
        
        % Lineal Intercept Data
        E112_13.li.Mu = TotLine / sum(NSegs);
        E112_13.li.LineLength = TotLine;
        E112_13.li.Sigma = std(li);
        E112_13.li.StErr = std(li) / sqrt(Sz);
        E112_13.li.CI =...
            bsxfun(@times, std(li)/sqrt(Sz), CI95li(:));
        
        % ASTM Grain Size Data
        E112_13.G.Nums = ASTM_Nums;
        E112_13.G.Mu = mean(ASTM_Nums);
        E112_13.G.StDev = std(ASTM_Nums);
        E112_13.G.ASTM.StError = std(ASTM_Nums) / sqrt(Sz);
        E112_13.G.ASTM.CI =...
            bsxfun(@times, std(ASTM_Nums)/sqrt(Sz), CI95ASTM(:));
        
        % Surface Volume Data
        E112_13.Sv.Mu = mean(Sv);
        E112_13.Sv.Sigma = std(Sv);
        E112_13.Sv.StError = std(ASTM_Nums) / sqrt(Sz);
        E112_13.Sv.CI =...
            bsxfun(@times, std(ASTM_Nums)/sqrt(Sz), CI95Sv(:));
        
        % Save the number of segmentations
        E112_13.Ns = NSegs;
        
        % Return the grain object
        varargout{2} = E112_13;
    end
    % First structure is always the myEBSD structure, second (if grains) is
    % a structure that will be added to the main Grain Structure
    varargout{1} = myEBSD;
end

