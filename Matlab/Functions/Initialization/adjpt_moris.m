function [adj_pts,Mori,flag] = adjpt_moris(varargin)
% compute adjacent measurements

% Create adjacency array based on specific interactions between points
% within the input grain

Ebsd = varargin{1};
nn = varargin{2};

len_v = length(varargin);

if length(Ebsd) < 10
    Dl = [];
    Dr = [];
else
    [~,~,I_FD] = spatialDecomposition([Ebsd.prop.x(:), Ebsd.prop.y(:)],Ebsd.unitCell,'unitCell');
    A_D = I_FD.' * I_FD;

    A_D1 = A_D;
    for i = 1:nn-1
        A_D = A_D + A_D*A_D1 + A_D1*A_D;
    end
    clear A_D1

    % extract adjacent pairs
    [Dl, Dr] = find(A_D);

    % take only ordered pairs of same, indexed phase
    use = Dl > Dr & Ebsd.phaseId(Dl) == Ebsd.phaseId(Dr) & Ebsd.isIndexed(Dl);

    Dl = Dl(use); Dr = Dr(use);

    if nn == 2 && len_v == 2
       Xlen = length(unique(Ebsd.x));
       del_x2 = Dr+2;
       del_y2 = (Xlen*2+Dr);
       del = find(Dl == del_x2 | Dl == del_y2);
       Dl(del) = [];
       Dr(del) = [];
    end

    phaseId = Ebsd.phaseId(Dl);


    for p=1:numel(Ebsd.phaseMap)


        currentPhase = phaseId == p;
        if any(currentPhase)

            o_Dl = orientation(Ebsd.rotations(Dl(currentPhase)),Ebsd.CSList{p});
            o_Dr = orientation(Ebsd.rotations(Dr(currentPhase)),Ebsd.CSList{p});
            %omega(currentPhase) = angle(o_Dl,o_Dr);

        end
    end
end

% If everything is empty, it prompts the termination flag
if isempty(Dl) || isempty(Dr)
    flag = 1;
    Mori = [];
    adj_pts = [];
else
    flag = 0;
    Mori=inv(o_Dl).*(o_Dr);
    adj_pts = [Dl,Dr];
    
end

end

