function PltPairedVars(varargin)
    % Plot even and odd variants related to the differing rotations of 
    % the Bain correspondence axes
    
    % Var ID
    myEBSD = varargin{1};
    
    if length(varargin)==2
        mbar = 'MicronBar';
        mkey = 'off';
    else
        mbar = [];
        mkey = [];
    end
    Vars = myEBSD.Variants.IDs;
    % Phase IDs
    TransID = myEBSD.Phase.ID{1};
    ReconID = myEBSD.Phase.ID{2};
    
    InitAus = find(myEBSD.Ebsd.phase==ReconID);
    Vars(InitAus) = [];
    
    Martensite = myEBSD.Ebsd(find(myEBSD.Ebsd.phase==TransID));
    RecEbsd = myEBSD.Recon.Ebsd;
    
    key = ipfHSVKey(Martensite);
    key.inversePoleFigureDirection=zvector;
    mc = key.orientation2color(Martensite.orientations);
    VarZero = find(Vars==0);
    VarEven = find(mod(Vars,2)==0);
    VarOdd = find(mod(Vars,2)>0);
    
    figure;
    plot(RecEbsd(VarEven),mc(VarEven),'FaceColor','green',mbar,mkey)
    hold on
    plot(RecEbsd(VarOdd),mc(VarOdd),'FaceColor','blue')
    hold on
    plot(RecEbsd(VarZero),mc(VarZero),'FaceColor','black')
    hold on
    
    % Plot block boundaries
    BlckBounds = myEBSD.Blocks.Boundaries;
    
    hold on
    plot(BlckBounds.boundary,'FaceColor','black')
end

