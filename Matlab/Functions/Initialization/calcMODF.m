function [myEBSD] = calcMODF(myEBSD)

    %% Define Crystal/Sample Symmetries and OR
    CS_T = myEBSD.CS{2};
    CS_R = myEBSD.CS{3};
    SS = myEBSD.SS;
    Material = myEBSD.Material;
    OR = myEBSD.OR;
    
    %% Reconstruction Parameters and MODF Computation
    
    halfwidth = myEBSD.noise.halfwidth;

    %measured noise
    psi=deLaValeePoussinKernel('halfwidth',halfwidth);
    kappa=psi.kappa;
    
    % Generate random martensite orientations and permutate a second set of
    % the same orientations
    TransOrs=generate_simulated_data(orientation(idquaternion,CS_R),OR,halfwidth,1000,CS_T);
    TransOrs2=TransOrs(randperm(length(TransOrs)));
    
    % Compute misorientations
    mori=inv(TransOrs)*TransOrs2;
    
    % Construct MODF and add to our structure
    warning('off','MTEX:EBSD:calcODF')
    trans_modf=calcODF(mori,'kernel',psi,'silent');
    
    % Add modf and psi value to our structure
    myEBSD.MODF = trans_modf;
    myEBSD.noise.psi = psi;

end