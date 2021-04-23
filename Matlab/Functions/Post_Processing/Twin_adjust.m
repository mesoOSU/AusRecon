function [myEBSD,Twin] = Twin_adjust(myEBSD,Twin,twinIP,twinOP)
    % Adjustment of twin interfaces with respect to the parent grains for each
    % identified twin
    transpts    = myEBSD.Recon.TransformedPoints;
    Uptwin_ebsd = myEBSD.Recon.Ebsd(transpts);
    Trans_ebsd  = myEBSD.TransEbsd(transpts);
    twinInds    = Twin.Merged;
    Par_Or      = Twin.Parent.Or;
    Twin_Or     = Twin.Single.Or;
    CS_T        = myEBSD.CS{2};
    CS_R        = myEBSD.CS{3};
    ksi         = myEBSD.OR.ksi;
    psi         = myEBSD.noise.psi;
    R2T         = myEBSD.TransMats{2};
% 
%     hexFlag = length(Uptwin_ebsd.unitCell)==6;
%     % Determine nearest neighbors based on whether the ebsd grid is square or
%     % hexagonal
%     if hexFlag == 1
%         nn = 1;
%     else
%         nn = 2;
%     end
    
nn = 2;
    % Compute variants and corresponding groupoid
    [V,flag] = YardleyVariants(ksi);
    [G,CT] = GroupoidVariantsCubic(V,[]);
    
    CT_ParTwin = CT(25:end,1:24);
    CT_PTmax = max(CT_ParTwin(:));
        
    % Single out parent and twin variant misorientations
    Gpar = G(1:16,:);
    Gtwin = G(17:CT_PTmax,:);
    parax = vector3d(Gpar(:,2),Gpar(:,3),Gpar(:,4));
    twinax = vector3d(Gtwin(:,2),Gtwin(:,3),Gtwin(:,4));
    GparEul = orientation('axis',parax,'angle',Gpar(:,1),CS_T);
    GtwinEul = orientation('axis',twinax,'angle',Gtwin(:,1),CS_T);
    
    % Compute corresponding MODFs
    par_modf = calcODF(GparEul,'kernel',psi);
    twin_modf = calcODF(GtwinEul,'kernel',psi);
    
    byangle = 0;
    
%%
    for ii = 1:length(twinInds)
        % Set up our individual ebsds for austenite and martensite, calculate
        % grains, and then allocate the parent and twin orientations as
        % individual arrays
        tmp_recon = Uptwin_ebsd(twinInds{ii});
        tmp_trans = Trans_ebsd(twinInds{ii});
        parent_Or = Par_Or{ii};
        twin_Or = Twin_Or{ii};
        
        % Adjacency array and misorientations
        [tmp_adj,tmp_mori] = adjpt_moris(tmp_trans,nn);
        
%         if byangle == 1
            tmp_grains = calcGrains(tmp_recon);

            % Create array for twin planes and directions for fcc metals
            twin_planes=(Miller(1,1,1,CS_R,'hkl'));
            twin_directions = Miller(1,1,-2,CS_R,'uvw');

            % Combine to form the system of twins
            twin_syms=slipSystem(twin_directions,twin_planes);
            twin_syms=twin_syms.symmetrise('antipodal'); 
            % Rotate parent by twin families
            twin_syms=parent_Or*twin_syms;   

            % Find trace based on individual twin systems (this will give us a
            % preferred inferface direction for each system)
            trace=twin_syms.trace;

            % Determine "Schmid" factor for each system and find the maximum value.
            % This effectively gives us the optimized twin interface direction
            % based on the individual parent orientation with respect to the
            % highlighted twin
            twin_Orax = vector3d(axis(twin_Or));
            m = twin_syms.SchmidFactor(twin_Orax);
            [m,id] = max(m);
            

    %     figure; plot(tmp_grains,'micronbar','off')
    %     hold on
    %     quiver(tmp_grains,trace(id))
%             tmp_wts = (eval(twin_modf,tmp_mori)+10);

            tmp_vec = vector3d(tmp_recon.x,tmp_recon.y,0);
            tmp_adjwts = vector3d(tmp_vec(tmp_adj(:,1))-tmp_vec(tmp_adj(:,2)),0);
            tmp_IPwts = angle(tmp_adjwts,trace(id))./degree;
            tmp_IPwts(find(isnan(tmp_IPwts)))=min(tmp_IPwts);
            tmpIP_wts = (twinIP./tmp_IPwts);
%         else
%             
%         end

        tmp_graph = digraph;

        mart_par = symmetrise(parent_Or)*R2T;
        mart_par_ODF = calcODF(mart_par,'kernel',psi);
        tmpc2g_wts = eval(mart_par_ODF,tmp_trans.orientations);

        [tmp_graph,tmp_endnode,tmp_sinknode] = graph_setup(tmp_graph,tmp_adj,tmpIP_wts,tmpc2g_wts,2);

        mart_twin = symmetrise(twin_Or)*R2T;
        mart_twin_ODF = calcODF(mart_twin,'kernel',psi);
        tmpg2t_wts = eval(mart_twin_ODF,tmp_trans.orientations).*twinOP;

        tmp_graph=rmedge(tmp_graph,(1:length(tmp_trans))+tmp_endnode,tmp_sinknode);
        tmp_graph=addedge(tmp_graph,(1:length(tmp_trans))+tmp_endnode,tmp_sinknode,tmpg2t_wts);

        [mf_c,gf,cs,ct]=maxflow(tmp_graph,1,tmp_sinknode);
        ctcopy = ct;
        ctcopy(end)=[];
        tmp_ind = ctcopy-tmp_endnode;
        tmp_recon(tmp_ind).orientations = twin_Or;
       
        Uptwin_ebsd(twinInds{ii}) = tmp_recon;
        
    end
    
    Recon_ebsd = myEBSD.Recon.Ebsd;
    Recon_ebsd(transpts) = Uptwin_ebsd;
    Twin.Ebsd = Uptwin_ebsd;
end
