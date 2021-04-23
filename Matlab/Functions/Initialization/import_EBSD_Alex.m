function [myEBSD] = import_EBSD(fname,myEBSD)

    %% Import the Data
    
    % Import EBSD data with an accurate description for each specific phase
    % based on how many exist in the post-transformation microstructure.
    % Assign crystal symmetry to each phase
    
    % Now set up crsytal symmetry with a first entry that relates to
    % non-indexed points
    myEBSD.CS{1} = 'notIndex';
    Material = myEBSD.Material;
    
    % Figure out the file type
    if isempty(strfind(fname,'.ang'))==0
        % Create an EBSD variable containing the data
        ftype = 'ang';
        ebsd = loadEBSD(fname,{},'interface',ftype,'convertEuler2SpatialReferenceFrame');
        ebsd=ebsd.fill;
        unIndpts = [];
    elseif isempty(strfind(fname,'.ctf'))==0
        ebsd = loadEBSD_ctf(fname,'convertEuler2SpatialReferenceFrame');
        ftype = 'ctf';
        % Find unindexed points
        unIndpts = find(ebsd.isIndexed==0);
    end
    
    % Number of unique phases (delete if unindexed)
    phases = unique(ebsd.phase);
    nph = length(phases);
        if isempty(unIndpts)==0
            unIndPh = ebsd(unIndpts).phase(1);
            phases(find(phases==unIndPh)) = [];
            nph = nph+1;
        end
    num_phases = length(phases);
    
    % Unit Cell Dimensions
    dims1 = [myEBSD.Celldims{1}];
    dims2 = [myEBSD.Celldims{2}];
    
    % Pre-transformation symmetry (bcc/fcc for Beta-Ti/Gamma-Steel)
    SymPre = 'm-3m';
    
    if strcmp(Material,'Steel')
        SymPost = 'm-3m';
    elseif strcmp(Material,'Titanium')
%         SymPost = '622';
        SymPost = '6/mmm';
    else
        Error('Algorithm can only handle Mart to Aus or Alpha-Ti to Beta-Ti Transformations')
    end
    
    % Define post-transformation crystal symmetry
    CS_post=crystalSymmetry(SymPost,dims1,'mineral',myEBSD.Phase.Name{1});
    CS_post.i=CS_post.i(:,1);
    CS_post.a=CS_post.a(:,1);
    CS_post.b=CS_post.b(:,1);
    CS_post.c=CS_post.c(:,1);
    CS_post.d=CS_post.d(:,1);

    % Pre-transformation crystal symmetry
    CS_pre=crystalSymmetry('m-3m',dims2,'mineral',myEBSD.Phase.Name{2});
    CS_pre.i=CS_pre.i(:,1);
    CS_pre.a=CS_pre.a(:,1);
    CS_pre.b=CS_pre.b(:,1);
    CS_pre.c=CS_pre.c(:,1);
    CS_pre.d=CS_pre.d(:,1);

    if num_phases > 1
        for i = 1:2
            len_ph(i) = length(find(ebsd.phase==phases(i)));
        end
        [~,posttrans] = max(len_ph);
        posttrans = posttrans+1;
        
        if posttrans == 2
            pretrans = 3;
        else
            pretrans = 2;
        end
    else
        posttrans=2;
        pretrans=3;
    end
    
    % Assign crystal symmetries
    myEBSD.CS{posttrans} = CS_post;
    myEBSD.CS{pretrans} = CS_pre;
    for ii = pretrans+1:nph
            myEBSD.Phase.Name{3} = 'Arbitrary';
            CS_arb=crystalSymmetry(SymPost,dims1,'mineral','Arbitrary');
            CS_arb.i=CS_arb.i(:,1);
            CS_arb.a=CS_arb.a(:,1);
            CS_arb.b=CS_arb.b(:,1);
            CS_arb.c=CS_arb.c(:,1);
            CS_arb.d=CS_arb.d(:,1);
            myEBSD.CS{ii} = CS_arb;
    end
    
    clear ebsd
    
    % Additional post transformation phases will be considered the
    % pre-transformation phase in terms of the reconstruction
%     if num_phases > 2
%         fprintf('Type in phase name surrounded by single quotations \n')
%         myEBSD.Ph_names{3} = input('Phase')
%         fprintf('Type in cell dimension surrounded by single quotations \n')
%         dim = prompt('Cell Dimension')
%         dim = str2num(dim)
%         Sym = prompt('Symmetry')
%         myEBSD.Celldims{3} = dim
%         dim = myEBSD.Celldims{3};
%         CS_post2=crystalSymmetry(Sym, [dim dim dim], 'mineral', myEBSD.Phase.Name{3});
%         CS_post2.i=CS_post2.i(:,1);
%         CS_post2.a=CS_post2.a(:,1);
%         CS_post2.b=CS_post2.b(:,1);
%         CS_post2.c=CS_post2.c(:,1);
%         CS_post2.d=CS_post2.d(:,1);
%         myEBSD.CS{4} = CS_post2;
%     end
%     if num_phases > 3
%         warning('Max number of phases present, any more will result in an error')
%         fprintf('Type in phase name surrounded by single quotations \n')
%         myEBSD.Ph_names{4} = input('Phase')
%         fprintf('Type in cell dimension surrounded by single quotations \n')
%         dim = prompt('Cell Dimension')
%         dim = str2num(dim)
%         Sym = prompt('Symmetry')
%         myEBSD.Celldims{4} = dim
%         CS_post3=crystalSymmetry(Sym, [dim dim dim], 'mineral',myEBSD.Phase.Name{4});
%         CS_post3.i=CS_post3.i(:,1);
%         CS_post3.a=CS_post3.a(:,1);
%         CS_post3.b=CS_post3.b(:,1);
%         CS_post3.c=CS_post3.c(:,1);
%         CS_post3.d=CS_post3.d(:,1);
%         myEBSD.CS{5} = CS_post3;
%     end
%     
%     if num_phases > 4
%         error('Too many phases present in microstructure')
%     end
    % Create an EBSD variable containing the data
%     if isempty(strfind(fname,'.ang'))==0
        % Create an EBSD variable containing the data
    ebsd = loadEBSD(fname,myEBSD.CS,'interface',ftype,'convertEuler2SpatialReferenceFrame');
    ebsd=ebsd.fill;
    if num_phases > 2
        tmpEb = ebsd;
        nontranspts = find(ebsd.phase~=(posttrans-1) & ebsd.phase~=(pretrans-1));
        tmpEb(nontranspts).phase=max(posttrans,pretrans);
        ebsd(tmpEb.id).phase = tmpEb.phase;
        if isempty(unIndpts)==0
            ebsd(unIndpts).phase=0;
        end
    end
    if isempty(strfind(fname,'.ang'))==0
        myEBSD.ci = ebsd.ci;
    elseif isempty(strfind(fname,'.ctf'))==0
        myEBSD.ci = ebsd.mad;
    end
    myEBSD.ci(myEBSD.ci<=0)=1e-6;

    if num_phases > 1
        for j = 1:length(phases)
            ph_Id = ebsd(myEBSD.Phase.Name{j}).phase(1);
            if isempty(ph_Id)==1
                myEBSD.Phase.ID{j}=[];
            else
                myEBSD.Phase.ID{j} = ph_Id;
            end
            myEBSD.Phase.Index{j} = j;
        end
    else
        init_ph = ebsd(ebsd.isIndexed).phase(1);
        myEBSD.Phase.ID{1} = init_ph;
        myEBSD.Phase.ID{2} = init_ph+1;
        myEBSD.Phase.Index{1} = 1;
        myEBSD.Phase.Index{2} = 2;
    end

    % Re-assign crystal symmetry to align with order of phase names
    myEBSD.CS{2} = CS_post;
    myEBSD.CS{3} = CS_pre;
    
    % include triclinic sample symmetry for completeness
    SS=specimenSymmetry('triclinic');
    SS.i=SS.i(:,1);
    SS.a=SS.a(:,1);
    SS.b=SS.b(:,1);
    SS.c=SS.c(:,1);
    SS.d=SS.d(:,1);
    
    myEBSD.SS = SS;
    myEBSD.Ebsd = ebsd;
    myEBSD.Phase.n_phases = num_phases;
    myEBSD = rmfield(myEBSD,'Celldims');
    myEBSD.origEbsd = ebsd;
   
end