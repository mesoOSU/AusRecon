function [Twin,TCount] = TwnAssign(Twin,TCount,Ebsd,Inds,POr,TOrs,TID)
    % Assign relevant data to the Twin Structure if grain was deemed to
    % contain a twin

    % Increase the twin counter
    TCount = TCount+1;
    
    % Keep track of the merged parent-twin system
    Twin.Merged{TCount,1} = Ebsd(Inds).id;
    Twin.Parent.Or{TCount,1} = POr;

    % Keep track of individual twin indices within the parent
    % grain along with orientations and update the likelihoods
    for j = 1:length(TID)
        Twin.Single.Or{TCount,j} = TOrs(j);
        Twin.Single.ID{TCount,j} = TID(j); 
    end
    Twin.Flg{TCount} = 1;

end

