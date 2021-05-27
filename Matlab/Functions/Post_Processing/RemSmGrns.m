function [Ebsd,grains,SmllGrns] = RemSmGrns(Ebsd,grains,PixTol)
%%                      Function Description

% Function that removes grains smaller than a user-defined tolerance and
% rearranges the grain structure such that only grains that surpass the
% size limit are considered standalone grains; all other grain boundaries
% are removed entirely and these grains are stored separately

% Input Variables: grain structure; Ebsd structure; pixel tolerance for
%                  grain size acceptance

% Output Variables: updated grain structure, updated Ebsd structure, set of
% smaller grains under pixel tolerance


%%                      Function Implementation

    % First fine computed grains smaller than the pixel size tolerance
    SmllGrns = find(grains.grainSize<PixTol);
    
    % If any of the grains are flagged as being too small, we need to cycle
    % through and reindex everything for proper bookkeeping
    if isempty(SmllGrns)==0
        % First, identify the larger grains (that we'll keep) and their
        % corresponding grain ids
        LrgGrns = grains.grainSize>=PixTol;
        tmpGrnIds = Ebsd.grainId;
        tmpGBedgs = grains.boundary.grainId;
        LrgGrnIds = grains.id(LrgGrns);
        FullGrnCnt = 0;

        % Now reassign the id values from 1 to the number of ``large''
        % grains, and increase our counter
        for ij = 1:length(LrgGrnIds)
            FullGrnCnt = FullGrnCnt+1;
            % Reindex grain ids and the boundaries connecting grains
            tmpGrnIds(Ebsd.grainId==LrgGrnIds(ij)) = FullGrnCnt;
            GBedg = ismember(tmpGBedgs,LrgGrnIds(ij));
            grains.boundary.grainId(GBedg) = FullGrnCnt;
%                 tmpGBedgs(find(tmpGBedgs==LrgGrnIds(ij))) = FullGrnCnt;

        end

        % Now do the same for small grains so that the ids assigned to the
        % ebsd data set are consistent
        SmllGrnIds = grains.id(SmllGrns);
        for ij = 1:length(SmllGrnIds)
            FullGrnCnt = FullGrnCnt+1;
            % Reindex grain ids and the boundaries connecting grains
            tmpGrnIds(Ebsd.grainId==SmllGrnIds(ij)) = FullGrnCnt;
            GBedg = ismember(tmpGBedgs,SmllGrnIds(ij));
            grains.boundary.grainId(GBedg) = FullGrnCnt;
        end
        % Now arrange the grains based on first a set of the large
        % grains and next the set of smaller grains that will be
        % omitted from subsequent grain calculations
        TotGrnIds = [LrgGrnIds;SmllGrnIds];
        grains = grains(TotGrnIds);
        grains.id = [1:FullGrnCnt]';
        SmllGrns = find(grains.grainSize<PixTol);
        grains(SmllGrns)=[];  
        Ebsd.grainId = tmpGrnIds;
    end


end

