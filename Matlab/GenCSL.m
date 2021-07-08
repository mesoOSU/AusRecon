function [CSL,SQ,SV,SN]=GenCSL(varargin)
% GenCSL generates the exact misorientations of the coincident site lattice
% theory as rotation matrices and as quaternions.
%
% INPUTS:  varargin can be a list of CSL values (as in, for example,
%          [3, 9, 27]). If left blank, it generates the CSLs up to Sigma49.
%
% OUTPUTS: the CSL in a structural format (CSL), or as lists of
%          quaternions (SQ), their corresponding sigma values (SV), and
%          their standard names (SN).
%
% REFERENCES:
% [1] H. Grimmer, W. Bollman, D. H. Warrington. "Coincident-Site Lattices
% and Complete Pattern-Shift Lattices in Cubic Crystals." Acta Cryst. A30
% (1974) p197.
% [2] D. H. Warrington, P. Bufalini. "The Coincident Site Lattice and Grain
% Boundaries." Scripta Metall. 5 (1971) p771.
%
% NOTES: Function exactly reproduces Table 1 from Grimmer et al, except for
% Sigma31b, which appears to be in error in the paper.
%
% -------------------------------------------------------------------------
% 2011-05-01 | Eric Payton, Ruhr-Universitaet Bochum (payton.28[at]osu.edu)
% -------------------------------------------------------------------------
% This program is provided without any guarantee of correctness.
% If you modify it and/or improve it, please kindly share with me your new
% and improved version to the email address above. Thanks!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin)
    N=1:2:49; % Same as Grimmer et al
else
    N=varargin{1};
end

if any(iseven(N))
    error(['GenCSL: Coincident Site Lattices' ...
        'are only defined for odd integers.'])
end

% This could probably be done more efficiently... but it works and it isn't
% time consuming enough to warrant further work.
SYM=RMat2Quat(RotationalSymmetries('m-3m'));TT=[];SQ=[];SV=[];ee=1;
for Ndex=1:1:length(N)
    S=N(Ndex);ind=1;
    for h1=S:-1:1
        for k1=0:1:h1
            l1=-sqrt(S^2-h1^2-k1^2);
            if isreal(l1) && round(l1)==l1 % must be an integer; proceed
                for k2=S:-1:1
                    for l2=0:1:k2
                        h2=-sqrt(S^2-k2^2-l2^2);
                        if isreal(h2) && round(h2)==h2 % must be an integer
                            hkl1=[h1 k1 l1];
                            hkl2=[h2 k2 l2];
                            hkl3=cross(hkl1,hkl2)/S;
                            if all(round(hkl3)==hkl3) % must be integers
                                X=[hkl1' hkl2' hkl3'];
                                % Criterion 1-2 in Grimmer et al
                                if all(sum(X.^2,2)==S^2) && (S==1 || ...
                                        ~all(reshape(X./S==eye(3),9,1)))
                                    SM{ind}=[hkl1' hkl2' hkl3']/S; % matrix
                                    ind=ind+1;
                                end % if meets criterion and not identity
                            end % if hkl3 is all integers
                        end % if h2 is an integer
                    end % for possible values of l2
                end % for possible values of k2
            end % if l1 is an integer
        end % for possible values k1
    end % for possible values h1
    
    % Remove symmetric redundancies by putting into fundamental region
    O=RMat2Quat(SM);RR=zeros(length(SM),4);
    for ind=1:length(SM)
        q=QuatProd(repmat(O(ind,:),size(SYM,1),1),SYM);
        [an,ax]=Quat2AngAx(q);
        ax=sort(abs(ax(an==min(an),:)),2,'descend');an=min(an);ax=ax(1,:);
        RR(ind,:)=[an,ax];
    end
    [~,b]=unique(sigdec(RR,5),'rows');RR=RR(b,:);
    
    % Sort the results by increasing angle in axis/angle description, such
    % that, in the event of multiple results, th.a<th.b<th.c
    [~,b]=sort(RR(:,1),'ascend');RR=RR(b,:);
    RR=AngAx2Quat(RR(:,1),RR(:,2:4));
    
    % Only allow minimum sigma solutions 
    % (otherwise sigma3 is also a sigma9 -- but not vice versa)
    [~,b]=ismember(sigdec(RR,5),sigdec(TT,5),'rows');
    if ~isempty(b),RR=RR(~b,:);end;TT=vertcat(TT,RR);
    
    % Store results in structural array of TVecs; store names for output
    if size(RR,1)>1
        for ind=1:size(RR,1)
            tmp=Quat2RMat(RR(ind,:));
            eval(['CSL.S' num2str(S) char(96+ind) ...
                '=tmp{1};']);
            SN{ee}=['\Sigma' num2str(S) char(96+ind)];ee=ee+1;
        end
    else
        tmp=Quat2RMat(RR);
        eval(['CSL.S' num2str(S) '=tmp{1};']);
        SN{ee}=['\Sigma' num2str(S)];ee=ee+1;
    end
    
    % Store results as a list of quaternions corresponding Sigma values
    SQ=vertcat(SQ,RR);
    SV=vertcat(SV,repmat(S,size(RR,1),1));
    
    clear SM O an ax RR b ind  % clean up 
end % loop over all desired N
clear Ndex S SYM TT X h1 k1 l1 h2 k2 l2 hkl1 hkl2 hkl3 q N ee % clean up

% This is it.