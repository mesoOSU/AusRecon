function thc=CSLcrit(varargin)
% Calculates the angular deviation from the ideal misorientation
% relationship allowed by different criteria of the form
%               \theta_c=\theta_0 / sigma^eta
% 
% INPUTS: An array of CSL values, theta_0, eta, and/or the criterion name.
%         Allowed names: Brandon, Pumphrey, Palumbo-Aust, and
%         Ishida-McLean. If the array argument is left out, then the
%         function returns a list for all CSL between 1 and 51.
%
%         If specifying CSL values, theta_0 and eta, they must be given in
%         exactly that order.
%
%         If no input arguments are provided, the Brandon criterion is
%         used.
%
% OUTPUTS: List of maximum deviation angles for the provided CSL values
%          under the specified criterion. Angle is returned in radians!
%
% EXAMPLE USAGE: thc=CSLcrit([3 9 27],'Brandon');
%                thc=CSLcrit([3 9 27 51],15,0.5)
%
% REFERENCES:
% [1] Warrington & Boon, Acta Metall 23 (1975) p599
% [2] Brandon, Acta Metall 14 (1966) p1479
% [3] Pumphrey, in Grain Boundary Structure & Properties, London: Academic
% Press (1976) p139
% [4] Palumbo & Aust, Acta Metall 38 (1990) p2343
% [5] Ishida & McLean, Philos Mag 27 (1973) p1125
%
% -------------------------------------------------------------------------
% 2011-05-02 | Eric Payton, Ruhr-Universitaet Bochum (payton.28[at]osu.edu)
% -------------------------------------------------------------------------
% This program is provided without any guarantee of correctness.
% If you modify it and/or improve it, please kindly share with me your new
% and improved version to the email address above. Thanks!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name=[];
if isempty(varargin)
    sig=1:2:51;name='Brandon';
else % parse input arguments
    if length(varargin)==1 % then it can be a list of CSL values or a name
        if isnumeric(varargin{1})
            sig=varargin{1};name='Brandon';
        elseif ischar(varargin{1})
            sig=1:2:51;name=varargin{1};
        else
            error('CSLcrit: check varargins');
        end
    elseif length(varargin)==2 % list+name OR theta_0+eta
        if isnumeric(varargin{1}) && isnumeric(varargin{2})
            sig=1:2:51;th0=varargin{1};eta=varargin{2};
        elseif isnumeric(varargin{1}) && ischar(varargin{2})
            sig=varargin{1};name=varargin{2};
        elseif isnumeric(varargin{2}) && ischar(varargin{1})
            sig=varargin{2};name=varargin{1};
        else
            error('CSLcrit: check varargins');
        end
    elseif length(varargin)==3 % everything is specified
        % check that arguments are in the correct order
        if length(varargin{2})>1 || length(varargin{3})>1
            error('CSLcrit: check varargins')
        else % assign values from varargins
            sig=varargin{1};th0=varargin{2};eta=varargin{3};
        end
    end
end
if ~isempty(name)
    if strcmp(name,'Brandon') || strcmp(name,'brandon') || strcmp(name,'B')
        th0=15;eta=1/2;
    elseif strcmp(name,'Pumphrey') || strcmp(name,'pumphrey') ...
            || strcmp(name,'P')
        th0=15;eta=2/3;
    elseif strcmp(name,'Palumbo-Aust') || strcmp(name,'palumbo-aust')...
            || strcmp(name,'PA')
        th0=15;eta=5/6;
    elseif strcmp(name,'Ishida-McLean') || strcmp(name,'ishida-mclean')...
            || strcmp(name,'Ishida-Mclean') || strcmp(name,'IM')
        th0=8;eta=1;
    else
        error('CSLcrit: unrecognized named criterion')
    end
end

th0=th0*pi/180;

% Calculate criteria
thc=th0./sig.^eta;         
