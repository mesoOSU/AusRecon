function pt = tstat3( v, tp, stat )
% TSTAT3 computes one of three t_statistics: one-sided t-probability, or
%        two-sided t-probability, or the inverse t-statistic in any single 
%        call.  It does not take vectors as arguments.  
% 
%   INPUT ARGUMENTS: 
%       v       Degrees-of freedom (integer)
% 
%       tp      T-statistic: to produce t-probability
%                   OR
%               Probability (alpha): to produce t-statistic
% 
%       stat    Desired test ? 
%                   'one'   One-tailed t-probability
%                   'two'   Two-tailed t-probability
%                   'inv'   Inverse t-test
% 
%   OUTPUT ARGUMENT:
%       pt      T-probability OR T-statistic, as requested in ?stat?
% 
%   USE:
%       FIND ONE-TAILED PROBABILITY GIVEN t-STATISTIC & DEGREES-OF-FREEDOM
%           p = tstat3(v, t, 'one')
% 
%       FIND TWO-TAILED PROBABILITY GIVEN t-STATISTIC & DEGREES-OF-FREEDOM
%           p = tstat3(v, t, 'two')
% 
%       FIND ONE-TAILED t-STATISTIC GIVEN PROBABILITY & DEGREES-OF-FREEDOM
%           t = tstat3(v, p, 'inv')
%       
%  
%  
% Star Strider ? 2016 01 24 ? 
% % % % T-DISTRIBUTIONS ? 
% Variables: 
% t: t-statistic
% v: degrees of freedom
tdist2T = @(t,v) (1-betainc(v/(v+t^2), v/2, 0.5));                              % 2-tailed t-distribution
tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;                                          % 1-tailed t-distribution
% This calculates the inverse t-distribution (parameters given the
%   probability ?alpha? and degrees of freedom ?v?: 
t_inv = @(alpha,v) fzero(@(tval) (max(alpha,(1-alpha)) - tdist1T(tval,v)), 5);  % T-Statistic Given Probability ?alpha? & Degrees-Of-Freedom ?v?
statcell = {'one' 'two' 'inv'};                                                 % Available Options
nc = cellfun(@(x)~isempty(x), regexp(statcell, stat));                          % Logical Match Array
n = find(nc);                                                                   % Convert ?nc? To Integer
if (length(v) > 1) || (length(tp) > 1)
    error('                    ?> TSTAT3 does not take vectorised inputs.')
elseif isempty(n)                                                                   % Error Check ?if? Block
    error('                    ?> The third argument must be either ''one'', ''two'', or ''inv''.')
elseif (n == 3) && ((tp < 0) || (tp > 1))
    error('                    ?> The probability for ''inv'' must be between 0 and 1.')
elseif (isempty(v) || (v <= 0))
    error('                    ?> The degrees-of-freedom (''v'') must be > 0.')
end
switch n                                                                        % Calculate Requested Statistics
    case 1
        pt = tdist1T(tp, v);
    case 2
        pt = tdist2T(tp, v);
    case 3
        pt = t_inv(tp, v);
    otherwise
        pt = NaN;
end
end
% ???????????????????????????  END: tstat3.m  ????????????????????????????
