function QF = Quat2FundZone(Q,S)
% ToFundamentalRegionQ - To quaternion fundamental region.
%   
%   USAGE:
%
%   q = ToFundamentalRegionQ(quat, qsym)
%
%   INPUT:
%
%   quat is 4 x n, 
%        an array of n quaternions
%   qsym is 4 x m, 
%        an array of m quaternions representing the symmetry group
%
%   OUTPUT:
%
%   q is 4 x n, the array of quaternions lying in the
%               fundamental region for the symmetry group 
%               in question
%
%   NOTES:  
%
%   *  This routine is very memory intensive since it 
%      applies all symmetries to each input quaternion.
%

% Get the number of points and number of symmetries
[n,chk1] = size(Q); [m,chk2] = size(S);

% Throw errors if not correct input formats
if chk1~=4,error('quat is not a n x 4 list of quaternions.');end;
if chk2~=4,error('qsym is not a n x 4 list of quaternions.');end;

% Compute all distances to fundamental regions
R=[1. 0. 0. 0.];% reference quaternion                                     
QR=QuatProd(S,repmat(R,[m,1]));
w=QR(:,1)*Q(:,1).';
x=QR(:,2)*Q(:,2).';
y=QR(:,3)*Q(:,3).';
z=QR(:,4)*Q(:,4).';
d=abs(w+x+y+z);

% Locate and extract fundamental region for each quaternion
[om,loc]=max(d,[],1);Q2F=S(loc,:);
Q2Fi=QuatConj(Q2F); % get quaternion conjugates

% Project to fundamental region
QF=QuatProd(Q,Q2Fi);