%clc
%q1=Q;
%q2=QQ;
%qsym=QS;


function [th, ax] = QuatMis(q1, q2, qsym)
% Misorientation - Return misorientation data for quaternions.
%
%   USAGE:
%
%   angle = Misorientation(q1, q2, sym)
%   [angle, mis] = Misorientation(q1, q2, sym)
%
%   INPUT:
%
%   q1 is 4 x n1, 
%      is either a single quaternion or a list of n quaternions
%   q2 is 4 x n,  
%      a list of quaternions
% 
%   OUTPUT:
%
%   angle is 1 x n, 
%         the list of misorientation angles between q2 and q1
%   axis  is 3 x n list of 
%   mis   is 4 x n, (optional) 
%         is a list of misorientations in the fundamental region 
%         (there are many equivalent choices)
%
%   NOTES:
%
%   *  The misorientation is the linear tranformation which
%      takes the crystal basis given by q1 to that given by
%      q2.  The matrix of this transformation is the same
%      in either crystal basis, and that is what is returned
%      (as a quaternion).  The result is inverse(q1) * q2.
%      In the sample reference frame, the result would be
%      q2 * inverse(q1).  With symmetries, the result is put
%      in the fundamental region, but not into the Mackenzie cell.
%

% Get sizes of quaternion arrays
[n1,qchk1]=size(q1);[n2,qchk2]=size(q2);

% Throw errors if not correct input formats
if qchk1~=4,error('q1 is not a n x 4 list of quaternions.');end;
if qchk2~=4,error('q2 is not a n x 4 list of quaternions.');end;

% We allow three options: n1 and n2 are equal, n1 is unity, or n2 is unity
if n1~=n2 && (n1~=1 || n2~=1) %EJP! This is not logically correct
    error(['If q1 or q2 do not contain the same number of quaternions,' ...
        ' one of them must contain only one quaternion.']);
end
if (n1==1), q1=repmat(q1,[n2 1]);end
if (n2==1), q2=repmat(q2,[n1 1]);end

% Get the quaternion that describes the misorientation
q1i=[-q1(:,1) q1(:,2:4)];

mis=Quat2FundZone(QuatProd(q2,q1i), qsym);
% Extract the misorientation angle and axis
[th,ax]=Quat2AngAx(mis);


