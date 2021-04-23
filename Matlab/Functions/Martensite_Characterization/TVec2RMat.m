function T=TVec2RMat(V)
% rvec2matcel converts a cell structure of 3x3 matrices into a matrix where
% the row vectors represent the matrices of the cell structure. Useful for
% storing groupings of symmetry matrices, for example.
%
% SEE ALSO: matcel2rvec (available from current author)
%
% -------------------------------------------------------------------------
% 2010-11-09 | Eric Payton, Ruhr-Universitaet Bochum (payton.28[at]osu.edu)
% -------------------------------------------------------------------------
% This program is provided without any guarantee of correctness.
% If you modify it and/or improve it, please kindly share with me your new
% and improved version to the email address above. Thanks!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,xxx]=size(V);clear xxx;
for i=1:n
    T{i,1}=[V(i,1:3);V(i,4:6);V(i,7:9)];
end