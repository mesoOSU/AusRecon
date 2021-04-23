function V=RMat2TVec(T)
% matcel2rvec converts a matrix where the row vectors represent 3x3
% matrices into a cell structure of 3x3 matrices. Useful for storing
% groupings of symmetry matrices, for example.
%
% SEE ALSO: rvec2matcel (available from current author)
%
% -------------------------------------------------------------------------
% 2010-11-09 | Eric Payton, Ruhr-Universitaet Bochum (payton.28[at]osu.edu)
% -------------------------------------------------------------------------
% This program is provided without any guarantee of correctness.
% If you modify it and/or improve it, please kindly share with me your new
% and improved version to the email address above. Thanks!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(T)
    V(i,:)=[T{i}(1,1:3) T{i}(2,1:3) T{i}(3,1:3)];
end