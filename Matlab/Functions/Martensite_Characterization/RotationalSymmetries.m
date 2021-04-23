function [S,PGN]=RotationalSymmetries(varargin)
% RotationalSymmetries returns the proper rotation symmetry elements of the
%                      supplied point group.
%
% REQUIRES:
% PointGroupElements, sigdec (both available from the current author)
%
% -------------------------------------------------------------------------
% 2010-11-26 | Eric Payton, Ruhr-Universitaet Bochum (payton.28[at]osu.edu)
% -------------------------------------------------------------------------
% This program is provided without any guarantee of correctness.
% If you modify it and/or improve it, please kindly share with me your new
% and improved version to the email address above. Thanks!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(varargin)==1
    [PGE,PGN]=PointGroupElements(varargin{1});
elseif length(varargin)==2
    [PGE,PGN]=PointGroupElements(varargin{1},varargin{2});
end

if isempty(PGE)
    error(['Invalid point group name. '...
        'See PointGroupElements documentation.'])
end

% Proper rotations have a determinant of 1; rotoinversions have a
% determinant of -1
j=1;n=length(PGE);
for i=1:n
    if sigdec(det(PGE{i}),5)==1 % in case determinant not exactly integer
        S{j}=PGE{i};j=j+1;
    end
end