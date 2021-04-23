% MATRIX2QUATERNION - Homogeneous matrix to quaternion
%
% Converts 4x4 homogeneous rotation matrix to quaternion
%
% Usage: Q = matrix2quaternion(T)
%
% Argument:   T - 4x4 Homogeneous transformation matrix
% Returns:    Q - a quaternion in the form [w, xi, yj, zk]
%
% See Also QUATERNION2MATRIX

% Copyright (c) 2008 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

function Q = RMat2Quat(T)

V=RMat2TVec(T);
Q=TVec2Quat(V);

   