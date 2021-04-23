function [V,varargout]=YardleyVariants(OR)
% YardleyVariants returns the matrices corresponding to the variants
%                 produced from the provided orientation relationship,
%                 specified in Kurjumov-Sachs angles.  
%
% In the case of the Kurdjumov-Sachs or Nishiyama-Wasserman orientation
% relationships, 'KS' or 'NW' can be passed to the function as OR.
%
% -------------------------------------------------------------------------
% 2010-11-25 | Victoria Yardley (victoria.yardley[at]rub.de)
%              Eric Payton (payton.28[at]osu.edu)
%              Ruhr-Universitaet Bochum
%
% %%%         Matlab function written by EP based on VY's             %%% %
% %%%           spreadsheet for calculation of variants               %%% %
%
% -------------------------------------------------------------------------
% This program is provided without any guarantee of correctness.
% If you modify it and/or improve it, please kindly share with me your new
% and improved version to the email address above. Thanks!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag=0;
%% Parse input orientation relationship
if ~isnumeric(OR)
    if OR=='KS'
        s6=sqrt(6.0);s3=sqrt(3.0);
        ksi1=real(acos((s6+1.0)/(2.0*s3)));
        ksi2=real(acos((s6+18.0)/(12.0*s3)));
        ksi3=real(acos((s6+12.0)/(6.0*s6)));
        OR=[ksi1 ksi2 ksi3];clear ksi1 ksi2 ksi3 s6 s3
    elseif OR=='NW'
        s6=sqrt(6);s2=sqrt(2);
        ksi0=real(acos((s2+1.0)/s6));
        OR=[0.0 ksi0 ksi0];clear ksi0 s6 s2
    else
        error('Unknown named orientation relationship')
    end
else
   % Convert OR specification into radians
   OR=OR.*(pi/180);   
end

%% Get the misorientation of the OR from the 1st Bain correspondence matrix
MB=cell(2,1); MB{1}=zeros([3 3]);
MB{1}(1,1)=cos(OR(1)); MB{1}(2,2)=cos(OR(2)); MB{1}(3,3)=cos(OR(3));
costh=0.5*(trace(MB{1})-1.0); mosth=1-costh;sinth=sqrt(1-costh^2);
r1=(sqrt((MB{1}(1,1)-costh)/(mosth)));
r2=(sqrt((MB{1}(2,2)-costh)/(mosth)));
r3=(sqrt((MB{1}(3,3)-costh)/(mosth)));
clear costh OR
r1r2=r1*r2*mosth; r1r3=r1*r3*mosth; r2r3=r2*r3*mosth;
r3st=r3*sinth;    r2st=r2*sinth;    r1st=r1*sinth;
clear r1 r2 r3 mosth sinth
MB{1}(2,3)= r2r3-r1st; MB{1}(3,2)= r2r3+r1st;
MB{2}=MB{1};
MB{1}(1,2)=-r1r2+r3st; MB{1}(1,3)=-r1r3-r2st;
MB{1}(2,1)=-r1r2-r3st; MB{1}(3,1)=-r1r3+r2st;
clear r1r2 r1r3 r2r3 r3st r2st r1st
MB{2}(1,2)=-MB{1}(1,2); MB{2}(1,3)=-MB{1}(1,3);
MB{2}(2,1)=-MB{1}(2,1); MB{2}(3,1)=-MB{1}(3,1);
% MB{1} is the positive solution; MB{2} is the negative solution

%% Bain correspondence matrices
B{ 1}=[ 1 -1  0; 1  1  0; 0  0  1];
B{ 2}=[ 0  1 -1; 0  1  1; 1  0  0];
B{ 3}=[-1  0  1; 1  0  1; 0  1  0];
B{ 4}=[ 0  1  1; 0 -1  1; 1  0  0];
B{ 5}=[-1 -1  0; 1 -1  0; 0  0  1];
B{ 6}=[ 1  0 -1; 1  0  1; 0 -1  0];
B{ 7}=[ 1  1  0;-1  1  0; 0  0  1];
B{ 8}=[-1  0 -1;-1  0  1; 0  1  0];
B{ 9}=[ 0 -1  1; 0  1  1;-1  0  0];
B{10}=[ 1  0  1; 1  0 -1; 0  1  0];		
B{11}=[ 0 -1 -1; 0  1 -1; 1  0  0];
B{12}=[-1  1  0; 1  1  0; 0  0 -1];
% Normalize correspondence matrices
for i=1:length(B)
    for j=1:3, B{i}(:,j)=B{i}(:,j)/norm(B{i}(:,j));end
end

%% Produce variants
j=1;for i=1:length(B),V{j}=MB{1}*B{i};j=j+1;V{j}=MB{2}*B{i};j=j+1;end
clear B MB i j

%% Reduce redundancies, if they exist (for example, as in NW)
T=RMat2TVec(V);clear V
[~,idx]=unique(sigdec(T,7),'rows','first');
T=T(sort(idx),:); % don't allow reordering!

%% Check if results are valid
if isreal(T)
    V=TVec2RMat(T);
else
    %warning('Ksi values produce some imaginary numbers.')
    V=TVec2RMat(T);
    flag=1;
end

varargout{1}=flag;
%% Notes
% For the 'true' KS first variant, we should have:
%T1KS=(1/(6*sqrt(6)))* ...
%    [2*(sqrt(6)+3) -4*sqrt(6) 2*(sqrt(6)-3); ...
%    12-sqrt(6) 2*(sqrt(6)+3) -sqrt(6); ...
%    sqrt(6) -2*(sqrt(6)-3) 12+sqrt(6)];
% On my machine, the above algorithm produces a result with a max round-off
% error of 2.5E-15 for this case.
%
% For the 'true' NW third variant, we should have
%T1KS=(1/(6*sqrt(2)))* ...
%    [0                  6                   -6; ...
%    sqrt(6)*(sqrt(2)-2) sqrt(6)*(sqrt(2)+1) sqrt(6)*(sqrt(2)+1);...
%    sqrt(6)*(sqrt(2)+2) sqrt(6)*(sqrt(2)-1) sqrt(6)*(sqrt(2)-1)];
% On my machine, the above algorithm produces a result with a max round-off
% error of 3.3E-16 for this case.