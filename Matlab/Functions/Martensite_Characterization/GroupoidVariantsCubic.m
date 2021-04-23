%restart
%V=YardleyVariants('NW');
%S=RMat2Quat(RotationalSymmetries('m-3m'));
%varargin={};

function [G,CompTable]=GroupoidVariantsCubic(varargin)
%keep('V')
%varargin={};
    
S=RMat2Quat(RotationalSymmetries('m-3m'));
res=7;

Vp = varargin{1};
Vtotal{1} = Vp;

if length(varargin)==2
    Twin_eul(1,1) = orientation('axis',vector3d(1,1,1),'angle',60*degree);
    Twin_eul(2,1) = orientation('axis',vector3d(-1,-1,1),'angle',60*degree);
    Twin_eul(3,1) = orientation('axis',vector3d(-1,1,1),'angle',60*degree);
    Twin_eul(4,1) = orientation('axis',vector3d(1,-1,1),'angle',60*degree);
    
    for ii = 1:length(Vp)
%         Vtrans_eul(ii,1) = orientation('matrix',transpose(Vp{ii}));
        Vtrans_eul(ii,1) = orientation('matrix',Vp{ii});
    end
    
    Twins = Vtrans_eul * Twin_eul;
    sz = size(Twins);
    
    for jj = 1:sz(2)
        for kk = 1:sz(1)
            Vtwinmat{kk,1} = transpose(matrix(Twins(kk,jj)));
%             Vtwinmat{kk,1} = matrix(Twins(kk,jj));
        end
        Vtotal{jj+1,1} = Vtwinmat;
        clear Vtwinmat
    end
    Vtlen = 5;
else
    Vtlen = 1;
end

numvar = length(Vp)*Vtlen;
numsym=length(S(:,1)); % number of symmetries
numint=0.5*((numvar^2)-numvar); % min # intersection calculations needed
% Create the matrix of variant intersections
k=1;
int=zeros(numint,6);


counter = 1;
Qtotal=[];
for ii = 1:Vtlen
    V1 = Vtotal{ii};
    if ii > 1
        Q1 = QuatConj(RMat2Quat(V1));
    else
        Q1 = RMat2Quat(V1);
    end
    for jj = 1:length(Q1)
        Qtotal(counter,:) = Q1(jj,:);
        counter=counter+1;
    end
end
    
for i=1:length(Qtotal(:,1))
    for j=i:length(Qtotal(:,1)) % matrix is symmetric; diagonal is all ID quat
        [an, ax]=QuatMis(Qtotal(i,:),Qtotal(j,:),S);
        int(k,:)=[i j an ax]; % intersections
        k=k+1;
    end
end

clear an ax tmp Q k i j
%% Get symmetrically-equivalent axes for each intersection
an=abs(int(:,3));
ax=sort(abs(int(:,4:6)),2);

%symvec=[an ax];
%symvec(abs(symvec)<eps('single'))=0;

%%
%ct=size(unique(sigdec(symvec,res),'rows','first'),1); %EJP! resolution can be adjusted here!

%%
int2=zeros(size(int,1),7);
int2(:,1:2)=int(:,1:2);
int2(:,4:7)=[an ax];
int2=sigdec(int2,res);

j=0;
for i=1:size(int2,1)
    if int2(i,3)==0
        j=j+1;
        symvec2(j,1:4)=int2(i,4:7);
        a=int2(:,4)==int2(i,4);
        b=int2(:,5)==int2(i,5);
        c=int2(:,6)==int2(i,6);
        d=int2(:,7)==int2(i,7);
        loc=all([a b c d],2);
        int2(loc,3)=j;
    end
end

%% Clean up a little
clear numsym sumint mt mx n nint S

% rearrange data as table
inttab=zeros(numvar,numvar);
inttab=inttab+max(max(int2))+1;
for i=1:numint
    inttab(int2(i,2),int2(i,1))=int2(i,3)-1;;
end
% for i=1:numvar
%     inttab(i,i)=0;
% end
CompTable=inttab;

% Assign outputs
G=symvec2;

G(1,:)=[];
%G(:,1)=intmin(:,1);
%G(:,2:4)=fliplr(intmin(:,2:4));
