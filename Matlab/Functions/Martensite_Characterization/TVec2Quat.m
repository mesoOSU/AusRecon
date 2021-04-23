function Q = TVec2Quat(V)
% A similar algorithm to MTEX

% Select the largest for best accuracy
uu=V(:,1);vv=V(:,5);ww=V(:,9);
d(:,1)=sqrt(1+uu+vv+ww);d(:,2)=sqrt(1-uu-vv+ww);
d(:,3)=sqrt(1+uu-vv-ww);d(:,4)=sqrt(1-uu+vv-ww);
[mh,j]=max(d,[],2); % maximum values and locations
mh=0.5.*mh; mi=0.25./mh; % convenience variables
Q=zeros(length(V(:,1)),4);

x=j==1; % First case, d1 is max
if any(x)
  Q(x,1)=-mh(x);
  Q(x,2)=( V(x,6)-V(x,8)).*mi(x);
  Q(x,3)=( V(x,7)-V(x,3)).*mi(x);
  Q(x,4)=( V(x,2)-V(x,4)).*mi(x);
end

x=j==2; % First case, d2 is max
if any(x)
  Q(x,4)= mh(x);
  Q(x,3)=( V(x,6)+V(x,8)).*mi(x);
  Q(x,2)=( V(x,7)+V(x,3)).*mi(x);
  Q(x,1)=-( V(x,2)-V(x,4)).*mi(x);
end

x=j==3; % First case, d3 is max
if any(x)
  Q(x,2)=mh(x);
  Q(x,1)= (-V(x,6)+V(x,8)).*mi(x);
  Q(x,4)=( V(x,7)+V(x,3)).*mi(x);
  Q(x,3)=( V(x,2)+V(x,4)).*mi(x);
end

x=j==4; % First case, d4 is max
if any(x)
  Q(x,3)=mh(x);
  Q(x,4)=( V(x,6)+V(x,8)).*mi(x);
  Q(x,1)=(-V(x,7)+V(x,3)).*mi(x);
  Q(x,2)=( V(x,2)+V(x,4)).*mi(x);
end

% Follow convention that quaternion rotation is always positive
loc=Q(:,1)<0;Q(loc,:)=-Q(loc,:);