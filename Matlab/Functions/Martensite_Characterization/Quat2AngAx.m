function [th,ax]=Quat2AngAx(Q)

th=2.*acos(min(1,Q(:,1)));

sqi=1./sin(th./2.);
vt=[Q(:,2).*sqi Q(:,3).*sqi Q(:,4).*sqi];

loc=th>pi;
th(loc)=2.*pi-th(loc);
vt(loc,:)=-vt(loc,:);

fix=th<eps('single');
if ~isempty(fix)
    th(fix)=0;
    vt(fix,1)=1;
    vt(fix,2)=0;
    vt(fix,3)=0;
end

fix=abs(pi-th)<eps('single');
if ~isempty(fix)
    th(fix)=pi;
end

% normalize the vector
vt2=1./sum(vt.^2,2);
ax=[vt(:,1).*vt2 vt(:,2).*vt2 vt(:,3).*vt2];