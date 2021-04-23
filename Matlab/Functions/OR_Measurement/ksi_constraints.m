function [c ceq]=ksi_constraints(samples);

ksi1=samples(1)*degree;
ksi2=samples(2)*degree;
ksi3=samples(3)*degree;

c=zeros(7,1);
c(1)=ksi1-ksi2;
c(2)=(cos(ksi2)+cos(ksi3))-(cos(ksi1)+1);
c(3)=(cos(ksi1)+cos(ksi3))-(cos(ksi2)+1);
c(4)=(cos(ksi1)+cos(ksi2))-(cos(ksi3)+1);
c(5)=cos(ksi1)+cos(ksi2)+cos(ksi3)-3;
c(6)=-(cos(ksi1)+cos(ksi2)+cos(ksi3)+1);
c(7)=abs(ksi2-ksi3)-2;

ceq=[];