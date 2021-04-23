function q = QuatProd(q1, q2)

q=q1;

a1 = q1(:,1); b1 = q1(:,2); c1 = q1(:,3); d1 = q1(:,4);
a2 = q2(:,1); b2 = q2(:,2); c2 = q2(:,3); d2 = q2(:,4);

a = a1 .* a2 - b1 .* b2 - c1 .* c2 - d1 .* d2;
b = a1 .* b2 + b1 .* a2 + c1 .* d2 - d1 .* c2;
c = a1 .* c2 + c1 .* a2 + d1 .* b2 - b1 .* d2;
d = a1 .* d2 + d1 .* a2 + b1 .* c2 - c1 .* b2;

q(:,1) = a; q(:,2) = b; q(:,3) = c; q(:,4) = d;

% Make quaternion rotation always positive
loc=q(:,1)<0; q(loc,:)=-q(loc,:);