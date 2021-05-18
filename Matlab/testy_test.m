x=1:256;
y=1:256;
[X,Y] = meshgrid(x,y);

%Problem was juxt X=X(:) made a column vector, but to get the writer to work,
% you need row vectors (AKA, the transopose
% Explicit way to make rows:
% X = reshape(X,[1,256*256]);;
% Y = reshape(Y,[1,256*256]);
%Speedy quick way to make rows
X = X(:).';
Y = Y(:).';

A = [X,Y];
fID = fopen('test.txt','wt');
fprintf(fID, '%6.0f   %6.0f \n', [X; Y]);
fclose(fID);