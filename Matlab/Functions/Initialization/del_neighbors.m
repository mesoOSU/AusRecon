function [neighbors] = del_neighbors(dtr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

x=dtr.Points(:,1);
y=dtr.Points(:,2);

attachedTriangles=vertexAttachments(dtr);

for ii=1:size(x)
    vertices=dtr.ConnectivityList(attachedTriangles{ii},:);
    neighbors{ii}=setdiff(unique(vertices),ii);
end


end

