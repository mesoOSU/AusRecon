function [ flow ] = calc_flow(node_index,gf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

flow=sum(gf.Edges.Weight(gf.Edges.EndNodes(:,2)==node_index));
end

