function [DiGraph,end_node,sink_node] = graph_setup(DiGraph,adj_pts,IP_wts,OP_wts,ngrid)
% Set up the graph with minimal inputs to clear up space in the main code
% and keep everything as organized as possible.

Dl = adj_pts(:,1);
Dr = adj_pts(:,2);

% Add source node
DiGraph=addnode(DiGraph,1); %apple

% Add in current nodes corresponding to the length of the cut out image
DiGraph=addnode(DiGraph,length(OP_wts));

if ngrid == 1
    sourcewts = ones(size(OP_wts))*10;
else
    sourcewts = ones(size(OP_wts))*1e9;
end

% Add source weights (connecting source to the current image--uncuttable)
DiGraph=addedge(DiGraph,1,2:(length(OP_wts)+1),sourcewts);

% Add inplane weights for current
DiGraph=addedge(DiGraph,Dl+1,Dr+1,IP_wts);
DiGraph=addedge(DiGraph,Dr+1,Dl+1,IP_wts); 

if ngrid == 2
    
    % Add endnode
    end_node=size(DiGraph.Nodes,1);

    % Add nodes for guess
    DiGraph=addnode(DiGraph,length(OP_wts));

    % Add in-plane weights for guess
    DiGraph=addedge(DiGraph,Dl+end_node,Dr+end_node,IP_wts);
    DiGraph=addedge(DiGraph,Dr+end_node,Dl+end_node,IP_wts);

    % Add out-of-plane weights from the current to the guess microstructure
    DiGraph=addedge(DiGraph,(1:length(OP_wts))+1,(1:length(OP_wts))+end_node,OP_wts);

    % Add penalty weight (can't cut in opposite direction)
    DiGraph=addedge(DiGraph,(1:length(OP_wts))+end_node,(1:length(OP_wts))+1,1e9);
else
    end_node = 1;
end

% Add sink node
DiGraph=addnode(DiGraph,1);
sink_node=size(DiGraph.Nodes,1);

% Finally add preliminary 2nd-tier OP weights
DiGraph=addedge(DiGraph,(1:length(OP_wts))+end_node,sink_node,1);
end
