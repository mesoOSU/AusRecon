function [myEBSD] = DataQuadrants(myEBSD)
    % Splits the EBSD dataset up into 4, overlapping quadrants for the
    % reconstruction to be run on each one individually. This substantially 
    % increases the computational efficiency of algorithm, especially with
    % regards to extremely large datasets (>1e6 data points). If input
    % dataset is small enough (<1e5 data points), algorithm leaves dataset
    % as is and performs the reconstruction on the entire dataset.
%%
    % Extract ebsd dataset and number of data points
    Ebsd = myEBSD.Ebsd;
    
    % Phase IDs
    TransID = myEBSD.Phase.ID{1};
    Ebsd(find(Ebsd.phase~=TransID))=[];
    len = length(Ebsd);
    
    % Flag if necessary and adjust the structure to relate to this
    if len < 5e6
        QuadFlag = 1;
        myEBSD.Quad.Ebsd = {Ebsd};
        myEBSD.Quad.Overlap = [];
        myEBSD.Quad.Flag = QuadFlag;
        myEBSD.Quad.OvEbsd = [];
        myEBSD.Quad.IDs = [];
        myEBSD.Quad.Edges = [];
    else
        QuadFlag = 0;
        
        % Coordinate Range
        xmin = min(Ebsd.x);
        xmax = max(Ebsd.x);
        ymin = min(Ebsd.y);
        ymax = max(Ebsd.y);
        Rx = xmax-xmin;
        Ry = ymax-ymin;
        
        % Establish overlap
        OvX = 0.1*Rx;
        OvY = 0.1*Ry;
        
        % Establish the halfpoints and the indexing step sizes used for
        % each coordinate
        unqx = unique(Ebsd.x);
        unqy = unique(Ebsd.y);
        xhlf = unqx(round(length(unqx)/2));
        yhlf = unqy(round(length(unqy)/2));
        
        % Quadrant True Coordinates
        X = [xhlf xmax];
        Y = [yhlf ymax];
        
        % Now split the dataset into 4 quadrants
        Eb = [];
        count = 1; 
        OvQuadIds = zeros(4,2);
        ActEb = Ebsd;
        Inds = ActEb.id;
        for i = 1:2
%             if i == 1
%                 Yst = Y(i);
%                 Ynd = Y(i+1);
%             else
%                 Yst = Y(i);
%                 Ynd = Y(i+1);
%             end
% 
            for j = 1:2
%                 if j == 1
%                     Xst = X(j);
%                     Xnd = X(j+1);
%                 else
%                     Xst = X(j);
%                     Xnd = X(j+1);
%                 end
                tmpInds = find(ActEb.x<=X(j) & ActEb.y<=Y(i));
                tmpEb = ActEb(tmpInds);
                Eb{count} = tmpEb;
                
                if X(j) ~= xmax
                    tmpEdg{count,1} = find(tmpEb.x==X(j) & tmpEb.y<=Y(i));
                else
                    tmpEdg{count,1} = find(tmpEb.x==min(tmpEb.x) & tmpEb.y<=Y(i));
                end
                
                if Y(i) ~= ymax 
                    tmpEdg{count,2} = find(tmpEb.x<=X(j) & tmpEb.y==Y(i));
                else
                    tmpEdg{count,2} = find(tmpEb.y==min(tmpEb.y) & tmpEb.y<=Y(i));
                end
                
                ActEdgs{count} = tmpEdg;
                ActEb(tmpInds) = [];
   
                % Now fill in the structures with Quadrant Ebsds and
                % the overlapping coordinates
                Eb{count} = tmpEb;
                count = count+1;                
            end
        end           
            
    % Now designate the overlap coordinates in a sloppy manner
    OvXCrds = find(Ebsd.x>=(xhlf-OvX) & Ebsd.x<=(xhlf+OvX));
    OvYCrds = find(Ebsd.y>=(yhlf-OvY) & Ebsd.y<=(yhlf+OvY));
    OvCrds  = [OvXCrds;OvYCrds];
    
    % Identify the overlapping regions and the region that connects them 
    % (1 for x; 2 for 1)
    OvQuadIds = [1,2,1;1,3,2;2,4,2;3,4,1];
    
    % Index the EBSD ids for each overlapping region (let's try to clean
    % this up later on, this is bad)
    OvCrdsy = find(Ebsd.y>=(yhlf-OvY) & Ebsd.y<=(yhlf+OvY));
    OvCrdsx = find(Ebsd(OvCrdsy).x<=(xhlf+OvX));
    OvEb{1} = Ebsd(OvCrdsy(OvCrdsx));
    OvCrdsx = find(Ebsd.x>=(xhlf-OvX) & Ebsd.x<=(xhlf+OvX));
    OvCrdsy = find(Ebsd(OvCrdsx).y<=(yhlf+OvY));
    OvEb{2} = Ebsd(OvCrdsx(OvCrdsy));
    OvCrdsy = find(Ebsd.y>=(yhlf-OvY) & Ebsd.y<=(yhlf+OvY));
    OvCrdsx = find(Ebsd(OvCrdsy).x>=(xhlf-OvX));
    OvEb{3} = Ebsd(OvCrdsy(OvCrdsx));
    OvCrdsx = find(Ebsd.x>=(xhlf-OvX) & Ebsd.x<=(xhlf+OvX));
    OvCrdsy = find(Ebsd(OvCrdsx).y>=(yhlf-OvY));
    OvEb{4} = Ebsd(OvCrdsx(OvCrdsy));   
    
    % Add to structure
    myEBSD.Quad.Ebsd = Eb;
    myEBSD.Quad.Overlap = OvCrds;
    myEBSD.Quad.Flag = QuadFlag;
    myEBSD.Quad.OvEbsd = OvEb;
    myEBSD.Quad.IDs = OvQuadIds;
    myEBSD.Quad.Edges = tmpEdg;
    end
end

