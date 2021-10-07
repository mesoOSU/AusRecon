clear all
clc

% Read Image into MATLAB and make it a vector
% Image_full = imread('Image_Segmentation/Ash_gprime/2B_ORtilt_T1/New_T1/2B_7.tif');
% Image_full = imread('Image_Segmentation/Ash_gprime/2B_ORtilt_T1/Ex_Sec/2B_10.tif');
% Image_full = imread('Image_Segmentation/Ash_gprime/3A_subsol_wpri/pGam_pri/SE/Tile01.tif');
% Image_full = imread('Image_Segmentation/Ash_gprime/12A_supsol_wopri/ETD_notilt/10KX3.TIF');
% -- doesn't really work
 Image_full = imread('Image_Segmentation/Ash_gprime/12A_supsol_wopri/ETD_tilt/Tilt2.tif');

Image = im2double(Image_full(100:1000,100:1000))*255;
Image_vec = Image(:);

%%
% Gaussian Mean distribution to fit data to a bimodal distribution, and
% pull out the mean and variance from the fitted data
% part = find(Image_vec<=50);
% datapart = gmdistribution.fit(Image_vec(part),1);
% mean_part = datapart.mu;
% sigma_part = sqrt(datapart.Sigma);
% 
% % Matrix
% mat = find(Image_vec>50);
% datamat = gmdistribution.fit(Image_vec(mat),1);
% mean_mat = datamat.mu;
% sigma_mat = sqrt(datamat.Sigma);

data = gmdistribution.fit(Image_vec(:),2);
mean_mat = max(data.mu);
mean_part = min(data.mu);
sigma_mat = max(sqrt(data.Sigma));
sigma_part = min(sqrt(data.Sigma));

mat_prob=normpdf(Image_vec,mean_mat,sigma_mat);
part_prob=normpdf(Image_vec,mean_part,sigma_part);

% Dimensions of input image
sizeI = size(Image);

% Index is the correct indexing of a matrix
index = zeros(sizeI);
counter = 1;
for j = 1:sizeI(2)
    for i = 1:sizeI(1)
        index(i,j) = counter;
        counter = counter+1;
    end
end

%%
% What are the values for the
topedge = index(1,:);
leftedge = index(:,1);

% Length of image vector
len_Ivec = length(Image_vec);

ap_txt = fopen('adj_pairs_seg.txt','w');

% Find nearest neighbors for each pixel
for i = 1:len_Ivec
        if any(topedge==i)==0
            adj_pair = [i,i-1];
            adj_wt = abs(double(Image_vec(i))-double(Image_vec(i-1)));
            fprintf(ap_txt,'%d %d %d \n',adj_pair,adj_wt);
        end
        
        if any(leftedge==i)==0
            adj_pair = [i,i-sizeI(1)];
            adj_wt = abs(double(Image_vec(i))-double(Image_vec(i-sizeI(1))));
            fprintf(ap_txt,'%d %d %d \n',adj_pair,adj_wt);
        end
        
%         if any(bottomedge==i)==0
%             adj_pair = [i,i+1];
%             adj_wt = abs(double(Image_vec(i))-double(Image_vec(i+1)));
%             fprintf(ap_txt,'%d %d %d \n',adj_pair,adj_wt);
%         end
%         
%         if any(rightedge==i)==0
%             adj_pair = [i,i+sizeI(1)];
%             adj_wt = abs(double(Image_vec(i))-double(Image_vec(i+sizeI(1))));
%             fprintf(ap_txt,'%d %d %d \n',adj_pair,adj_wt);
%         end
end

fclose(ap_txt);

% Now read that giant bloody text file into matlab
adj_data = importdata('adj_pairs_seg.txt');
adj_pairs = adj_data(:,1:2);
adjpr_wts = adj_data(:,3);

% Delete the text file so we rewrite it every time
delete adj_pairs_seg.txt

clear adj_data i ap_txt topedge leftedge counter
%% Commence Graph Cutting
% close all

adj_wts = 1./adjpr_wts;
Infwts = adj_wts==Inf;
adj_wts(Infwts) = 2e0;
adjwt_scale = 5e0;

adj_wts = adj_wts.*adjwt_scale;

% Add source node
segment=digraph;
segment=addnode(segment,1);

% Add in nodes corresponding to matrix
segment=addnode(segment,len_Ivec);
endnode = size(segment.Nodes,1);

m2p_reg = 0;
m2p_scale = 1e1;
% m2p_wts = 1./((abs(Image_vec - mean_mat)+m2p_reg).*m2p_scale);
m2p_wts = (mat_prob+m2p_reg).*m2p_scale;
% segment = addedge(segment,(1:len_Ivec)+1,(1:len_Ivec)+endnode,m2p_wts);

% Add source weights
segment=addedge(segment,1,(1:len_Ivec)+1,m2p_wts);

% Add inplane weights for matrix
segment=addedge(segment,(adj_pairs(:,1)+1),(adj_pairs(:,2)+1),adj_wts);

% Add sinknode
segment = addnode(segment,1);
sinknode = size(segment.Nodes,1);

% Add the image to sink weights
p2t_reg = 0;
p2t_scale = 1e2;
% p2t_wts = (abs(Image_vec - mean_part)+10).*5.25;
% p2t_wts = 1./(abs(Image_vec - mean_part)+p2t_reg).*p2t_scale;
p2t_wts = (part_prob+p2t_reg).*p2t_scale;
segment = addedge(segment,(1:len_Ivec),sinknode,p2t_wts);

% Perform fist cut 
[mf,gf,cs,ct]=maxflow(segment,1,sinknode);
ct(end)=[];

% Copy the cut for convenience
ctcopy = ct-1;

% figure
% imshow(Image./255);

% figure
% imshow(imread('Rene88DT_GT/Segmented/slice-l-010.tiff'));

figure
seg_I = ones(length(Image_vec),1);
seg_I(ctcopy,1)=0;
seg_I = reshape(seg_I,[size(Image,1),size(Image,2)]);
imshow(seg_I)