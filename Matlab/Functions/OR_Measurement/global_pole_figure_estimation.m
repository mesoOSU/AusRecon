function [ austenite_proposal ] = global_pole_figure_estimation(austenite,CS,SS,smooth)

austenite=project2FundamentalRegion(austenite,orientation(idquaternion,CS,SS));
Euler_angles=niceEuler(austenite,'Bunge')/degree;
bin_edges=(ceil(Euler_angles));

bin_edges(find(bin_edges==0))=1;
accumulator=accumarray(bin_edges,1,[360,90,90]);

if smooth
    accumulator=smooth3(accumulator,'box',3);
end

possible_austenite=find(accumulator==max(accumulator(:)));

[approx_phi1,approx_Phi,approx_phi2]=ind2sub([360,90,90],possible_austenite);


 test_orientations=orientation(euler2quat(approx_phi1*degree,approx_Phi*degree,approx_phi2*degree,'Bunge'),CS,SS);
for ii=1:length(test_orientations) 
 deviation=angle(austenite,test_orientations(ii));
 deviation(deviation>pi)=0;
 weights=cos(deviation);
 weights(deviation>5*degree)=0;
 try
 austenite_proposal(ii)=mean(austenite,'weights',weights);
 catch
 keyboard
     
 end
 %Error(ii)=angle(austenite_proposal(ii),austenite_initial)/degree;
end
end

