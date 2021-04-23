function TransOrs=generate_simulated_data(g_orient,OR,halfwidth,Num_TransOrs,CS_T)

% Define the crystal symmetry
CS_R = g_orient.CS;

% First define misorientation array based on the variant
% transformation. This accounts for the symmetry difference in the
% Titanium case, and still symmetrically accounts for the cubic crystal
% structure for austenite and martensite.
R2T = calc_R2T(OR,CS_R,CS_T);

%sample noisy austenite
Rec_noise=sample_halfwidth(g_orient,halfwidth,Num_TransOrs,CS_R,g_orient.SS);

%generate all martensite variants for each
potential_TransOrs=symmetrise(Rec_noise)*R2T;

keep=randperm(length(potential_TransOrs));
TransOrs=potential_TransOrs(keep(1:Num_TransOrs));

end

