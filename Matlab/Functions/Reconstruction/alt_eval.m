function f = alt_eval(odf_component,orientations,varargin)
% Austin's port of mtex's eval function. Allows better control over NFSOFT
% commands. evaluates an odf at orientation(s) g
%
% Input
%  odf - @ODF
%  g   - @orientation
%
% Flags
%  even       - calculate even portion only
%
% Output
%  f   - values of the ODF at the orientations g
%
% See also
% kernel/sum_K kernel/K

persistent plans;
%NOTE TO FUTURE READERS: I don't think this needs to be added to the path
%every time, but also don't think it slows down the code, so I'm leaving it
%as is. feel free to change.
nfft_path = [mtex_path '\extern\nfft'];
addpath(nfft_path)
%('C:/Users/arger/workspace/AusRecon/Matlab/mtex-5.1.1/extern/nfft')
%addpath('mtex-5.1.1/extern/nfft')

% Command for if the user wants to clear all plans;
if check_option(varargin,'terminate_all_plans')
    for i = 1:(3+length(plans))
        try
            nfsoftmex('finalize',i-1);
        end
    end
    plans = strings(0);
else
    
    
    % Initial checks and values
    nfsoft_bandwidth = min(odf_component.bandwidth,65);
    ori_size = length(orientations);
    plan = 0;
    f = zeros(ori_size,1);
    if isempty(f), return; end
    
    % Check for existing plan that will work:
    if isempty(plans)
        plans = strings(0);
    else
        for i = 1:length(plans)
            try
                plan_parameters = split(plans(i),"_");
                plan_size = plan_parameters(1);
                plan_size = str2double(plan_size);
                plan_bandwidth = plan_parameters(2);
                plan_bandwidth = str2double(plan_bandwidth);
                if ori_size <= plan_size && plan_bandwidth == nfsoft_bandwidth
                    plan = i;
                    break
                end
            end
        end
    end
    
    % make a new plan if none already exists.
    if plan == 0
        % fprintf('new plan started\n')
        plan = initialize(orientations,odf_component)+1;
        plan_size =ori_size ;
        plan_bandwidth = nfsoft_bandwidth;
        plan_parameters = string(plan_size)+'_'+string(plan_bandwidth);
        plans(plan) = plan_parameters;
        % fprintf('plan %i made\n',plan)
    end
    
    
    if length(orientations) <plan_size
        orientations = repmat(orientations,ceil(plan_size/length(orientations)),1);
        orientations = orientations(1:plan_size,:,1);
    end
    
    %Run it
    nfft_plan = plan-1;
    Ldim = deg2dim(double(nfsoft_bandwidth+1));
    nfsoftmex('set_x',nfft_plan,Euler(orientations,'nfft').');
    nfsoftmex('precompute',nfft_plan);
    nfsoftmex('set_f_hat',nfft_plan,odf_component.components{1}.f_hat(1:Ldim));
    nfsoftmex('trafo',nfft_plan); % <-should probably figure out what this does
    nfsoft_results = real(nfsoftmex('get_f',nfft_plan));
    nfsoft_results = nfsoft_results(1:ori_size);
    
    % Final values
    f = f + reshape(nfsoft_results,ori_size,1);
end

end


function plan = initialize(orientations,odf_component)
% initializes
nfsoft_bandwidth = min(odf_component.bandwidth,65);
nfsoft_flags     = 2^4;
nfsoft_size      = 2*ceil(1.5*nfsoft_bandwidth);
ori_size = length(orientations);
plan = nfsoftmex('init',nfsoft_bandwidth,ori_size,nfsoft_flags,0,4,1000,nfsoft_size);
nfsoftmex('set_x',plan,Euler(orientations,'nfft').');
nfsoftmex('precompute',plan);
end

function plan = terminate(plan,plans)
% removes plan
nfsoftmex('finalize',plan);
plans(plan) = [];
end