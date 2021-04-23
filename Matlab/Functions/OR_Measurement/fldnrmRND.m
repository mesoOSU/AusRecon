function x=fldnrmRND(mu,sigma,varargin)

% parse varargins
if ~exist('varargin','var')
    error('Too few input parameters.')
else
    str=[];
    for i=1:length(varargin)
        if ischar(varargin{i})
            str=horzcat(str,[varargin{i} ',']);
        elseif isnumeric(varargin{i})
            str=horzcat(str,[num2str(varargin{i}) ',']);
        else
            error('Verify that varargins are valid inputs for randn')
        end
    end
    str=['x=randn(' str(1:end-1) ');'];
end

eval(str); % generate normally distributed values
x=mu+sigma.*x; % apply location and scale
x=abs(x); % fold the distribution