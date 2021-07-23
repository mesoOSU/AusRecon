function Plot_Blocks(myEBSD,varargin)
%Plots Blocks. Plotting is done in Plot_Microstructure_Hierarchy
% mbar      : True/False, turns Micronbar on/off
% key       : True/False, turns color key on/off
% A_outline : True/False, turns Austenite border drawing on/off
% P_outline : True/False, turns Packet border drawing on/off
% B_outline : True/False, turns Block border drawing on/off
% V_outline : True/False, turns Variant border drawing on/off

% set all options to true by default, then update from varargin
p = inputParser;
addRequired(p,'myEBSD')
addOptional(p,'mbar',1)
addOptional(p,'cbar',1)
addOptional(p,'A_outline',1)
addOptional(p,'P_outline',1)
addOptional(p,'B_outline',0)
addOptional(p,'V_outline',0)
parse(p,myEBSD,varargin{:})

Plot_Microstructure_Hierarchy(myEBSD,'V',0,'B',1,'P',0,...
    'mbar',p.Results.mbar,...
    'cbar',p.Results.cbar,...
    'A_outline',p.Results.A_outline,...
    'P_outline',p.Results.P_outline,...
    'B_outline',p.Results.B_outline,...
    'V_outline',p.Results.V_outline);
end