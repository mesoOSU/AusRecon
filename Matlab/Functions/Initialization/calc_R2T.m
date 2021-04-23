function [R2T,flag]=calc_R2T(OR,CS_R,CS_T)

% Compute the variants: R2T Defines the orientation relationship as an 
% MTEX misorientation. This maps the crystal axes of austenite to the 
% crystal axes of martensite or from beta to alpha
if length(OR) == 3

    [OR,flag]=YardleyVariants(OR);
    %Note in MTEX convention is opposite Yardley convention USE Variants Transpose
    % for aus->Martensite using Variant'*aus_orientation->martensite
    Var1 = (OR{1})';
    
elseif length(OR) == 12
    flag = 0;
    Var1 = OR{1};
    
end

R2T=orientation('matrix',(Var1),CS_T,CS_R);
%M2A=orientation('matrix',OR{1},CS_A,CS_M);

end

