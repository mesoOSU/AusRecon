function [T2R,flag]=calc_T2R(OR,CS_R,CS_T)

% Compute the variants: R2T Defines the orientation relationship as an 
% MTEX misorientation. This maps the crystal axes of austenite to the 
% crystal axes of martensite or from beta to alpha
if length(OR) == 3

    [OR,flag]=YardleyVariants(OR);
    %Note in MTEX convention is opposite Yardley convention USE Variants Transpose
    % for aus->Martensite using Variant'*aus_orientation->martensite
    
elseif length(OR) == 12
    flag = 0;
    OR = BurgersVariantsA2B;
end

Var1 = OR{1};
T2R=orientation('matrix',Var1,CS_R,CS_T);

end
