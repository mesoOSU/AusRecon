function [Parent,PCount] = ParAssign(Parent,PInds,PCount,T_ebsd,Or)
    % Assign relevant data to Parent structure if grain was deemed to be
    % a PAG
    
    PCount = PCount+1;
    Parent.Indices{PCount,1} = T_ebsd(PInds).id;
    Parent.Or{PCount,1} = Or;
    Parent.Flg{PCount,1} = 1;
end

