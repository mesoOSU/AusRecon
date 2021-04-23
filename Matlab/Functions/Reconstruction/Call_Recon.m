function [myEBSD,Parent,Twin] = Call_Recon(myEBSD,IP,OP,numOrs)
% Function that calls either the mixed space or the single space (austenite
% and martensite) reconstruction functions
    rec_space = myEBSD.rec_space;       % Reconstruction space
    
    % If mixed space, call specific function
    if strcmp(rec_space,'Mixed')
        if isempty(numOrs)
            numOrs = 1e3;
        end
        [myEBSD,Parent,Twin] = Recon_Mixed_SymOP(myEBSD,numOrs,OP,IP);
    else
    % Call function that uses either martensite or austenite space
        [myEBSD] = Recon_SingSpace(myEBSD,IP,OP,numOrs);
        Twin = [];
        Parent = [];
    end
end
