% Append a structure without overwriting
% Michael Nickerson 2022
function S1 = appendstruct(S1, S2)
    if ~isstruct(S1) || ~isstruct(S2)
        return;
    end
    
    for ff = fieldnames(S2)'
        ff = ff{1};
        if isfield(S1, ff) && isstruct(S1.(ff)) && isstruct(S2.(ff))
            % If substruct, append
            S1.(ff) = appendstruct(S1.(ff), S2.(ff));
        elseif isfield(S1, ff)
            % If different, rename
            if ~isequaln(S1.(ff), S2.(ff))
                S1.(ff+"_") = S2.(ff);
            end
        else
            % Does not already exist; just add
            S1.(ff) = S2.(ff);
        end
    end
end
