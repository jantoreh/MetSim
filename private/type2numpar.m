% Define type2numpar function
function out = type2numpar(type)
    switch type
        case {'n','ln','w2','vm'}
            out = 2;
        case {'w3'}
                out = 3;
        case {'n2','vm2'}
            out=6;
        case {'vm3'}
            out=9;
    end
end
