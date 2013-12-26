function split_string(string, sep)
% function split_string(string, sep)
    
    if size(ver('Octave'), 1),
        return strsplit(string, sep);
    else
        return regexp(string, sep, 'split');
    end;
