function s = split_string(string, sep)
% function split_string(string, sep)
    
    if size(ver('Octave'), 1) == 1,
        s = strsplit(string, sep);
    else
        s = regexp(string, sep, 'split');
    end;
