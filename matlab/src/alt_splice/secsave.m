function secsave(file, var, varname)

    if iscell(var),
        for i = 1:length(var),
            eval(sprintf('%s = var{i};', varname{i}));
        end;
        s = warning('error', 'MATLAB:save:sizeTooBigForMATFile');
        try
            save(file, varname{:});
        catch
            save(file, varname{:}, '-v7.3');
        end;
        warning(s);
    else
        eval(sprintf('%s = var;', varname));

        s = warning('error', 'MATLAB:save:sizeTooBigForMATFile');
        try
            save(file, varname);
        catch
            save(file, varname, '-v7.3');
        end;
        warning(s);
    end;
