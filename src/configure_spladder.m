function CFG = configure_spladder(bam_files, anno_file, out_dir, log_file, user_settings, confidence)

    %%% set custom confidence level
    if nargin <= 5,
        confidence = 3;
    end;
    if confidence < 0 || confidence > 3,
        error('ERROR: Pre-defined confidence levels must be within the range of 0 to 3\n');
    end;

    %%% load default settings
    default_settings

    %%% add user settings
    if nargin > 4,
        if exist(user_settings, 'file'),
            if length(user_settings) > 2 && ~strcmp(user_settings(end-1:end), '.m'),
                eval(user_settings(end-1:end));
            else
                eval(user_settings);
            end;
        else
            error(sprintf('ERROR: Can not find file %s to load user settings!\n', user_settings));
        end;
    end;

    if nargin > 3,
        CFG.fd_log = fopen(log_file, 'w');
    else
        CFG.fd_log = 1;
    end;

    if nargin > 2,
        CFG.out_dir = out_dir;
    else
        CFG.out_dir = pwd();
    end;

    if nargin > 1,
        CFG.anno = CFG.anno_file;
        CFG.bam_files = separate(bam_files, ',');
    else
        error('ERROR: BAM file(s) and annotation need to be provided!\n');
    end;
