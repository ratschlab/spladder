disp('rproc finishing') ;

global MATLAB_RETURN_VALUE

if ~isempty(MATLAB_RETURN_VALUE)
  fprintf('exit code %i\n', MATLAB_RETURN_VALUE) ;
  if ~isempty(getenv('MATLAB_RETURN_FILE'))
    fd=fopen(getenv('MATLAB_RETURN_FILE'),'w+') ;
    if fd~=-1,
      fprintf(fd, '%i\n', MATLAB_RETURN_VALUE(1)) ;
      fclose(fd) ;
    end ;
  else
    warning('environment MATLAB_RETURN_FILE not defined') ;
  end ;
end ;

%disp('rproc finished') ;
