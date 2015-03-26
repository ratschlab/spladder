function [retval1, retval2] = rproc_result(jobinfo, read_attempts)
% [retval1, retval2] = rproc_result(jobinfo, [read_attempts])
  
if ~fexist(jobinfo.result_fname),
  att = 1;
  while ~fexist(jobinfo.result_fname),
    disp('Job not finished yet. Waiting for result file to appear.');
    if exist('fflush', 'builtin'),
      fflush(1);
    end
    if nargin > 1 && att > read_attempts,
      error('Unable to load result from %s', jobinfo.result_fname);
    end
    pause(10);
    att = att + 1;
  end
end

L = load(jobinfo.result_fname);
retval1 = L.retval1;
retval2 = L.retval2;
