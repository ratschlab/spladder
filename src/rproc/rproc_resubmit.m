function jobinfo2 = rproc_resubmit(jobinfo, force);
% jobinfo2 = rproc_resubmit(jobinfo);

if nargin<2, force=1 ; end ;

if length(jobinfo)==0
  jobinfo2 = jobinfo;
  return
elseif length(jobinfo)>1,
  jobinfo2 = jobinfo ;
  for i=1:length(jobinfo)
    jobinfo2(i) = rproc_resubmit(jobinfo(i)) ;
  end 
  return ;
end ;
  
if (jobinfo.retries>=0) & (rproc_time_since_submission(jobinfo)<1),
  %warning('job was submitted less than a minute ago. not resubmitted.') ;
  jobinfo2=jobinfo ;
  return ;
end ;

if jobinfo.retries>=3 & rand(1)<0.5,
  if jobinfo.options.verbosity>=0,
    warning(sprintf('job has already been submitted %i times', ...
                    jobinfo.retries)) ;
  end ;
  if jobinfo.options.verbosity>0,
    fprintf('check file %s\n', jobinfo.log_fname) ;
  end ;
  jobinfo2=jobinfo ;
  return ;
end ;

if (jobinfo.retries>=0)
  still_running = rproc_still_running(jobinfo);
  if still_running, 
    if jobinfo.options.verbosity>0,
      fprintf('.') ;
    end ;
    jobinfo2=jobinfo ;
    jobinfo2.time_of_loss = [] ;
    return ;
  end ;
  if isempty(jobinfo.time_of_loss)
    jobinfo.time_of_loss=now ;
  end ;
  % more than a minute lost?
  if ~force & (24*60*(jobinfo.time_of_loss-now)<1)
    %warning('do not resubmit yet ... ') ;
    jobinfo2=jobinfo ;
    return ;
  end ;
  %rproc_cleanup(jobinfo) ;
end ;

%fprintf('\nresubmitting job\n') ;
jobinfo2 = rproc(jobinfo.ProcName, jobinfo.P1, jobinfo.Mem, ...
                 jobinfo.options, jobinfo.time) ;

if jobinfo.jobid~=-1,
   % increase only, if it has not been resubmitted before
  jobinfo2.retries = jobinfo.retries + 1 ;
end ;
