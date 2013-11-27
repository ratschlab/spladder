function time = rproc_time_since_submission(jobinfo)
% time = rproc_time_since_submission(jobinfo)
% returns time in minutes since submission  
  
time = 24*60*(now - jobinfo.submission_time) ; 
  