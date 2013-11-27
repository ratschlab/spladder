function isfinished = rproc_finished(jobinfo) ;
% isfinished = rproc_finished(jobinfo) ;

isfinished = 0 ;
if jobinfo.jobid==-1, return ; end ;

if fexist(jobinfo.result_fname),
  isfinished = 1 ;
  rproc_register('finished', jobinfo) ;
end ;
