function rproc_kill(jobinfo) ;
% rproc_kill(jobinfo)

if isequal(jobinfo,'wait'),
  global rproc_wait_jobinfo
  jobinfo = rproc_wait_jobinfo ;
end ;

for i=1:length(jobinfo)  
  if ~isempty(jobinfo(i).jobid) & jobinfo(i).jobid>0
    unix(sprintf('qdel %i 2> /dev/null', jobinfo(i).jobid)) ;

    rproc_register('kill', jobinfo(i)) ;
  end ;
end ;
