function rproc_register(action, jobinfo) ;
% rproc_register(action, jobinfo) ;

this_jobid = getenv('JOB_ID') ;
if ~isempty(this_jobid),
  this_jobid = str2num(this_jobid) ;
else
  this_jobid = -1 ;
end ;

[ret, timedate] = system('date') ;
assert(ret==0) ;
timedate(timedate<20) = [] ;

rproc_log_fname = '~/tmp/rproc.log' ;

if ~fexist(rproc_log_fname),
  fd = fopen(rproc_log_fname, 'a+') ;
  fprintf(fd, '# prefix\taction\tparent jobid\tjobid\tfunction\ttime\n') ;
  fclose(fd) ;
end ;

fd = fopen(rproc_log_fname, 'a+') ;
if fd>=0,
  fprintf(fd, '%i\t%s\t%s\t%i\t%s\t%s\n', jobinfo.jobid, jobinfo.prefix, action, this_jobid, jobinfo.ProcName, timedate) ;
  fclose(fd) ;
end ;
