function [result, jobwalltime] = rproc_reached_timelimit(jobinfo)
% [result, jobwalltime] = rproc_reached_timelimit(jobinfo)

  result = false ;

  %[s,w]=system(sprintf('qacct -j %i | grep ru_wallclock|sed ''s/ru_wallclock//g''', jobinfo.jobid));
  [s,w]=system(sprintf('qstat -ta | grep %i | sed -e ''s/  */ /g'' | cut -f 11 -d '' ''', jobinfo.jobid));
  if ~isempty(strfind(w, 'error')),
    jobwalltime = -1 ;
    return ;
  end ;
  ii=find(w==10) ; % extract last line .. ignore the error lines
  ii=[1 ii] ;
  try
      jobwalltime=str2num(w(ii(end-1):ii(end))) ; % in seconds
  catch
    jobwalltime = -1 ;
  end;
  if isempty(jobwalltime), 
    jobwalltime = -1 ; 
    return ; 
  end ;
  if ~(jobwalltime>0 & jobwalltime<10000*3600) ; % sanity checks
      warning('invalid output from qacct\n')
      jobwalltime = -1 ;
      return ;
  end 
  if jobwalltime > jobinfo.time*60
    result = true ;
  end ;
