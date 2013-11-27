function rproc_cleanup(jobinfo) ;
% rproc_cleanup(jobinfo)
  
for ix = 1:length(jobinfo)
  command = sprintf('rm -f %s %s %s %s %s\n', jobinfo(ix).mat_fname, jobinfo(ix).result_fname, ...
		    jobinfo(ix).m_fname, jobinfo(ix).log_fname, jobinfo(ix).qsublog_fname) ;
  unix(command) ;

  rproc_register('cleanup', jobinfo(ix)) ;

end

%unix('rm -f ' jobinfo.mat_fname ' ' jobinfo.result_fname ' ' jobinfo.m_fname ' ' jobinfo.log_fname) ;
