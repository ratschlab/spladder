function rproc_cmd(unix_cmd, jobinfo)
% rproc_cmd(unix_cmd, jobinfo)

for i=1:length(jobinfo(:)),
  if ~isempty(jobinfo(i).jobid) && jobinfo(i).jobid~=-1,
    unix(sprintf('%s %i', unix_cmd, jobinfo(i).jobid)) ;
    %assert(s==0) ;
  end ;
end ;