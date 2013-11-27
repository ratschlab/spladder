function [jobinfo, meta_jobinfo] = rproc_submit_batch(jobinfo, blocksize) 
% [jobinfo, meta_jobinfo] = rproc_submit_many(jobinfo, blocksize) 

%blocksize = 1 

meta_jobinfo=rproc_empty(0) ;  

time_per_submission = 1/60 ; % 1 seconds

meta_i = 1 ;
time_per_metajob=zeros(ceil(length(jobinfo)/blocksize),1) ;
metablockassignment=zeros(1,length(jobinfo));
[tmp,sidx]=sort(-[jobinfo.time]) ;
for i=sidx ;
  [tmp,idx] = min(time_per_metajob + linspace(-time_per_submission*length(time_per_metajob), 0, length(time_per_metajob))') ;
  metablockassignment(i)=idx ;
  time_per_metajob(idx) = time_per_metajob(idx) + jobinfo(i).time ;
end ;

meta_i = 1 ;
for i=1:ceil(length(jobinfo)/blocksize)
  %idx=i:min(i+blocksize-1, length(jobinfo)) ;
  idx = find(metablockassignment==i) ;
  if isempty(idx),
    continue ;
  end ;
  for j=1:length(idx),
    options = jobinfo(idx(j)).options ;
    options.submit_now = 0 ;
    jobinfo(idx(j))=rproc(jobinfo(idx(j)).ProcName, jobinfo(idx(j)).P1, jobinfo(idx(j)).Mem, ...
                          options, jobinfo(idx(j)).time) ;
  end ;
  jobinfo_ = jobinfo(idx) ;
  options = jobinfo(idx(1)).options ;
  options.submit_now = 1 ;
  options.verbosity = 1 ;
  %options.immediately = 0 ;
  memory_MB=max([jobinfo_.Mem]) ;
  minutes=sum([jobinfo_.time]) ;
  fprintf('submitting job %i/%i (%i subjobs) \r', i, ceil(length(jobinfo)/blocksize), length(idx)) ;
  meta_jobinfo(meta_i) = rproc('rproc_submit_batch_helper', jobinfo_, memory_MB, options, minutes) ;
  
  for j=1:length(idx),
    jobinfo(idx(j)).log_fname = meta_jobinfo(meta_i).log_fname ;
    jobinfo(idx(j)).jobid = meta_jobinfo(meta_i).jobid ;    
    jobinfo(idx(j)).submission_time = meta_jobinfo(meta_i).submission_time ;    
  end ;
  meta_i=meta_i+1 ;
end ;
fprintf('\n') ;

return

% a test
jobinfo = rproc_create('ls', 'a*', 200, [], 1) ;
jobinfo(2) = rproc_create('ls', 'b*', 200, [], 1) ;

[jobinfo_, meta_jobinfo]=rproc_submit_batch(jobinfo,2) ;
