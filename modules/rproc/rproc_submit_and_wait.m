function [jobinfo,frac_finished]=rproc_submit_and_wait(jobinfo, finish_frac, jobtimeout) 
% frac_finished=rproc_submit_and_wait(jobinfo, finish_frac, jobtimeout) 
  
num_jobs = 0 ;
for i=1:length(jobinfo),
  if jobinfo(i).created==1 ;
    num_jobs=num_jobs+1 ;
  end ;
end ;

%j=0 ;
%for i=1:length(jobinfo),
%  if jobinfo(i).created==1 ;
%    j=j+1 ;
%    fprintf('submitting %i jobs: %i\r', num_jobs, j) 
%    jobinfo(i)=rproc_resubmit(jobinfo(i)) ;
%  end ;
%end ;
%fprintf('\n')

num_finished=0 ;
while(num_finished/num_jobs<finish_frac)
  num_finished=0 ;
  for id=1:length(jobinfo)
    if rproc_finished(jobinfo(id)) 
      num_finished = num_finished +1 ;
    else
      if jobinfo(id).created==1
        if rproc_time_since_submission(jobinfo(id))>jobtimeout,
          warning('job took longer than timeout. Killing and restarting it');
          rproc_kill(jobinfo(id)) ;
        end ;
        jobinfo(id) = rproc_resubmit(jobinfo(id), 1) ;
      end ;
    end ;
  end ;
  fprintf('waiting for jobs to finish: %i/%i  \r', num_finished, num_jobs) ;
  if (num_finished/num_jobs<finish_frac), pause(10) ; end ;
end ;

fprintf('\n') ;
