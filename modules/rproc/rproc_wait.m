function [jobinfo, num_crashed] = rproc_wait(jobinfo, pausetime, frac_finished, resub_on, verbosity) 
% jobinfo=rproc_wait(jobinfo, pausetime,frac_finished,resub_on) 


% dummy change

global rproc_wait_jobinfo

rproc_wait_jobinfo=jobinfo ;

if nargin < 2
  pausetime = 120;
end
if nargin < 3
  frac_finished = 1 ;
end
if nargin < 4
  resub_on = 1 ;
end
if resub_on==1
  fprintf('\n\ncrashed jobs will be resubmitted by rproc_wait\n') ;
elseif resub_on==-1,
  fprintf('\n\ncrashed jobs may be resubmitted by rproc_wait\n') ;
else
  fprintf('\n\ncrashed jobs will not be resubmitted by rproc_wait\n') ;
end
if nargin<5,
  verbosity=2 ;
end ;


num_jobs = 0 ;
num_crashed = 0 ;
for i=1:length(jobinfo),
  if jobinfo(i).created==1 ;
    if isempty(jobinfo(i).time)
      warning('job created but not submitted yet. ignoring') ;
      jobinfo(i).created=0 ;
    else
      num_jobs=num_jobs+1 ;
    end ;
  end ;
end ;

num_finished=0 ; first_iter = 1 ;
while (num_finished<num_jobs*frac_finished) || (num_crashed>0)
  if ~first_iter,
    pause(pausetime) ;
  end ;
  first_iter = 0 ;
  num_finished=0 ;
  num_crashed=0 ;
  crashed_files='log files of crashed jobs:' ;
  for id=1:length(jobinfo)
    cur_finished = rproc_finished(jobinfo(id)) ;
    [still_running, qstat_line, start_time, status]=rproc_still_running(jobinfo(id));
	if status == -1
		return
	end
    jobinfo(id).start_time = start_time;
    if cur_finished
      num_finished = num_finished +1 ;
    elseif ~still_running 
      num_finished = num_finished + 1 ;
      num_crashed = num_crashed + 1 ;
      crashed_files = sprintf('%s\n%s', crashed_files, jobinfo(id).log_fname) ;
      if isempty(jobinfo(id).crashed_time)
        jobinfo(id).crashed_time = now ;
      elseif 24*60*(now-jobinfo(id).crashed_time) > max(3*pausetime/60, 0.1)  && (resub_on==1 || (resub_on==-1 && jobinfo(id).resubmit>=jobinfo(id).retries+1)),
        if resub_on==1,
          [reachedlimit, jobwalltime] = rproc_reached_timelimit(jobinfo(id)) ;
          if reachedlimit, % check whether the job has been killed because it reached the time limit
            if verbosity>=1,
              fprintf('job has been canceled because it used %1.0fs, but time limit was %1.0fs walltime.\nhence, we increase the time limit to %1.0fs.\n', ...
                      jobwalltime, jobinfo(id).time*60, max(jobinfo(id).time, jobwalltime)*2);
            end ;
            jobinfo(id).time = max(jobinfo(id).time, jobwalltime/60)*2 ;
          end ;
        elseif resub_on==-1,
          jobinfo(id).time = jobinfo(id).time_req_resubmit(jobinfo(id).retries+1) ;
          jobinfo(id).Mem = jobinfo(id).mem_req_resubmit(jobinfo(id).retries+1) ;
          jobinfo(id).start_time = [];
          if verbosity>=1,
			fprintf('resubmitting job (%i) with new time and memory limitations: %iMb and %i minutes (retry #%i)\n', jobinfo(id).jobid, jobinfo(id).Mem, jobinfo(id).time, jobinfo(id).retries+1) ;
          end ;
        end ;
        if verbosity>=2,
          fprintf('log file of previous attempt %s\n', jobinfo(id).log_fname) ;
        end ;
        jobinfo(id) = rproc_resubmit(jobinfo(id)) ;
        jobinfo(id).crashed_time = [] ;
        num_finished = num_finished - 1 ;
      end ;
    else
      if verbosity>=2,
        fprintf('%s', qstat_line) ;
      end ;
    end ;
	%% hard_time_limit in minutes
	if ~isempty(jobinfo(id).start_time) && 24*60*(now-jobinfo(id).start_time) > jobinfo(id).hard_time_limit
      fprintf('delete job (%i) because hard time limit (%imin) was reached\n', jobinfo(id).jobid, jobinfo(id).hard_time_limit);
      unix(sprintf('qdel %i', jobinfo(id).jobid));
	end
  end ;
  if verbosity>=1,
    fprintf('\n%i of %i jobs finished (%i of them crashed) \n', ...
            num_finished, num_jobs, num_crashed) ;
  end ;
  if verbosity>=2,
    if sum(crashed_files==sprintf('\n'))>0,
      fprintf('%s\n', crashed_files) ;
    end ;
  end ;
  if resub_on==0 && num_finished==num_jobs*frac_finished
    break
  end
  if resub_on==-1 && num_finished==num_jobs*frac_finished
    all_tried = 1 ;
    for i=1:length(jobinfo)
      fin = rproc_finished(jobinfo(i));
      if jobinfo(i).resubmit>=jobinfo(i).retries+1 && ~fin,
        all_tried = 0 ;
      end ;
    end ;
    if all_tried,
      break
    end ;
  end
end ;

pause(1);

