function [still_running, line, start_time, status] = rproc_still_running(jobinfo);
% [still_running, line] = rproc_still_running(jobinfo);

status = 0;
still_running = 0 ;
global rproc_nqstat_output
global rproc_nqstat_time
global engine
start_time = [];
%id_len = 5;
if isequal(jobinfo.jobid, 0),
  % runs locally in background
  still_running = ~rproc_finished(jobinfo) ;
  line = sprintf('local job %s is still running: %s\n', jobinfo.prefix, jobinfo.log_fname) ;
  return ;
end ;

curtime = now;
if isempty(rproc_nqstat_time) || curtime-rproc_nqstat_time>0.5e-4
  %[ret, text]=unix(sprintf('qstat -u %s 2> /dev/null', whoami)) ;
  [ret, text]=unix(sprintf('/opt/torque/bin/qstat -u %s', whoami)) ;
  if ret==0 ;
    rproc_nqstat_output=text ;
    rproc_nqstat_time=curtime ;
  else
  	if ret==130
		fprintf('rproc_still_running interupted by user\n')
		status = -1;
		line = '';
		start_time = '';
	end
    warning('qstat failed') ;
  end ;
else
  text=rproc_nqstat_output ;
end ;

if ~any(text==sprintf('\n')),
  idx = [length(text)]';
else
  idx = [find(text==sprintf('\n')), length(text)]';
end ;
for i=1:length(idx)-1,
  line = text(idx(i):idx(i+1)-1) ;
  if length(line)> 0,
    if isequal(engine, 'octave'),
        items = strsplit(deblank(line), ' ') ;
    else
        items = regexp(deblank(line), ' ', 'split');
    end;
    for j=1:length(items)% assume that first non-empty item is the jobid
      if ~isempty(items{j}),
        break;
      end ;
    end
    %p=str2num(items{j}) ;
    if isequal(engine, 'octave'),
        p = strsplit(items{j}, '.') ;
    else
        p = regexp(items{j}, '\.', 'split');
    end;

    p = str2num(p{1});
    if p==jobinfo.jobid,
      still_running=1 ;
	  status = get_status(items);
	  still_running = check_status(status);

      if isempty(jobinfo.start_time)&&strcmp(status, 'R')
        start_time = now;
      else
        start_time = jobinfo.start_time;
	  end
      return ;
    end ;
  end ;
end ;
line = [] ;

return
function status = get_status(items)
  status = '';
  num = 0;
  for j=1:length(items)
    if ~isempty(items{j}),
      num = num+1;
    end
	if num == 10
	  status = items{j};
      break;
	end
  end
return

function ret = check_status(status)
  ret = 1;
  if strcmp(status, 'E') % Job is exiting after having run.
  	ret = 0;
  elseif strcmp(status, 'C') % Job is completed after having run.
  	ret = 0;
  elseif strcmp(status, 'S') % (Unicos only) Job is suspended.
  	ret = 0;
  end
return
