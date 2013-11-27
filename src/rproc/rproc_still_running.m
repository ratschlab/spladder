function [still_running, line, start_time, status] = rproc_still_running(jobinfo);
% [still_running, line] = rproc_still_running(jobinfo);

status = 0;
still_running = 0 ;
global rproc_nqstat_output
global rproc_nqstat_time
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
  [ret, text]=unix(sprintf('qstat -u %s', whoami)) ;
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
  %elems = separate(line(2:end),' ',1) ;
  %if length(line)>=9,
  %  if str2num(elems{1})==jobinfo.jobid,
  if length(line)> 0,
    items=separate(deblank(line),' ') ;
    for j=1:length(items)% assume that first non-empty item is the jobid
      if ~isempty(items{j}),
        break;
      end ;
    end
    %p = textscan(line,'%d%s',1);
    p=str2num(items{j}) ;
    %if str2num(line(1:id_len))==jobinfo.jobid,
    if p==jobinfo.jobid,
      still_running=1 ;
      %try
      %  %[ret, mem]=unix(sprintf('qstat -j %i 2> /dev/null | grep maxvmem | cut -d "," -f 5 | cut -d "=" -f 2', p));
      %  [ret, mem]=unix(sprintf('qstat -j %i | grep maxvmem | cut -d "," -f 5 | cut -d "=" -f 2', p));
      %catch
        mem = '';
      %end
	  status = get_status(items);
	  still_running = check_status(status);

      if isempty(jobinfo.start_time)&&strcmp(status, 'r')
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
	if num == 5
	  status = items{j};
      break;
	end
  end
return

function ret = check_status(status)
  ret = 1;
  if strcmp(status, 't')
  	ret = 0;
  elseif strcmp(status, 'Eqw')
  	ret = 0;
  elseif strcmp(status, 'dt')
  	ret = 0;
  elseif strcmp(status, 'dr')
  	ret = 0;
  end
return
