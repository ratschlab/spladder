function rproc_clean_register()

fname = '~/tmp/rproc.log' ;

global engine

jobids = [] ;
parent_jobids = [] ;

fd=fopen(fname, 'r') ;
while ~feof(fd),
  line = fgetl(fd) ;
  if ~ischar(line), break ; end ;
  if isempty(line), continue ; end ;


  if isequal(engine, 'octave'),
      items = strsplit(line, '\t') ;
  else
      items = regexp(line, '\t', 'split');
  end;
  if length(items)<4, continue ; end ;
  if isempty(items{1}), continue ; end ;
  if isempty(items{4}), continue ; end ;
  if items{1}(1)>='0' && items{1}(1)<='9' && ((items{4}(1)>='0' && items{4}(1)<='9') || items{4}(1)=='-'),
    jobids(end+1) = str2num(items{1}) ;
    parent_jobids(end+1) = str2num(items{4}) ;
  end ;
end ;

[ret,text]=system('/opt/torque/bin/qstat') ;
idx=find(text==sprintf('\n')) ;

running_jobids = [] ; 
for i=1:length(idx)-1,
  line = text(idx(i)+1:idx(i+1)-1) ;
  if isequal(engine, 'octave'),
      items = strsplit(line, ' ') ;
  else
      items = regexp(line, ' ', 'split');
  end;
  if isempty(items{1}),
    items(1)=[] ;
  end ;
  if items{1}(1)>='0' && items{1}(1)<='9',
    running_jobids(end+1) = str2num(items{1}) ;
  end ;
end ;

[tmp,idx1,idx2]=intersect(jobids, running_jobids) ;

for i=1:length(idx1)
  if ~any(running_jobids == parent_jobids(idx1(i))) && parent_jobids(idx1(i))~=-1 ;
    fprintf('job %i is still running, but the parent job %i not\n', jobids(idx1(i)), parent_jobids(idx1(i)))
  %else
  %  fprintf('job %i and parent job %i are still running\n', jobids(idx1(i)), parent_jobids(idx1(i)))
  end ;
end ;

%keyboard
