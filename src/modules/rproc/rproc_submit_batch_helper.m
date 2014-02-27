function x = rproc_submit_batch_helper(parameters)
% rproc_submit_batch_helper(parameters)

fprintf('Executing a batch of %i jobs in a super-job\n', length(parameters)) 
pp=cwd ;
for i=1:length(parameters)
  cd(pp) ;
  fprintf('starting job %i in file %s\n', i, parameters(i).mat_fname) ;
  fprintf('=========================================\n') ;
  try
    start_proc(parameters(i).mat_fname, 0) ;
  catch
    disp('execution of start_proc failed') ;
  end ;
end ;

% remove files
for i=1:length(parameters)
  fname = parameters(i).mat_fname ;
  unix(['rm -rf ' fname]) ; % mat file
  unix(['rm -rf ' fname(1:end-2)]) ; % m file ;-)
end ;

% dummy return argument
x=0;
