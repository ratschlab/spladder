function jobinfo=rproc_empty(N) ;
%  jobinfo=rproc_empty ;
  
if nargin<1, N=1; end ;

NN=N ;
if N==0, N=1 ; end ;

for i=1:N,
  jobinfo(i).ProcName = [] ;
  jobinfo(i).P1 = [] ;
  jobinfo(i).Mem = [] ;
  jobinfo(i).options = [] ;
  jobinfo(i).time = [] ;
  jobinfo(i).prefix = [] ;
  jobinfo(i).mat_fname = [] ;
  jobinfo(i).result_fname = [] ;
  jobinfo(i).m_fname = [] ;
  jobinfo(i).log_fname = [] ;
  jobinfo(i).qsublog_fname = [] ;
  jobinfo(i).jobid = -1 ;  
  jobinfo(i).submission_time = [] ;
  jobinfo(i).retries = 0 ;
  jobinfo(i).created = 0 ;
  jobinfo(i).time_of_loss = [] ;
  jobinfo(i).crashed_time=[] ;
  jobinfo(i).maxvmem = '';
  jobinfo(i).resubmit = 0 ;
  jobinfo(i).time_req_resubmit = [] ;
  jobinfo(i).mem_req_resubmit = [] ;
  jobinfo(i).data_size = [] ;
  jobinfo(i).hard_time_limit = inf ;
  jobinfo(i).start_time = [] ;
end 

if NN==0, jobinfo(1)=[] ; end ;
