function rproc_rerun(mess)
% rproc_rerun

if nargin==0, mess=''; end ;
  
global MATLAB_RETURN_VALUE

MATLAB_RETURN_VALUE=99 ;

global THIS_IS_A_RPROC_PROCESS  
  
if ~isempty(THIS_IS_A_RPROC_PROCESS) && THIS_IS_A_RPROC_PROCESS==1
  %error('RPROC:rerun', mess) ;
  %pause(10)
  exit
else
  warning('RPROC:rerun', mess) ;
end ;


