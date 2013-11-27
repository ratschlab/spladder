function [jobinfo]=rproc_create(ProcName, P1, Mem, options, time)
% [jobinfo]=rproc(ProcName, P1, Mem, options, time)
%
% time in minutes
% mem in mb
  
if nargin<3, Mem=100 ; end ;
if nargin<4, options=[] ; end ;
if nargin<5, time=100*24*60 ; end ;

jobinfo=rproc_empty ;

jobinfo.ProcName = ProcName ;
jobinfo.P1 = P1 ;
jobinfo.Mem = Mem ;
jobinfo.options = options ;
jobinfo.time = time ;

jobinfo.created=1 ;
jobinfo.retries=-1 ;
