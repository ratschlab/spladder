function myassert(Val, Text)
%
% checks a condition and aborts if it has been failed (Val==0)


global debug_assert ;

%dbstack ;
if isempty(debug_assert),
  debug_assert=1 ;
elseif debug_assert==0,
  return ;
end ;

if ~Val,
	if nargin<2,
		rep=sprintf('assertion failed (%i)', Val) ;
	else
		rep=Text ;
	end ;
	dbstack ;
	error(rep) 
end ;