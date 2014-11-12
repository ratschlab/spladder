function user_name = whoami
% user_name = whoami
  
global whoami_user_name  
  
if isempty(whoami_user_name)
  
  [ret, whoami_user_name]=unix('whoami') ; 
  if ret~=0,
    if keyboard_allowed(),
      keyboard ;
    end ;
    whoami_user_name='unknown' ;
  end ;
end ;

user_name= whoami_user_name ;
idx=find(user_name<30, 1, 'first') ;
user_name = user_name(1:idx-1) ;
