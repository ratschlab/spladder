function x = fexist(filename,gz)
% FEXIST Check if file exists
%        FEXIST('filename') returns:
%        0 if filename does not exist
%        1 if filename exists in the given path

if nargin<2,
  gz=0;
end ;

fd = fopen(filename);
if fd ~= -1
  x = 1;
  fclose(fd);
elseif gz,
  fd = fopen([filename '.gz']);
  if fd~=-1,
    x=1 ;
    fclose(fd) ;
  else
    x=0 ;
  end ;
else
  x = 0;
end
x=logical(x) ;