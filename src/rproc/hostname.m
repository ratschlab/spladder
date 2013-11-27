function s=hostname

fname = tempname ;

unix(['hostname > ' fname]) ;
fh=fopen(fname,'r') ;
if fh==-1
  warning(['Couldnt open ' fname]) ;
end ;
s=fscanf(fh, '%c') ;
fclose(fh) ;
unix(['rm -f ' fname]) ;
s=s(1:end-1) ;