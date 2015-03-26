function [sum_alt, sum_const, idx_alt, idx_con, genes] = detect_altsplicing(genes) ;
% function [sum_alt, sum_const, idx_alt, idx_con, genes] = detect_altsplicing(genes) ;

level_e = zeros(1,length(genes));
level_i = zeros(1,length(genes));
for ix = 1:length(genes)
  if (mod(ix,50)==0)
    fprintf(1,'.');
  end
  [level_e(ix), level_i(ix)] = detectsplicegraph(genes(ix));
  if level_e(ix)>1 | level_i(ix)>1,
    genes(ix).is_alt=1 ;
  else
    genes(ix).is_alt=0 ;
  end ;
end

sum_alt   = sum(or(level_i>1,level_e>1)) ;
sum_const = sum(and(level_i<=1,level_e<=1)) ;

fprintf(1,'\n\nTotal alternatively spliced:\t\t\t\t\t%d\n',...
	sum_alt);
fprintf(1,'Total constitutively spliced:\t\t\t\t\t%d\n',...
	sum_const);

[dum,idx_alt] = find(or(level_i>1,level_e>1));
[dum,idx_con] = find(and(level_i<=1,level_e<=1));

