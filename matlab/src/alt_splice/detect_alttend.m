function [idx_alt_tend, exon_alt_tend] = detect_alttend(genes) ;

idx_alt_tend = {};
exon_alt_tend = {};

for ix=1:length(genes)
  if (mod(ix,1000)==0)
    fprintf(1,'.');
  end
  num_exons = size(genes(ix).splicegraph{1},2) ;
  vertices = genes(ix).splicegraph{1} ;
  edges = triu(genes(ix).splicegraph{2}) ;
  
  term_alt_idx = [] ;
  for i=1:num_exons,
    if sum(edges(i,:),2)==0,
      term_alt_idx(end+1) = i ;
    end ;
  end ;
  if length(term_alt_idx)>1
    idx_alt_tend{end+1} = ix;
    exon_alt_tend{end+1}.exons = term_alt_idx;
  end
end

fprintf(1,'\nTotal alternative transcription ends:\t\t\t\t%d\n',...
        length(idx_alt_tend));
