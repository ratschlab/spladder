function [idx_alt_tstart, exon_alt_tstart] = detect_alttstart(genes) ;

idx_alt_tstart = {};
exon_alt_tstart = {};

for ix=1:length(genes)
  if (mod(ix,1000)==0)
    fprintf(1,'.');
  end
  num_exons = size(genes(ix).splicegraph{1},2) ;
  vertices = genes(ix).splicegraph{1} ;
  edges = triu(genes(ix).splicegraph{2}) ;
  
  init_alt_idx = [] ;
  for i=1:num_exons,
    if sum(edges(:,i),1)==0, % initial
      init_alt_idx(end+1) = i ;
    end ;
  end ;

  if length(init_alt_idx)>1,
    idx_alt_tstart{end+1} = ix;
    exon_alt_tstart{end+1}.exons = init_alt_idx;
  end
end

fprintf(1,'\nTotal alternative transcription starts:\t\t\t\t%d\n',...
        length(idx_alt_tstart));
