function [idx_exon_skips, exon_exon_skips] = detect_exonskips(genes, idx_alt) ;
% [idx_exon_skips, exon_exon_skips] = detect_exonskips(genes, idx_alt) ;

idx_exon_skips = [];
exon_exon_skips = [];
for ix=idx_alt
  if (mod(ix,50)==0)
    fprintf(1,'.');
  end
  num_exons = size(genes(ix).splicegraph{1},2);
  edges = genes(ix).splicegraph{2};
  for exon_idx = 1:num_exons-2
    for exon_idx1 = (exon_idx+1):(num_exons-1)
      for exon_idx2 = (exon_idx1+1):num_exons
        if (edges(exon_idx,exon_idx1))&&...
              (edges(exon_idx,exon_idx2))&&...
              (edges(exon_idx1,exon_idx2))
          idx_exon_skips = [idx_exon_skips,ix];
          exon_exon_skips = [exon_exon_skips,[exon_idx; exon_idx1; exon_idx2]];
        end
      end
    end
  end
end
sum_skips=length(idx_exon_skips) ;
fprintf(1,'\n\nNumber of single exon skips:\t\t\t\t\t%d\n', sum_skips);


