function [idx_exon_intron5, exon_exon_intron5, idx_exon_intron3, exon_exon_intron3] = detect_exonintron(genes, idx_alt) ;

% Exon in intronic region, may be intron retention, but not complete.
idx_exon_intron5 = [];
exon_exon_intron5 = [];
idx_exon_intron3 = [];
exon_exon_intron3 = [];
for ix=idx_alt
  if (mod(ix,50)==0)
    fprintf(1,'.');
  end
  num_exons = size(genes(ix).splicegraph{1},2);
  vertices = genes(ix).splicegraph{1};
  edges = genes(ix).splicegraph{2};
  for exon_idx1 = 1:num_exons-1
    for exon_idx2 = (exon_idx1+1):num_exons
      if edges(exon_idx1, exon_idx2),
        for exon_idx3=1:num_exons
          %keyboard ;
          if (vertices(2,exon_idx3) > vertices(2,exon_idx1)) & ...
             (vertices(2,exon_idx3) < vertices(1,exon_idx2)) & ...
             (min(edges(exon_idx3,exon_idx3+1:end)==0)==1) & ...
             (exon_idx1~=exon_idx3) & (exon_idx2~=exon_idx3) & ...
             (vertices(1,exon_idx3) < vertices(2,exon_idx1))
            if genes(ix).strands(1)=='+',
              idx_exon_intron5 = [idx_exon_intron5, ix];
              exon_exon_intron5 = [exon_exon_intron5, [exon_idx1; exon_idx2; exon_idx3]];
            elseif genes(ix).strands(1)=='-',
              idx_exon_intron3 = [idx_exon_intron3, ix];
              exon_exon_intron3 = [exon_exon_intron3, [exon_idx1; exon_idx2; exon_idx3]];
            end
          end
          if (vertices(1,exon_idx3) > vertices(2,exon_idx1)) & ...
             (vertices(1,exon_idx3) < vertices(1,exon_idx2)) & ...
             (min(edges(1:exon_idx3-1, exon_idx3)==0)==1) & ...
             (exon_idx1~=exon_idx3) & (exon_idx2~=exon_idx3) & ...
             (vertices(2,exon_idx3) > vertices(1,exon_idx2))
            if genes(ix).strands(1)=='+',
              idx_exon_intron3 = [idx_exon_intron3, ix];
              exon_exon_intron3 = [exon_exon_intron3, [exon_idx1; exon_idx2; exon_idx3]];
            elseif genes(ix).strands(1)=='-',
              idx_exon_intron5 = [idx_exon_intron5, ix];
              exon_exon_intron5 = [exon_exon_intron5, [exon_idx1; exon_idx2; exon_idx3]];
            end
          end
        end
      end
    end
  end
end
fprintf(1,'\n\nNumber of incomplete exons in intronic regions (5ps):\t\t%d\n',...
	length(idx_exon_intron5));
fprintf(1,'Number of incomplete exons in intronic regions (3ps):\t\t%d\n',...
	length(idx_exon_intron3));

