function [idx_xor_exons, exon_xor_exons,id_xor_exons] = detect_xorexons(genes, idx_alt);

idx_xor_exons = zeros(1,0); 
id_xor_exons = zeros(1,0); 
exon_xor_exons = zeros(4,0); %%%% 5primesite of first exon, the 2 skipped
                     %exons, 3primesite of last exon %%%
id = 0 ;
for ix=idx_alt
  if (mod(ix,50)==0)
    fprintf(1,'.');
  end
  num_exons = size(genes(ix).splicegraph{1},2);
  edges = genes(ix).splicegraph{2};
  vertices = genes(ix).splicegraph{1};
  for exon_idx1 = 1:num_exons-3
    id=id+1 ;
    for exon_idx2 = (exon_idx1+1):(num_exons-2)
      if edges(exon_idx1, exon_idx2)
        for exon_idx3=(exon_idx2+1):(num_exons-1)
          if edges(exon_idx1, exon_idx3) & edges(exon_idx2, exon_idx3)==0 & ...
             (vertices(1,exon_idx3) > vertices(2,exon_idx2)),
            for exon_idx4=(exon_idx3+1):num_exons
              if edges(exon_idx2, exon_idx4) & edges(exon_idx3, exon_idx4) & ...
                    edges(exon_idx1, exon_idx4)==0
                idx_xor_exons = [idx_xor_exons, ix];
                id_xor_exons = [id_xor_exons, id];
                exon_xor_exons = [exon_xor_exons,[exon_idx1; exon_idx2; exon_idx3; exon_idx4]] ;
              end
            end
          end
        end
      end
    end
  end
end
fprintf(1,'\n\nNumber of XOR exons:\t\t\t\t\t\t%d\n',...
	length(idx_xor_exons));


