function [idx_alt_intron, introns_alt_intron] = detect_altintrons(genes, idx_alt) ;

idx_alt_intron = [];
introns_alt_intron = [];
for ix=idx_alt
  if (mod(ix,50)==0)
    fprintf(1,'.');
    %fprintf('%i\n', ix);
  end
  num_exons = size(genes(ix).splicegraph{1},2) ;
  vertices = genes(ix).splicegraph{1} ;
  edges = genes(ix).splicegraph{2} ;
  
  for exon_idx1=1:num_exons-3,
    for exon_idx2=exon_idx1+1:num_exons-2,
      
      version=0 ;
      if (vertices(2,exon_idx1) > vertices(1,exon_idx2)) && ...
         (vertices(2,exon_idx1) < vertices(2,exon_idx2)),
        version=1 ;
      elseif (vertices(2,exon_idx2) > vertices(1,exon_idx1)) && ...
             (vertices(2,exon_idx2) < vertices(2,exon_idx1)),
        version=2 ;
      end
      
      if version
        for exon_idx3=find(edges(exon_idx1,exon_idx1+1:end)==1)+exon_idx1
          for exon_idx4=find(edges(exon_idx2,exon_idx2+1:end)==1)+exon_idx2
            if exon_idx3==exon_idx4,
              continue ;
            end
            exists=0 ;
            if version==1
              if (vertices(1,exon_idx3) > vertices(1,exon_idx4)) && ...
                 (vertices(1,exon_idx3) < vertices(2,exon_idx4))
                version=11 ;
              elseif (vertices(1,exon_idx4) > vertices(1,exon_idx3)) && ...
                     (vertices(1,exon_idx4) < vertices(2,exon_idx3))
                version=12 ;
              else
                version=0 ;
              end
            end
            if version==2
              if (vertices(1,exon_idx3) > vertices(1,exon_idx4)) && ...
                 (vertices(1,exon_idx3) < vertices(2,exon_idx4))
                version=21 ;
              elseif (vertices(1,exon_idx4) > vertices(1,exon_idx3)) && ...
                     (vertices(1,exon_idx4) < vertices(2,exon_idx3))
                version=22 ;
              else
                version=0 ;
              end
            end
            if version
              for i=1:length(idx_alt_intron),
                if idx_alt_intron(i)==ix,
                  if (vertices(2,exon_idx1)==vertices(2,introns_alt_intron(1,i)) && ...
                      vertices(2,exon_idx2)==vertices(2,introns_alt_intron(3,i)) && ...
                      vertices(1,exon_idx3)==vertices(1,introns_alt_intron(2,i)) && ...
                      vertices(1,exon_idx4)==vertices(1,introns_alt_intron(4,i))) || ...    
                     (vertices(2,exon_idx1)==vertices(2,introns_alt_intron(3,i)) && ...
                      vertices(2,exon_idx2)==vertices(2,introns_alt_intron(1,i)) && ...
                      vertices(1,exon_idx3)==vertices(1,introns_alt_intron(4,i)) && ...
                      vertices(1,exon_idx4)==vertices(1,introns_alt_intron(2,i)))
                    exists=1 ;
                  end    
                end
              end
              if ~exists
                if edges(exon_idx1, exon_idx4) && edges(exon_idx2, exon_idx3),
                  %fprintf('\nhere: %i, idx1: %i idx2: %i idx3: %i idx4: %i\n', ...
                  %        ix, exon_idx1, exon_idx2, exon_idx3, exon_idx4) ;
                end
                idx_alt_intron=[idx_alt_intron, ix] ;
                introns_alt_intron=[introns_alt_intron, [exon_idx1; exon_idx3; ...
                                                         exon_idx2 ; exon_idx4; ...
                                                         version]] ;
              end
            end
          end
        end
      end
    end
  end
end
        

fprintf(1,'\n\nNumber of alternative introns:\t\t\t\t\t%d\n',...
	length(idx_alt_intron));

