function detect_rule(genes, idx_alt) ;

addpath ../../utils
addpath ../../splicegraphs
addpath /fml/ag-raetsch/share/software/matlab_tools/utils
  
for ix=idx_alt
  if (mod(ix,50)==0)
      fprintf(1,'.');
  end
  num_exons = size(genes(ix).splicegraph{1},2) ;
  vertices = genes(ix).splicegraph{1} ;
  edges = genes(ix).splicegraph{2} ;
  if genes(ix).strands(1)=='+'
    for i=1:num_exons-1,
      for j=i+1:num_exons,
        if vertices(2,i)==vertices(2,j)
          %disp('found 3 primesites');
          %disp(fprintf('Gene %i, strand: %s, node %i and %i', ix, ...
          %                 genes(ix).strands(1), i, j)),
          %viewsplicegraph(genes(ix)) ;
          %pause ;
          for k=i+1:num_exons
            if (edges(i,k)==1 & edges(j,k)~=1) & (edges(i,k)~=1 & edges(j,k)==1)
              disp(fprintf('Gene %i, strand: %s, node %i and %i', ix, ...
                           genes(ix).strands(1), i, j)),
              keyboard ;
            else
              %disp('ok')
            end
          end
        end
      end
    end
  end
  if genes(ix).strands(1)=='-'
    for i=2:num_exons,
      for j=1:i-1,
        if vertices(1,i)==vertices(1,j)
          %disp('found 5 primesites');
          %viewsplicegraph(genes(ix)) ;
          %disp(fprintf('Gene %i, strand: %s, node %i and %i', ix, ...
          %                 genes(ix).strands(1), i, j)),
          %pause ;
          for k=1:i-1
            if (edges(k,i)==1 & edges(k,j)~=1) & (edges(k,i)~=1 & edges(k,j)==1)
              disp(fprintf('Gene %i, strand: %s, node %i and %i', ix, ...
                           genes(ix).strands(1), i, j)),
              keyboard ;
            else
              %disp('ok')
            end
          end
        end
      end
    end 
  end
  %pause ;
end
