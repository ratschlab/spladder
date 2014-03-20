function genes = build_splice_graph(genes, g_idx) ;
% genes = build_splice_graph(genes) ;
%
% This script takes the given gene structure and
% constructs a splice graph for each transcript
% based on its annotated exons.
%
% !!! An already existing splicegraph in genes
% !!! is replaced by the newly constructed one!
%

if nargin < 2,
    g_idx = 1:length(genes);
end;

for gene_idx = g_idx,
  % construct a new empty splice graph, based on
  % http://proline.bic.nus.edu.sg/dedb/
  vertices =  [];
  edges = [];
  terminals = [];
  if (mod(gene_idx,100)==0)
    fprintf(1,'.');
    if (mod(gene_idx, 1000)==0),
        fprintf(1, '%i\n', gene_idx);
    end;
  end

  for transcript_idx = 1:length(genes(gene_idx).transcripts)
    exon_start_end = genes(gene_idx).exons{transcript_idx};
    
    %%% only one exon in the transcript
    if (size(exon_start_end,1) == 1)
      exon1_start = exon_start_end(1,1);
      exon1_end = exon_start_end(1,2);

      if isempty(vertices)
        vertices(1,1) = exon1_start;
        vertices(2,1) = exon1_end;
        edges = 0;
        num_exons = 1;
      else
        vertices(1,num_exons+1) = exon1_start;
        vertices(2,num_exons+1) = exon1_end;
        edges(1,num_exons+1) = 0;
        edges(num_exons+1,1) = 0;
        num_exons = num_exons + 1;
      end
    %%% more than one exon in the transcript
    else
      for exon_idx = 1:(size(exon_start_end,1)-1)
        exon1_start = exon_start_end(exon_idx,1);
        exon1_end = exon_start_end(exon_idx,2);
        exon2_start = exon_start_end(exon_idx+1,1);
        exon2_end = exon_start_end(exon_idx+1,2);
  
        if isempty(vertices)
          vertices(1,1) = exon1_start;
          vertices(2,1) = exon1_end;
          vertices(1,2) = exon2_start;
          vertices(2,2) = exon2_end;
          edges = zeros(2);
          edges(1,2) = 1;
          edges(2,1) = 1;
          num_exons = 2;
        else
          exon1_idx = 0;
          exon2_idx = 0;
          %%% check if current exon already occurred
          for idx = 1:num_exons,
            if ((vertices(1,idx)==exon1_start) &&...
               (vertices(2,idx)==exon1_end))
               exon1_idx = idx;
            end
            if ((vertices(1,idx)==exon2_start) &&... 
               (vertices(2,idx)==exon2_end))
               exon2_idx = idx;
            end
          end
          %%% both exons already occured -> only add an edge
          if (exon1_idx~=0) && (exon2_idx~=0)
            edges(exon1_idx,exon2_idx) = 1;
            edges(exon2_idx,exon1_idx) = 1;
          else
            %%% 2nd exon occured
            if ((exon1_idx==0) && (exon2_idx~=0))
              vertices(1,num_exons+1) = exon1_start;
              vertices(2,num_exons+1) = exon1_end;
              edges(exon2_idx,num_exons+1) = 1;
              edges(num_exons+1,exon2_idx) = 1;
              num_exons = num_exons + 1;
            %%% 1st exon occured
            elseif ((exon2_idx==0) && (exon1_idx~=0))
              vertices(1,num_exons+1) = exon2_start;
              vertices(2,num_exons+1) = exon2_end;
              edges(exon1_idx,num_exons+1) = 1;
              edges(num_exons+1,exon1_idx) = 1;
              num_exons = num_exons + 1;
            %%% no exon occured
            else
              assert((exon1_idx==0) && (exon2_idx==0));
              vertices(1,num_exons+1) = exon1_start;
              vertices(2,num_exons+1) = exon1_end;
              num_exons = num_exons + 1;
              vertices(1,num_exons+1) = exon2_start;
              vertices(2,num_exons+1) = exon2_end;
              num_exons = num_exons + 1;    
              
              edges(num_exons-1,num_exons) = 1;
              edges(num_exons,num_exons-1) = 1;
            end
          end
        end
      end
      terminals(1, find(vertices(1, :) == min(min(exon_start_end)))) = 1;
      terminals(2, find(vertices(2, :) == max(max(exon_start_end)))) = 1;
    end
  end
  %%% take care of the sorting by exon start
  [tmp, s_idx] = sort(vertices(1, :));
  vertices = vertices(:, s_idx);
  edges = edges(s_idx,s_idx);
  if ~isequal(size(vertices), size(terminals)),
    terminals(size(vertices, 1), size(vertices, 2)) = 0;
  end;
  terminals = terminals(:, s_idx);
  genes(gene_idx).splicegraph = {vertices,edges, terminals};
end

fprintf(1,'\n');

return
