function genes = reduce_splice_graph(genes) ;
% genes = reduce_splice_graph(genes) ;
%
% Iterates over all genes and removes complexity
% from the splice graph by:
% - collapsing identical exons into a single node
% - merging one single exon transcripts containing another 
%   single exon transcript into one
% - collapsing alternative transcript starts into a single
%   start if they share 3' exon boundary
% - collapsing alternative transcript ends into a single
%   end if they share 5' exon boundary

for gene_idx = 1:length(genes)

  vertices = genes(gene_idx).splicegraph{1};
  edges = genes(gene_idx).splicegraph{2};

  if (mod(gene_idx,1000) == 0)
    fprintf(1,'%d\n', gene_idx);
  end

  %%% no vertices in the splice graph
  if isempty(vertices), continue ; end ;

  %%% find all the intron locations
  [dummy, exon_order] = sort(vertices(1,:),2,'ascend');
  vertices = vertices(:, exon_order);
  edges = edges(exon_order, exon_order);
  intron_loc = [];
  for ix1 = 1:size(vertices, 2) - 1,
    for ix2 = 2:size(vertices, 2),
      if edges(ix1, ix2),
        intron_loc = union(intron_loc, [vertices(2, ix1) + 1 : vertices(1, ix2) - 1]);
      end;
    end;
  end;
  
  %%% if one or two vertices (one or no edge) exist 
  if (size(edges, 1) < 2)
    changed = 1;
    while changed
      changed = 0;
      exon_idx = 1;
      while exon_idx <= size(vertices, 2),
        test_exon_idx = exon_idx+1;
        while test_exon_idx <= size(vertices,2),
          %%%  <------- exon_idx -------> 
          %%%    <-- test_exon_idx -->
          if (vertices(1, exon_idx) <= vertices(1, test_exon_idx)) &&...
             (vertices(2, exon_idx) >= vertices(2, test_exon_idx))

             %%% keep longer exon
             new_index = [1:test_exon_idx-1, test_exon_idx+1:size(vertices,2)];
             vertices = vertices(:,new_index);
      
             changed = 1;
          end;
          test_exon_idx = test_exon_idx + 1;
        end;
        exon_idx = exon_idx + 1;
      end;
    end;
  %%% more than two vertices exist
  else
    changed = 1;
    while changed,
      changed = 0;
      exon_idx = 1;
      while exon_idx <= size(vertices,2),
        %% re-sort vertices if necessary
        [vertices_sorted, exon_order] = sort(vertices(1,:), 2, 'ascend');
        if ~isequal(vertices(1,:), vertices_sorted),
          vertices = vertices(:, exon_order);
          edges = edges(exon_order, exon_order);
        end;
    
        test_exon_idx = exon_idx + 1;
        while test_exon_idx <= size(vertices,2)
          reduce_now = 0;

          if (test_exon_idx < exon_idx) && keyboard_allowed(), keyboard; end;
        
          %%% count incoming and outgoing edges for exon and test_exon
          cur_edge_left = sum(edges(1:exon_idx,exon_idx));
          test_edge_left = sum(edges(1:test_exon_idx,test_exon_idx));
          cur_edge_right = sum(edges(exon_idx:end,exon_idx));
          test_edge_right = sum(edges(test_exon_idx:end,test_exon_idx));
      
          %%% 0000
          %%% no incoming or outgoing edges in exon and test_exon
          if (~cur_edge_left && ~cur_edge_right && ~test_edge_left && ~test_edge_right),
            %%%     <------ exon ------->>>>>>>>>>>>>  OR          <------ exon ------>
            %%%              <--- test_exon ---->           <---- test_exon ---->>>>>>>>>>>>
            if ((vertices(2, exon_idx) >= vertices(1, test_exon_idx)) && (vertices(1, exon_idx) <= vertices(1, test_exon_idx))) || ...
               ((vertices(2, test_exon_idx) >= vertices(1, exon_idx)) && (vertices(1, test_exon_idx) <= vertices(1, exon_idx))) && ...
               (sum(ismember([min(vertices(1, exon_idx), vertices(1, test_exon_idx)) : max(vertices(2, exon_idx), vertices(2, test_exon_idx))], intron_loc)) == 0),

              %%% merge exons if they overlap and they do not span any intronic position
              vertices(1, exon_idx) = min(vertices(1, exon_idx), vertices(1, test_exon_idx));
              vertices(2, exon_idx) = max(vertices(2, exon_idx), vertices(2, test_exon_idx));
              new_index = [1:test_exon_idx - 1, test_exon_idx + 1:size(vertices, 2)];
            
              vertices = vertices(:, new_index);
              edges = edges(new_index, new_index); % no need to combine any adges, as both exons have a degree of 0
            
              reduce_now = 1;
              changed = 1;
            end
        
          %%% 0101
          %%% outgoing edges in exon and test_exo, no incoming edges
          elseif (~cur_edge_left && cur_edge_right && ~test_edge_left && test_edge_right),
            %%%   ----- exon -----<
            %%%   --- test_exon --<
            if (vertices(2, exon_idx) == vertices(2, test_exon_idx)) && ...
               (sum(ismember([min(vertices(1, exon_idx), vertices(1, test_exon_idx)) : vertices(2, exon_idx)], intron_loc)) == 0),
              
              %%% merge exons if they share the same right boundary and do not span intronic positions
              vertices(1, exon_idx) = min(vertices(1, exon_idx), vertices(1,test_exon_idx));
              new_index = [1:test_exon_idx - 1, test_exon_idx + 1:size(vertices, 2)];
            
              vertices = vertices(:,new_index);
              edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
              edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
              edges = edges(new_index,new_index);
            
              reduce_now = 1;
              changed = 1;
            end

          %%% 1010
          %%% incoming edges in exon and test_exon, no outgoing edges
          elseif (cur_edge_left && ~cur_edge_right && test_edge_left &&~ test_edge_right)
            %%%   >---- exon ------
            %%%   >-- test_exon ---
            if (vertices(1, exon_idx) == vertices(1, test_exon_idx)) && ...
               (sum(ismember([vertices(1, exon_idx) : max(vertices(2, exon_idx), vertices(2, test_exon_idx))], intron_loc)) == 0)
              
              %%% merge exons if they share the same left boundary and do not span intronic positions
              vertices(2, exon_idx) = max(vertices(2, exon_idx), vertices(2, test_exon_idx));
              new_index = [1 : test_exon_idx - 1, test_exon_idx + 1 : size(vertices, 2)];
            
              vertices = vertices(:, new_index);
              edges(exon_idx, :) = or(edges(exon_idx, :), edges(test_exon_idx, :));
              edges(:, exon_idx) = or(edges(:, exon_idx), edges(:, test_exon_idx));
              edges = edges(new_index, new_index);
            
              reduce_now = 1;
              changed = 1;
            end

          %%% 1111
          %%% exon and test_exon have both incoming and outgoing edges
          elseif (cur_edge_left && cur_edge_right && test_edge_left && test_edge_right)
            %%%  >------ exon -----<
            %%%  >--- test_exon ---<
            if (vertices(1, exon_idx) == vertices(1, test_exon_idx)) && ...
               (vertices(2, exon_idx) == vertices(2, test_exon_idx))
            
              %%% collapse identical exons into one node
              new_index = [1 : test_exon_idx - 1, test_exon_idx + 1 : size(vertices, 2)];
            
              vertices = vertices(:, new_index);
              edges(exon_idx, :) = or(edges(exon_idx, :), edges(test_exon_idx, :));
              edges(:, exon_idx) = or(edges(:, exon_idx), edges(:, test_exon_idx));
              edges = edges(new_index, new_index);
              
              reduce_now = 1;
              changed = 1;
            end
          end % cases

          if ~reduce_now
            test_exon_idx = test_exon_idx + 1;
          end
        end
        exon_idx = exon_idx + 1;
      end %while exon_idx <= size(vertices,2)
    end %while changed
  end % if (size(edges,1) < 2) 
  genes(gene_idx).splicegraph = {};
  genes(gene_idx).splicegraph = {vertices,edges};
end % end for

fprintf(1,'\n');

return
