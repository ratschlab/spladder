function genes = infer_splice_graph_caller(genes) ;
% genes = infer_splice_graph_caller(genes) ;
  
%infer_splice_graph

MISALIGN = 5;

fprintf(1,'Number of genes:%d\n',length(genes));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some simple merging to reduce exons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'Merging some genes\n');
genes = reduce_splice_graph_caller(genes) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% infer_splice_graph infers the relevant splice graph given a set of ESTs. 
%
% It is based on two assumptions, neither which are true in reality, but
% are necessary for any reasonable processing.
% 1. All splicing events are independent. That is a splice site
%    at location A in ESTA does not depend on splice site B in ESTA, and
%    neither does it depend on splice site C in ESTB.
% 2. There exists a splicing event at a particular location if and only if
%    there exists an EST that is spliced at this point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf(1,'performing inference\n');

for gene_idx = 1:length(genes)

  if (mod(gene_idx,1000)==0)
    fprintf(1,'%d\n',gene_idx);
  end
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For each gene, we consider each pair of exons in the splicegraph.
  % After sorting the exons by the location of the left side, we make
  % a local copy of the splicegraph.
  %
  % Note that 'left' means the 5 prime end of a positive strand, but
  % the 3 prime end of a negative strand.
  %
  % The two exons considered are indexed by exon_idx, and test_exon_idx.
  % The general idea of the algorithm is as follows:
  % For each pair:
  %   check which case it is in (more on this later)
  %   if it is possible to merge:
  %     form a new exon
  %   end merge.
  %   if new exon is not in current set
  %     add to splicegraph
  %   end new exon insert
  % end for
  % 
  % In words:
  % Each time two exons are considered, and they can be merged to form
  % a new exon, they are added to the list of exons. Then, each of the
  % original parent exons are considered for deletion (variable to_keep).
  % All the exons are kept till the very end, when all exons labelled
  % to_keep == 0 are removed.
  %
  % Hence the number of exons could potentially grow quadratically.
  % Since we have to also consider merging the pairs of exons which are
  % newly created, the number of pairs considered can be larger than
  % m^2 where m is the number of original exons. However, I believe that
  % the algorithm should terminate because there are only finitely many
  % splice sites and hence at worst, we have all possible left and right
  % exon ends used in our splicegraph.
  %
  % I am not sure what the best termination condition for the loops are.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if isempty(genes(gene_idx).splicegraph{1}),
    continue ; 
  end ;
  [dummy,exon_order] = sort(genes(gene_idx).splicegraph{1}(1,:),2,'ascend');
  vertices = genes(gene_idx).splicegraph{1}(:,exon_order);
  edges = genes(gene_idx).splicegraph{2}(exon_order,exon_order);
  to_keep = ones(1,size(vertices,2));

  changed = 0;
  exon_idx = 1;
  first_merge_exon = 1;
  test_exon_idx = 0;
  while exon_idx <= size(vertices,2)

    if changed
      changed = 0;
      exon_idx = first_merge_exon;
      continue;
    end
    
    num_exons = size(vertices,2);
    %%% sort exons ascending if neccessary
    [dummy,exon_order] = sort(vertices(1,:),2,'ascend');
    if ~isequal(vertices(1,:),dummy)
      vertices = vertices(:,exon_order);
      edges = edges(exon_order,exon_order);
      to_keep = to_keep(exon_order);
    end ;
    
    %%% take all exons overlapping to exon_idx +- MISALIGN window
    exon_take_idx = find(((vertices(1,exon_idx)-2*MISALIGN)<vertices(2,:))&...
       ((vertices(2,exon_idx)+2*MISALIGN)>vertices(1,:)));
    if ~isempty(exon_take_idx)
      first_merge_exon = exon_take_idx(1);
    end
    %%% consider only exons downstream of exon_idx
    exon_take_idx = exon_take_idx(exon_take_idx>exon_idx);
    internal_idx = 1;

    while ((internal_idx <= length(exon_take_idx)) && ~changed)
      test_exon_idx = exon_take_idx(internal_idx);
      if (test_exon_idx < exon_idx) && keyboard_allowed(), keyboard;end;
      cur_edge_left = sum(edges(1:exon_idx,exon_idx));
      test_edge_left = sum(edges(1:test_exon_idx,test_exon_idx));
      cur_edge_right = sum(edges(exon_idx:end,exon_idx));
      test_edge_right = sum(edges(test_exon_idx:end,test_exon_idx));
    
      new_vertex = zeros(2,1);
  
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % In the following, the cases are labelled with a binary string
      % [cur_edge_left, cur_edge_right, test_edge_left, test_edge_right]
      %
      % cur refers to exon_idx
      % test refers to test_exon_idx
      % edge means introns.
      % e.g. cur_edge_left == 1 means that there exists an intron on
      % the left of exon_idx
      %
      % Hence there are 16 cases, and we treat each one individually.
      % For each case, the two logical tests are:
      % 1. should we merge?
      % 2. is the new exon known?
      %
      % Note that a splice site is defined by an intron-exon boundary.
      % Furthermore, we allow misalignment by MISALIGN nucleotides,
      % that is we shorten an exonAendA if it doesn't have a splice site
      % and there exists a splice site exonBendB within MISALIGN nucleotides.
      % 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      % 0000
      if (~cur_edge_left&&~cur_edge_right&&~test_edge_left&&~test_edge_right)
        %%% if exon_idx and test_exon_idx differ in at least one border and
        %%% at least one of both is not yet deleted and 
        %%%   --- exon ---X               X >=Y
        %%%         Y--- test_exon ---  
        if ((~isequal(vertices(1,exon_idx),vertices(1,test_exon_idx)))||...
            (~isequal(vertices(2,exon_idx),vertices(2,test_exon_idx))))&&...
              ((to_keep(test_exon_idx)~=0)||(to_keep(exon_idx)~=0))&&...
              (vertices(2,exon_idx)>=vertices(1,test_exon_idx))
          
          new_vertex(1) = min(vertices(1,exon_idx),vertices(1,test_exon_idx));
          new_vertex(2) = max(vertices(2,exon_idx),vertices(2,test_exon_idx));
          new_edges = zeros(1,num_exons);
        
          known = 0;
          for ix = 1:size(vertices,2)
            if isequal(new_vertex,vertices(:,ix))
              known = ix;
              break;
            end
          end
        
          if ~known
            vertices = [vertices,new_vertex];
            edges = [edges;new_edges];
            edges = [edges,[new_edges';0]];
            to_keep(end+1) = 1;
            to_keep(exon_idx) = 0;
            to_keep(test_exon_idx) = 0;
            changed = 1;
          else
            if (known ~= exon_idx) && to_keep(exon_idx)
              to_keep(exon_idx) = 0;
              changed = 1;
            end
            if (known ~= test_exon_idx) && to_keep(test_exon_idx)
              to_keep(test_exon_idx) = 0;
              changed = 1;
            end
          end
        end
        
      % 0001
      elseif (~cur_edge_left&&~cur_edge_right&&~test_edge_left&&test_edge_right)
        %%%   --- exon ---M
        %%%         --- test_exon ---<
        if (vertices(2,exon_idx)-MISALIGN<=vertices(2,test_exon_idx))&&...
            (vertices(2,exon_idx)>=vertices(1,test_exon_idx))
      
          new_vertex(1) = min(vertices(1,exon_idx),vertices(1,test_exon_idx));
          new_vertex(2) = vertices(2,test_exon_idx);
          new_edges = or(edges(exon_idx,:),edges(test_exon_idx,:));
        
          known = 0;
          for ix = 1:size(vertices,2)
            if isequal(new_vertex,vertices(:,ix)) && issubset(new_edges,edges(ix,:)) 
              known = ix;
              break;
            end
          end
          assert(known~=exon_idx);
        
          if to_keep(exon_idx)
            to_keep(exon_idx) = 0;
            changed = 1;
          end
        
          if ~known
            vertices = [vertices,new_vertex];
            edges = [edges;new_edges];
            edges = [edges,[new_edges';0]];
            to_keep(end+1) = 1;
            to_keep(exon_idx) = 0;
            to_keep(test_exon_idx) = 0;
            changed = 1;
          elseif (known ~= test_exon_idx) && to_keep(test_exon_idx)
            to_keep(test_exon_idx) = 0;
            changed = 1;
          end
        end

      % 0100
      elseif (~cur_edge_left&&cur_edge_right&&~test_edge_left&&~test_edge_right)
        %%%            --- exon ---<
        %%%  --- test_exon ---M
        if (vertices(2,exon_idx)>=vertices(2,test_exon_idx)-MISALIGN)&&...
           (vertices(1,exon_idx)<=vertices(2,test_exon_idx))
      
        new_vertex(1) = min(vertices(1,exon_idx),vertices(1,test_exon_idx));
        new_vertex(2) = vertices(2,exon_idx);
        new_edges = or(edges(exon_idx,:),edges(test_exon_idx,:));
      
        known = 0;
        for ix = 1:size(vertices,2)
          if isequal(new_vertex,vertices(:,ix)) && issubset(new_edges,edges(ix,:))
            known = ix;
            break;
          end
        end
        assert(known~=test_exon_idx);
      
        if to_keep(test_exon_idx)
          to_keep(test_exon_idx) = 0;
          changed = 1;
        end
      
        if ~known
          vertices = [vertices,new_vertex];
          edges = [edges;new_edges];
          edges = [edges,[new_edges';0]];
          to_keep(end+1) = 1;
          to_keep(exon_idx) = 0;
          to_keep(test_exon_idx) = 0;
          changed = 1;
        elseif (known ~= exon_idx) && to_keep(exon_idx)
          to_keep(exon_idx) = 0;
          changed = 1;
        end
      end
  
      % 0010
      elseif (~cur_edge_left&&~cur_edge_right&&test_edge_left&&~test_edge_right)
        %%%            M--- exon ---
        %%% >--- test_exon ---
        if (vertices(1,exon_idx)+MISALIGN>=vertices(1,test_exon_idx))&&...
              (vertices(1,exon_idx)<=vertices(2,test_exon_idx))
            
          new_vertex(1) = vertices(1,test_exon_idx);
          new_vertex(2) = max(vertices(2,exon_idx),vertices(2,test_exon_idx));
          new_edges = or(edges(exon_idx,:),edges(test_exon_idx,:));
        
          known = 0;
          for ix = 1:size(vertices,2)
            if isequal(new_vertex,vertices(:,ix)) && issubset(new_edges,edges(ix,:))
              known = ix;
              break;
            end
          end
          assert(known~=exon_idx);
        
          if to_keep(exon_idx)
            to_keep(exon_idx) = 0;
            changed = 1;
          end
        
          if ~known
            vertices = [vertices,new_vertex];
            edges = [edges;new_edges];
            edges = [edges,[new_edges';0]];
            to_keep(end+1) = 1;
            to_keep(exon_idx) = 0;
            to_keep(test_exon_idx) = 0;
            changed = 1;
          elseif (known ~= test_exon_idx) && to_keep(test_exon_idx)
            to_keep(test_exon_idx) = 0;
            changed = 1;
          end
        end
  
      % 1000
      elseif (cur_edge_left&&~cur_edge_right&&~test_edge_left&&~test_edge_right)
        %%%  >--- exon ---
        %%%      M--- test_exon ---
        if (vertices(1,exon_idx)<=vertices(1,test_exon_idx)+MISALIGN)&&...
              (vertices(2,exon_idx)>=vertices(1,test_exon_idx))
      
          new_vertex(1) = vertices(1,exon_idx);
          new_vertex(2) = max(vertices(2,exon_idx),vertices(2,test_exon_idx));
          new_edges = or(edges(exon_idx,:),edges(test_exon_idx,:));
        
          known = 0;
          for ix = 1:size(vertices,2)
            if isequal(new_vertex,vertices(:,ix)) && issubset(new_edges,edges(ix,:))
              known = ix;
              break;
            end
          end
          assert(known~=test_exon_idx);
        
          if to_keep(test_exon_idx)
            to_keep(test_exon_idx) = 0;
            changed = 1;
          end
        
          if ~known
            vertices = [vertices,new_vertex];
            edges = [edges;new_edges];
            edges = [edges,[new_edges';0]];
            to_keep(end+1) = 1;
            to_keep(exon_idx) = 0;
            to_keep(test_exon_idx) = 0;
            changed = 1;
          elseif (known ~= exon_idx) && to_keep(exon_idx)
            to_keep(exon_idx) = 0;
            changed = 1;
          end
        end
      
      % 0011
      elseif (~cur_edge_left&&~cur_edge_right&&test_edge_left&&test_edge_right)
        %%%     M--- exon ---M
        %%%  >--- test_exon ---<
        if (vertices(1,exon_idx)+MISALIGN>=vertices(1,test_exon_idx))&&...
            (vertices(2,exon_idx)-MISALIGN<=vertices(2,test_exon_idx))
      
          new_vertex = vertices(:,test_exon_idx);
          new_edges = or(edges(exon_idx,:),edges(test_exon_idx,:));
          
          known = 0;
          for ix = 1:size(vertices,2)
            if isequal(new_vertex,vertices(:,ix)) && issubset(new_edges,edges(ix,:))
              known = ix;
              break;
            end
          end
          assert(known~=exon_idx);
        
          if to_keep(exon_idx)
            to_keep(exon_idx) = 0;
            changed = 1;
          end
        
          if ~known
            vertices = [vertices,new_vertex];
            edges = [edges;new_edges];
            edges = [edges,[new_edges';0]];
            to_keep(end+1) = 1;
            to_keep(exon_idx) = 0;
            to_keep(test_exon_idx) = 0;
            changed = 1;
          elseif (known ~= test_exon_idx) && to_keep(test_exon_idx)
            to_keep(test_exon_idx) = 0;
            changed = 1;
          end
        end

      % 1100
      elseif (cur_edge_left&&cur_edge_right&&~test_edge_left&&~test_edge_right)
        %%%    >-------- exon ---------<
        %%%       M--- test_exon ---M
        if (vertices(1,exon_idx)<=vertices(1,test_exon_idx)+MISALIGN)&&...
            (vertices(2,exon_idx)>=vertices(2,test_exon_idx)-MISALIGN)
      
          new_vertex = vertices(:,exon_idx);
          new_edges = or(edges(exon_idx,:),edges(test_exon_idx,:));
          
          known = 0;
          for ix = 1:size(vertices,2)
            if isequal(new_vertex,vertices(:,ix)) && issubset(new_edges,edges(ix,:))
              known = ix;
              break;
            end
          end
          assert(known~=test_exon_idx);
        
          if to_keep(test_exon_idx)
            to_keep(test_exon_idx) = 0;
            changed = 1;
          end
        
          if ~known
            vertices = [vertices,new_vertex];
            edges = [edges;new_edges];
            edges = [edges,[new_edges';0]];
            to_keep(end+1) = 1;
            to_keep(exon_idx) = 0;
            to_keep(test_exon_idx) = 0;
            changed = 1;
          elseif (known ~= exon_idx) && to_keep(exon_idx)
            to_keep(exon_idx) = 0;
            changed = 1;
          end
        end
  
      % 0101
      elseif (~cur_edge_left&&cur_edge_right&&~test_edge_left&&test_edge_right)
        %%%    -------- exon ------<
        %%%       --- test_exon ---<
        if (vertices(2,exon_idx)==vertices(2,test_exon_idx))
          %%% same edges
          if (~isequal(edges(exon_idx,:),edges(test_exon_idx,:)))
            edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
            edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
            edges(test_exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
            edges(:,test_exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
            %%% keep longer exon
            if (vertices(1,exon_idx) > vertices(1,test_exon_idx))
              to_keep(exon_idx) = 0;
            else
              to_keep(test_exon_idx) = 0;
            end
            changed = 1;
          %%% different edges
          else
            %%%        ---- exon ------<
            %%%   ------- test_exon ---<
            if (vertices(1,exon_idx) > vertices(1,test_exon_idx)) && (to_keep(exon_idx))
              to_keep(exon_idx) = 0;
              changed = 1;
            end
            %%%   --------- exon ------<
            %%%       --- test_exon ---<
            if (vertices(1,exon_idx) <= vertices(1,test_exon_idx)) && (to_keep(test_exon_idx))
              to_keep(test_exon_idx) = 0;
              changed = 1;
            end
          end
        end
  
      % 1010
      elseif (cur_edge_left&&~cur_edge_right&&test_edge_left&&~test_edge_right)
        %%%    >-------- exon ------
        %%%    >--- test_exon ---
        if (vertices(1,exon_idx)==vertices(1,test_exon_idx))
          %%% same edges
          if (~isequal(edges(exon_idx,:),edges(test_exon_idx,:)))
            edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
            edges(:,exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
            edges(test_exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
            edges(:,test_exon_idx) = or(edges(:,exon_idx),edges(:,test_exon_idx));
            %%% keep longer exon
            if (vertices(2,exon_idx) < vertices(2,test_exon_idx))
              to_keep(exon_idx) = 0;
            else
              to_keep(test_exon_idx) = 0;
            end
            changed = 1;
          %%% different edges
          else
            %%%    >----- exon ------
            %%%    >------ test_exon ---
            if (vertices(2,exon_idx) < vertices(2,test_exon_idx)) && (to_keep(exon_idx))
              to_keep(exon_idx) = 0;
              changed = 1;
            end
            %%%    >-------- exon ------
            %%%    >--- test_exon ---
            if (vertices(2,exon_idx) >= vertices(2,test_exon_idx)) && (to_keep(test_exon_idx))
              to_keep(test_exon_idx) = 0;
              changed = 1;
            end
          end
        end
  
      % 0110  
      elseif (~cur_edge_left&&cur_edge_right&&test_edge_left&&~test_edge_right)
        %%%          MX----- exon -----<    X<=Y
        %%%    >--- test_exon ---YM
        if (vertices(1,exon_idx)+MISALIGN>=vertices(1,test_exon_idx))&&...
            (vertices(2,exon_idx)>=vertices(2,test_exon_idx)-MISALIGN)&&...
            (vertices(1,exon_idx)<=vertices(2,test_exon_idx))

          new_vertex = [vertices(1,test_exon_idx);vertices(2,exon_idx)];
          new_edges = or(edges(exon_idx,:),edges(test_exon_idx,:));
      
          known = 0;
          for ix = 1:size(vertices,2)
            % assert(isequal(find(new_edges),intersect(find(new_edges),find(edges(ix,:)))) == ~any(new_edges & ~edges(ix,:)));
            if isequal(new_vertex,vertices(:,ix)) && ~any(new_edges & ~edges(ix,:))
              known = ix;
              break;
            end
          end
          assert((known~=exon_idx)&&(known~=test_exon_idx));
        
          if ~known
            vertices = [vertices,new_vertex];
            edges = [edges;new_edges];
            edges = [edges,[new_edges';0]];
            to_keep(end+1) = 1;
            to_keep(exon_idx) = 0;
            to_keep(test_exon_idx) = 0;
            changed = 1;
          elseif to_keep(exon_idx) || to_keep(test_exon_idx)
            to_keep(exon_idx) = 0;
            to_keep(test_exon_idx) = 0;
            changed = 1;
          end
        end

      % 1001
      elseif (cur_edge_left&&~cur_edge_right&&~test_edge_left&&test_edge_right)
        %%%   >----- exon -----XM            X>=Y
        %%%          MY--- test_exon ---<
        if (vertices(1,exon_idx)<=vertices(1,test_exon_idx)+MISALIGN)&&...
            (vertices(2,exon_idx)-MISALIGN<=vertices(2,test_exon_idx))&&...
            (vertices(2,exon_idx)>=vertices(1,test_exon_idx))

          new_vertex = [vertices(1,exon_idx);vertices(2,test_exon_idx)];
          new_edges = or(edges(exon_idx,:),edges(test_exon_idx,:));
          known = 0;
          for ix = 1:size(vertices,2)
            % assert(isequal(find(new_edges),intersect(find(new_edges),find(edges(ix,:))))==~any(new_edges & ~edges(ix,:)));
            if isequal(new_vertex,vertices(:,ix)) && ~any(new_edges & ~edges(ix,:))
              known = ix;
              break;
            end
          end
          assert((known~=exon_idx)&&(known~=test_exon_idx));
        
          if ~known
            vertices = [vertices,new_vertex];
            edges = [edges;new_edges];
            edges = [edges,[new_edges';0]];
            to_keep(end+1) = 1;
            to_keep(exon_idx) = 0;
            to_keep(test_exon_idx) = 0;
            changed = 1;
          elseif to_keep(exon_idx) || to_keep(test_exon_idx)
            to_keep(exon_idx) = 0;
            to_keep(test_exon_idx) = 0;
            changed = 1;
          end
        end
  
      % 0111 
      elseif (~cur_edge_left&&cur_edge_right&&test_edge_left&&test_edge_right)
        %%%        ----- exon -----<           
        %%%      >--- test_exon ---<
        if (vertices(2,exon_idx)==vertices(2,test_exon_idx))
        
          %    (vertices(1,exon_idx)+MISALIGN>=vertices(1,test_exon_idx))
          % index to the right
          right_edge = max(vertices(2,exon_idx),vertices(2,test_exon_idx));
          idx = find(vertices(1,:)>right_edge, 1, 'first');
        
          if (~isequal(edges(exon_idx,idx:end),edges(test_exon_idx,idx:end)))
          
            edges(exon_idx,idx:end) = or(edges(exon_idx,idx:end),edges(test_exon_idx,idx:end));
            edges(test_exon_idx,idx:end) = edges(exon_idx,idx:end);
            edges(idx:end,exon_idx) = or(edges(exon_idx,idx:end),edges(test_exon_idx,idx:end));
            edges(idx:end,test_exon_idx) = edges(idx:end,exon_idx);
          
            %to_keep(exon_idx) = 0;
            changed = 1;
          end
          %%%        M---- exon -----<           
          %%%      >--- test_exon ---<
          if (vertices(1,exon_idx)+MISALIGN >= vertices(1,test_exon_idx)) && (to_keep(exon_idx))
            to_keep(exon_idx) = 0;
            changed = 1;
          end
        end
  
      % 1101
      elseif (cur_edge_left&&cur_edge_right&&~test_edge_left&&test_edge_right)
        %%%        >------ exon -----<           
        %%%         --- test_exon ---<
        if (vertices(2,exon_idx)==vertices(2,test_exon_idx))
        
          %    (vertices(1,exon_idx)-MISALIGN<=vertices(1,test_exon_idx))
          % index to the right
          right_edge = max(vertices(2,exon_idx),vertices(2,test_exon_idx));
          idx = find(vertices(1,:)>right_edge,1,'first');
      
          if (~isequal(edges(exon_idx,idx:end),edges(test_exon_idx,idx:end)))
            
            edges(exon_idx,idx:end) = or(edges(exon_idx,idx:end),edges(test_exon_idx,idx:end));
            edges(test_exon_idx,idx:end) = edges(exon_idx,idx:end);
            edges(idx:end,exon_idx) = or(edges(exon_idx,idx:end),edges(test_exon_idx,idx:end));
            edges(idx:end,test_exon_idx) = edges(idx:end,exon_idx);
          
            %to_keep(test_exon_idx) = 0;
            changed = 1;
          end
          %%%    M-------- exon -----<           
          %%%      >--- test_exon ---<
          if (vertices(1,exon_idx)-MISALIGN <= vertices(1,test_exon_idx)) && (to_keep(test_exon_idx))
            to_keep(test_exon_idx) = 0;
            changed = 1;
          end
        end

      % 1011  
      elseif (cur_edge_left&&~cur_edge_right&&test_edge_left&&test_edge_right)
        %%%      >------ exon ---           
        %%%      >--- test_exon ---<
        if (vertices(1,exon_idx)==vertices(1,test_exon_idx))
        
          %(vertices(2,exon_idx)-MISALIGN<=vertices(2,test_exon_idx))
          % index to the left
          left_edge = min(vertices(1,exon_idx),vertices(1,test_exon_idx));
          idx = find(vertices(2,:)<left_edge,1,'last');
        
          if (~isequal(edges(exon_idx,1:idx),edges(test_exon_idx,1:idx)))
          
            edges(exon_idx,1:idx) = or(edges(exon_idx,1:idx),edges(test_exon_idx,1:idx));
            edges(test_exon_idx,1:idx) = edges(exon_idx,1:idx);
            edges(1:idx,exon_idx) = or(edges(exon_idx,1:idx),edges(test_exon_idx,1:idx));
            edges(1:idx,test_exon_idx) = edges(1:idx,exon_idx);
          
            %to_keep(exon_idx) = 0;
            changed = 1;
          end
          %%%     >------ exon ---M           
          %%%     >--- test_exon ----<
          if (vertices(2,exon_idx)-MISALIGN <= vertices(2,test_exon_idx)) && (to_keep(exon_idx))
            to_keep(exon_idx) = 0;
            changed = 1;
          end
        end

      % 1110
      elseif (cur_edge_left&&cur_edge_right&&test_edge_left&&~test_edge_right)
        %%%      >------ exon ------<           
        %%%      >--- test_exon ---
        if (vertices(1,exon_idx)==vertices(1,test_exon_idx))
        
          %    (vertices(2,exon_idx)+MISALIGN>=vertices(2,test_exon_idx))
          % index to the left
          left_edge = min(vertices(1,exon_idx),vertices(1,test_exon_idx));
          idx = find(vertices(2,:)<left_edge,1,'last') ;
      
          if (~isequal(edges(exon_idx,1:idx),edges(test_exon_idx,1:idx)))
            edges(exon_idx,1:idx) = or(edges(exon_idx,1:idx),edges(test_exon_idx,1:idx));
            edges(test_exon_idx,1:idx) = edges(exon_idx,1:idx);
            edges(1:idx,exon_idx) = or(edges(exon_idx,1:idx),edges(test_exon_idx,1:idx));
            edges(1:idx,test_exon_idx) = edges(1:idx,exon_idx);
        
            %to_keep(test_exon_idx) = 0;
            changed = 1;
          end
          %%%      >------ exon -------<           
          %%%      >--- test_exon ---M
          %if (vertices(2,exon_idx) >= vertices(2,test_exon_idx)-MISALIGN) && (to_keep(test_exon_idx))
          if (vertices(2,exon_idx)+MISALIGN >= vertices(2,test_exon_idx)) && (to_keep(test_exon_idx))
            to_keep(test_exon_idx) = 0;
            changed = 1;
          end
        end

      % 1111
      elseif (cur_edge_left&&cur_edge_right&&test_edge_left&&test_edge_right)
        %%%      >------ exon ------<           
        %%%      >--- test_exon ----<
        if (vertices(1,exon_idx)==vertices(1,test_exon_idx))&&...
            (vertices(2,exon_idx)==vertices(2,test_exon_idx))
          if (~isequal(edges(exon_idx,:),edges(test_exon_idx,:)))
            edges(exon_idx,:) = or(edges(exon_idx,:),edges(test_exon_idx,:));
            edges(test_exon_idx,:) = edges(exon_idx,:);
            edges(:,exon_idx) = or(edges(exon_idx,:),edges(test_exon_idx,:));
            edges(:,test_exon_idx) = edges(:,exon_idx);
        
            to_keep(test_exon_idx) = 0;
            changed = 1;
          elseif (to_keep(test_exon_idx))
            to_keep(test_exon_idx) = 0;
            changed = 1;
          end
        %%%      >------ exon ----<   OR   >------ exon ------<          
        %%%      >--- test_exon ----<      >--- test_exon --<
        elseif (vertices(1,exon_idx)==vertices(1,test_exon_idx))&&...
              (vertices(2,exon_idx)~=vertices(2,test_exon_idx))&&...
              (to_keep(exon_idx)||to_keep(test_exon_idx))
          % index to the left
          %left_edge = min(vertices(1,exon_idx),vertices(1,test_exon_idx));
          %idx = find(vertices(2,:)<left_edge, 1, 'last');
          idx = find(vertices(2,:)<vertices(1,exon_idx), 1, 'last');

          if (~isequal(edges(exon_idx,1:idx),edges(test_exon_idx,1:idx)))
            edges(exon_idx,1:idx) = or(edges(exon_idx,1:idx),edges(test_exon_idx,1:idx));
            edges(test_exon_idx,1:idx) = edges(exon_idx,1:idx);
            edges(1:idx,exon_idx) = or(edges(exon_idx,1:idx),edges(test_exon_idx,1:idx));
            edges(1:idx,test_exon_idx) = edges(1:idx,exon_idx);
                
            changed = 1;
          end
        %%%      >------ exon ----<   OR   >------ exon ------<          
        %%%    >--- test_exon ----<          >--- test_exon --<
        elseif (vertices(1,exon_idx)~=vertices(1,test_exon_idx))&&...
            (vertices(2,exon_idx)==vertices(2,test_exon_idx))&&...
            (to_keep(exon_idx)||to_keep(test_exon_idx))
          % index to the right
          %right_edge = max(vertices(2,exon_idx),vertices(2,test_exon_idx));
          %idx = find(vertices(1,:)>right_edge,1,'first');
          idx = find(vertices(1,:)>vertices(2,exon_idx),1,'first');

          if (~isequal(edges(exon_idx,idx:end),edges(test_exon_idx,idx:end)))
        
            edges(exon_idx,idx:end) = or(edges(exon_idx,idx:end),edges(test_exon_idx,idx:end));
            edges(test_exon_idx,idx:end) = edges(exon_idx,idx:end);
            edges(idx:end,exon_idx) = or(edges(exon_idx,idx:end),edges(test_exon_idx,idx:end));
            edges(idx:end,test_exon_idx) = edges(idx:end,exon_idx);
          
            changed = 1;
          end
        end
      else
        fprintf(1,'Unknown case!!\n');
        if keyboard_allowed(), keyboard; end ;
      end  % cases
      %test_exon_idx = test_exon_idx+1;
      internal_idx = internal_idx + 1;
    end
    exon_idx = exon_idx + 1;
  end

  genes(gene_idx).splicegraph = {};
  to_keep = find(to_keep);
  genes(gene_idx).splicegraph = {vertices(:,to_keep),edges(to_keep,to_keep)};
end

fprintf(1,'\n');

