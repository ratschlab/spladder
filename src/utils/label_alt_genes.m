function genes = label_alt_genes(genes, CFG) ;
% genes = label_alt_genes(genes, CFG) ;
% 
% This script labels each gene 'is_alt_spliced' that
% has any transcript position that can be intronic or
% exonic.
%
% This script labels each gene 'is_alt', that has a
% simple alternative transcript start or end. Simple meaning
% that the remaining transcript structure is the same and no 
% other exons overlap to the alternative start or end.

tot_exons = 0;
for ix = 1:length(genes),
  if CFG.verbose && (mod(ix,1000) == 0),
    fprintf(CFG.fd_log, '.');
  end;

  num_exons = size(genes(ix).splicegraph{1}, 2) ;
  tot_exons = tot_exons + num_exons ;

  vertices = genes(ix).splicegraph{1} ;
  edges = genes(ix).splicegraph{2} ;
 
  %%% no edges in graph --> continue
  if isempty(edges),
    genes(ix).is_alt_spliced = 0 ;
    genes(ix).is_alt = 0 ;
    continue ;
  end 
  
  init_alt_idx = [] ;
  term_alt_idx = [] ;
  for i = 1:num_exons,
    %%% handle start terminal exons -> no incoming edges
    if ~any(edges(1:i-1, i)),
      %%% find other exons with same end
      idx = setdiff(find(vertices(2, i) == vertices(2,:)), i) ;
      if ~isempty(idx),
        is_simple = 1 ; 
        for j = idx,
          %%% not simple, if different introns follow --> break
          if ~isequal(find(edges(j, j+1:end)) + j, find(edges(i, i+1:end)) + i)
            is_simple = 0 ;
            break ;
          end ;
          %%% not simple, if exons previous to j overlap the start of i --> break
          idx_prev = find(edges(1:j-1, j)) ;
          for k = idx_prev,
            if vertices(2, k) > vertices(1, i)
              is_simple = 0 ;
              break ;
            end ;
          end ;
        end ;
        if ~is_simple, continue ; end ;
        %%% if simple, report alternative initial exon
        init_alt_idx(end + 1) = i ;
      end ;
    end ;
    %%% handle end terminal exons -> no outgoing edges
    if ~any(edges(i, i + 1:num_exons))
      %%% find other exons with the same start
      idx = setdiff(find(vertices(1, i) == vertices(1, :)), i) ;
      if ~isempty(idx) 
        is_simple = 1 ; 
        for j = idx,
          %%% not simple, if different introns precede --> break
          if ~isequal(find(edges(1:j - 1, j)),find(edges(1:i - 1, i)))
            is_simple = 0 ;
            break ;
          end ;
          %%% not simple, if exons following to j start before i ends --> break
          idx_next = find(edges(j, j+1:end)) + j ;
          for k = idx_next,
            if vertices(1, k) < vertices(2, i)
              is_simple = 0 ;
              break ;
            end ;
          end ;
        end ;
        if ~is_simple, continue ; end ;
        %%% if simple, report alternative terminal exon
        term_alt_idx(end + 1) = i ;
      end ;
    end ;
  end ;

  %%% further only consider exons that are neither init_alt nor term_alt
  take_idx = setdiff(1:num_exons, [init_alt_idx, term_alt_idx]) ;
  vertices = genes(ix).splicegraph{1}(:, take_idx) ;
  edges = genes(ix).splicegraph{2}(take_idx, take_idx) ;

  start = min(vertices(1, :)) ;
  stop = max(vertices(2, :)) ;
  
  exon_loc = zeros(1, stop - start + 1);
  
  %%% map all edges (introns) to genomic coordinates
  for i = 1 : size(edges, 1),
    for j = i + 1 : size(edges, 1),
      if edges(i, j) == 1,
        cur_edge = [vertices(2, i) + 1, vertices(1, j) - 1] - start;
        exon_loc(cur_edge(1):cur_edge(2)) = exon_loc(cur_edge(1):cur_edge(2)) + 1;
      end ;
    end ;
  end ;

  %%% map all vertices (exons) to genomic coordinates 
  for i = 1:size(vertices, 2),
    cur_vertex = vertices(:, i) - [start; start] + 1;
    exon_loc(cur_vertex(1):cur_vertex(2)) = exon_loc(cur_vertex(1):cur_vertex(2)) + 1;
  end ;

  %%% if at any position more than one exon or intron -> is_alt__spliced
  if max(exon_loc) > 1,
    genes(ix).is_alt_spliced = 1 ;
    genes(ix).is_alt = 1 ;
  else
    genes(ix).is_alt_spliced = 0 ;
    %%% if not alt_spliced but term_alt or init_alt --> is_alt
    if isempty([init_alt_idx, term_alt_idx])
      genes(ix).is_alt = 0 ;
    else
      genes(ix).is_alt = 1 ;
    end ;
  end ;
end ;

if CFG.verbose,
    fprintf(CFG.fd_log,'\n\nTotal genes:\t\t\t\t\t\t\t%d\n',...
            length(genes));
    fprintf(CFG.fd_log,'Total exons:\t\t\t\t\t\t\t%d\n',...
            tot_exons);
    fprintf(CFG.fd_log,'Total genes with alternative isoforms:\t\t\t\t%d\n',...
            sum([genes(:).is_alt]));
    fprintf(CFG.fd_log,'Total genes alternatively spliced:\t\t\t\t%d\n',...
            sum([genes(:).is_alt_spliced]));
    fprintf(CFG.fd_log,'Total constitutively spliced:\t\t\t\t\t%d\n',...
            sum(~[genes(:).is_alt_spliced]));
end;
