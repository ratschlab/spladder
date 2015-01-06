function [idx_multiple_skips, exon_multiple_skips, id_multiple_skips] = detect_multipleskips(genes, idx_alt) ;

id = 0 ;
idx_multiple_skips = zeros(1,0);
id_multiple_skips = zeros(1,0);
exon_multiple_skips = zeros(1,0) ;
for ix = idx_alt,
  if (mod(ix, 50) == 0),
    fprintf(1, '.');
  end;
  num_exons = size(genes(ix).splicegraph{1}, 2) ;
  edges = genes(ix).splicegraph{2} ;
  labels = repmat([1:num_exons]',1,num_exons) ;
  
  %adjecency matrix: upper half only
  A = sparse(zeros(num_exons)) ;
  for i = 1:num_exons-1,
    for j = i+1:num_exons,
      A(i,j) = edges(i,j) ;
    end;
  end;
  
  %possible starting and ending exons of a multiple exon skip
  Pairs = sparse(zeros(length(A))) ;
  Ai = (A * A) * A ; %paths of length 3
  while any(Ai(:)>0)
    [x_coords, y_coords] = find(A>0 & Ai>0) ; %multiple skip
    Pairs(x_coords, y_coords) = 1 ;
    Ai = Ai * A ; %paths of length ..+1
  end
  
  [edge1, edge2] = find(Pairs) ;
  
  if length(find(Pairs==1)) > 10000
    fprintf(1,'Warning: not processing gene %d, because there are more than 10000 potential hits.\n',ix);
    continue ;
  end
  
  for cnt = 1:size(edge1,1)
    exon_idx_first = edge1(cnt) ;
    exon_idx_last = edge2(cnt) ;

    if (edges(exon_idx_first,exon_idx_last))
      
      % find all pairs shortest path
      exist_path = double(triu(edges));
      exist_path(exon_idx_first,exon_idx_last)=0;
      exist_path(exist_path==0)=inf ;
      
      for i=1:num_exons; exist_path(i,i)=0 ; end ;
      
      long_exist_path = double(-triu(edges));
      long_exist_path(exon_idx_first,exon_idx_last)=0;
      long_exist_path(long_exist_path==0)=inf ;
      for i=1:num_exons; long_exist_path(i,i)=0 ; end ;
      
      
      path = ~isinf(exist_path).*labels;
      long_path = ~isinf(long_exist_path).*labels;
      
      for k=1:num_exons,
        for i=1:num_exons,
          idx = find((exist_path(i,k) + exist_path(k,:)) < exist_path(i,:));
          exist_path(i,idx) = exist_path(i,k) + exist_path(k,idx) ;
          path(i,idx) = path(k,idx) ;
          
          idx = find((long_exist_path(i,k) + long_exist_path(k,:)) < long_exist_path(i,:));
          long_exist_path(i,idx) = long_exist_path(i,k) + long_exist_path(k,idx);
          long_path(i,idx) = long_path(k,idx);
        end ;
      end ;
      
      temp_ix = ~isinf(long_exist_path);
      long_exist_path(temp_ix) = -long_exist_path(temp_ix);
      
      if (exist_path(exon_idx_first,exon_idx_last) > 2) && ~isinf(exist_path(exon_idx_first,exon_idx_last)),
        backtrace = path(exon_idx_first,exon_idx_last);
        while (backtrace(end) > exon_idx_first)
          backtrace = [backtrace,path(exon_idx_first,backtrace(end))];
        end
        backtrace = backtrace(1:end-1);
        backtrace = backtrace(length(backtrace):-1:1);
        idx_multiple_skips = [idx_multiple_skips, repmat(ix,1,length(backtrace)+2)];
        exon_multiple_skips = [exon_multiple_skips, exon_idx_first backtrace exon_idx_last];
        id = id + 1 ;
        id_multiple_skips = [id_multiple_skips,id*ones(1,length(backtrace)+2)];
      elseif (long_exist_path(exon_idx_first,exon_idx_last) > 2) && ~isinf(long_exist_path(exon_idx_first,exon_idx_last)),
        backtrace = long_path(exon_idx_first,exon_idx_last);
        while (backtrace(end) > exon_idx_first)
          backtrace = [backtrace, long_path(exon_idx_first, backtrace(end))];
        end
        backtrace = backtrace(1:end-1);
        backtrace = backtrace(length(backtrace):-1:1);
        idx_multiple_skips = [idx_multiple_skips,repmat(ix, 1, length(backtrace)+2)];
        exon_multiple_skips = [exon_multiple_skips, exon_idx_first backtrace exon_idx_last];
        id = id + 1 ;
        id_multiple_skips = [id_multiple_skips, id * ones(1, length(backtrace)+2)];
      end
      %if exist_path(exon_idx_first,exon_idx_last) < ...
      %      long_exist_path(exon_idx_first,exon_idx_last)
      %  disp(ix);
      %end
      %if exist_path(exon_idx_first,exon_idx_last) > ...
      %      long_exist_path(exon_idx_first,exon_idx_last)
      %  keyboard
      %end
      
    end
  end
end
fprintf(1,'\nNumber of multiple exon skips:\t\t\t\t\t%d\n',...
        length(unique(id_multiple_skips)));



%for ix=1:length(idx_multiple_skips);
%  disp(idx_multiple_skips(ix));
%  viewsplicegraph_compare(genes(idx_multiple_skips(ix)),genes(idx_multiple_skips(ix)).splicegraph{1}(:,exon_multiple_skips(ix)));
%  pause;
%end;
