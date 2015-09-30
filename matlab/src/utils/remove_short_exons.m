function [genes] = remove_short_exons(genes, CFG)
% [genes] = remove_short_exons(genes, terminal_short_extend, terminal_short_len, short_exon_skipped, short_exon_removed)

short_exon_removed = 0 ;
short_exon_skipped = 0 ;

rm_idx = [] ;
for i = 1:length(genes)
  if CFG.verbose && mod(i,1000) == 0, 
    fprintf(CFG.log_fd, '%i\n', i);
  end ;

  s = (genes(i).strand == '-') + 1 ;
  %%% remove genes with no exons
  if isempty(genes(i).splicegraph{2}), rm_idx(end + 1) = i ; continue ; end ;

  if genes(i).splicegraph{1}(2, 1) - genes(i).splicegraph{1}(1, 1) < CFG.remove_exons.terminal_short_len,
    genes(i).splicegraph{1}(1, 1) = genes(i).splicegraph{1}(2, 1) - CFG.remove_exons.terminal_short_extend ;
    genes(i).start = min(genes(i).start, genes(i).splicegraph{1}(1, 1)) ;
  end ;

  if genes(i).splicegraph{1}(2, end) - genes(i).splicegraph{1}(1, end) < CFG.remove_exons.terminal_short_len,
    genes(i).splicegraph{1}(2, end) = genes(i).splicegraph{1}(2, end) + CFG.remove_exons.terminal_short_extend ;
    genes(i).stop = max(genes(i).stop, genes(i).splicegraph{1}(2, end)) ;
  end ;

  % check for very short exons and insert an edge that allows skipping them
  exons_remove_idx = [] ;
  j = 2 ;
  while (j <= size(genes(i).splicegraph{1}, 2) - 1),
    genes(i).splicegraph{1}(2, j) - genes(i).splicegraph{1}(1, j)
    if genes(i).splicegraph{1}(2, j) - genes(i).splicegraph{1}(1, j) < CFG.remove_exons.min_exon_len,
      foundp = 0 ;
      jp = j + 1 ;
      for jp = j + 1:size(genes(i).splicegraph{1}, 2)
        if genes(i).splicegraph{1}(2, jp) - genes(i).splicegraph{1}(1, jp) >= CFG.remove_exons.min_exon_len_remove && genes(i).splicegraph{2}(j, jp),
          foundp = 1 ;
          break ;
        end ;
      end ;
      foundn = 0 ;
      for jn = j - 1:-1:1,
        if genes(i).splicegraph{1}(2, jn) - genes(i).splicegraph{1}(1, jn) >= CFG.remove_exons.min_exon_len_remove && genes(i).splicegraph{2}(jn, j),
          foundn = 1 ;
          break ;
        end ;
      end ;
      if foundp == 1 && foundn == 1,
        genes(i).splicegraph{2}(jn, jp) = 1 ;
        genes(i).splicegraph{2}(jp, jn) = 1 ;

        if genes(i).splicegraph{1}(2, j) - genes(i).splicegraph{1}(1, j) < CFG.remove_exons.min_exon_len_remove,
          short_exon_removed = short_exon_removed + 1 ;
          exons_remove_idx(end + 1) = j ;
        else
          short_exon_skipped = short_exon_skipped + 1 ;
        end ;
      end ;
    end ;
    j = j + 1 ;
  end ;

  idx = setdiff(1:size(genes(i).splicegraph{1}, 2), exons_remove_idx) ;
  genes(i).splicegraph{1} = genes(i).splicegraph{1}(:, idx) ;
  genes(i).splicegraph{2} = genes(i).splicegraph{2}(idx, idx) ;
  genes(i).splicegraph{3} = genes(i).splicegraph{3}(:, idx) ;

end ;
rm_idx
genes(rm_idx) = [] ;

short_exon_removed
short_exon_skipped
