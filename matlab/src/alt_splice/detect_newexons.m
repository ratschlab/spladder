function [genes, new_exon_count]=detect_newexons(genes, idx_alt) ;

new_exon_count = 0 ;
for ix=idx_alt
  if (mod(ix,50)==0)
    fprintf(1,'.');
  end
  genes(ix).splicegraph{3}=[] ; % for each exon: 1 if new, 0 if not
  show_splicegraph=0 ;
  for k=1:size(genes(ix).splicegraph{1},2)
    genes(ix).splicegraph{3}(end+1)=1 ; % exon is new?
    start_end=genes(ix).splicegraph{1}(:,k) ;
    exon_start=start_end(1) ;
    exon_end=start_end(2) ;
    for j=1:size(genes(ix).exons,2)
      for l=1:size(genes(ix).exons{j},1)
        exon=genes(ix).exons{j}(l,:) ;
        if exon(1)==exon_start & exon(2)==exon_end
          %fprintf('found exon: %i %i\n', exon(1), exon(2)) ;
          genes(ix).splicegraph{3}(end)=0 ; % exon not new
        end 
      end 
    end
    if genes(ix).splicegraph{3}(end)==1,
      show_splicegraph=1 ;
      new_exon_count = new_exon_count+1 ;
      %ix
    end
  end
  if 0 & show_splicegraph==1 
    gene=genes(ix) ;
    viewsplicegraph_newexons(gene) ;
    genes(ix).splicegraph{3}
    pause
  end
end

fprintf(1,'\n\nNumber of new exons:\t\t\t\t\t\t%d\n',...
	new_exon_count);
