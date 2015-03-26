function [idx_intron_reten,intron_intron_reten] = detect_intronreten(genes,idx_alt) ;

idx_intron_reten = [];
intron_intron_reten = [];
for ix=idx_alt
  %viewsplicegraph(genes(ix)) ;
  if (mod(ix,50)==0)
    fprintf(1,'.');
  end
  num_exons = size(genes(ix).splicegraph{1},2);
  vertices = genes(ix).splicegraph{1};
  edges = genes(ix).splicegraph{2};
  introns=[] ;
  for exon_idx = 1:num_exons-1 %start of intron
    %exon_idx
    %edges(exon_idx,exon_idx+1:length(edges))
    idx=find(edges(exon_idx,exon_idx+1:length(edges))==1) ;
    if isempty(idx), continue; end
    idx=idx+exon_idx ;
    for exon_idx2=idx, %end of intron
      is_intron_reten=0 ;
      for exon_idx1=1:num_exons %exon
	% check that the exon covers the intron
        if (vertices(2,exon_idx) >= vertices(1,exon_idx1)) && (vertices(1,exon_idx2) <= vertices(2,exon_idx1))
	
	% check that the exon covers the intron, but doesn't extend beyond
        %if (vertices(2,exon_idx) > vertices(1,exon_idx1)) && (vertices(1,exon_idx2) < vertices(2,exon_idx1))...
	%      && (vertices(1,exon_idx) <= vertices(1,exon_idx1)) && (vertices(2,exon_idx2) >= vertices(2,exon_idx1))
          is_intron_reten=1 ;
          long_exon=exon_idx1 ;
          for len=1:size(introns,2),
            if (vertices(2,exon_idx)==introns(1,len)) && (vertices(1,exon_idx2)==introns(2,len))
              is_intron_reten=0 ;
            end  
          end
        end
      end
      if is_intron_reten==1,
        idx_intron_reten = [idx_intron_reten, ix];
        intron_intron_reten = [intron_intron_reten, [exon_idx ; exon_idx2 ...
                   ; long_exon]] ;
        introns=[introns, [vertices(2,exon_idx) ; vertices(1,exon_idx2)]] ;
      end
    end
    %pause ;
  end
  
  %pause ;
end
fprintf(1,'\n\nNumber of intron retentions:\t\t\t\t\t%d\n', ...
	length(idx_intron_reten));

