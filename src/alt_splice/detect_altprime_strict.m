function [idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] ...
    = detect_altprime_strict(genes, idx_alt);
% function [idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] ...
%    = detect_altprime(genes, idx_alt);
%
% detect the alternative 5 and 3 prime ends of the intron. Note that 5 prime refers to the left
% and 3 prime to the right for a positive strand
%
% The difference from detect_altprime is that we filter out intron retentions.
% 
% Due to the limitation of intersect, if there are more than one intron retentions
% in the gene, and one of them causes the alternative detection, this function
% may fail to catch that case.

MIN_OVERLAP = 11; % two exons that are alternative should overlap by at least the BLAT seed

idx_alt_5prime = [];
idx_alt_3prime = [];
exon_alt_5prime = {};
exon_alt_3prime = {};


for ix=idx_alt
  if (mod(ix,50)==0)
    fprintf(1,'.');
  end
  num_exons = size(genes(ix).splicegraph{1},2) ;
  vertices = genes(ix).splicegraph{1} ;
  edges = genes(ix).splicegraph{2} ;
  if isfield(genes(ix),'strands')
    strand = genes(ix).strands(1);
  else
    strand = genes(ix).strand;
  end
  % Find alternative sites on the right of the intron,
  % same site on the left of the intron.
  % Store the splice site at the left of the intron in 'starts'.
  starts=[] ;
  for exon_idx = 1:num_exons-2
    if sum(find(starts==vertices(2,exon_idx)))==0,
      starts=[starts, vertices(2,exon_idx)] ;
      nr_exons=sum(edges(exon_idx,exon_idx+1:num_exons)) ;
      rightsites=[] ;
      rightidx=[] ;
      if nr_exons>=2,
	which_exons = find([zeros(1,exon_idx),edges(exon_idx,exon_idx+1:num_exons)]);
        exons=[vertices(:, which_exons)] ;
        for i=1:nr_exons-1
          for j=i+1:nr_exons
	    % check that the left splice site of the exons are different
	    % make sure exons overlap - either:
	    % - left splice site of exon(i) is in exon(j)
	    % - left splice site of exon(j) is in exon(i)
	    % note that the 'overlap' relationship is not transitive
	    if ( (exons(1,i) ~= exons(1,j)) && ...
	    	 (( (exons(1,i) >= exons(1,j)) && (exons(1,i) <= exons(2,j)) ) ...
	    	  || ( (exons(1,j) >= exons(1,i)) && (exons(1,j) <= exons(2,i)) ) ) && ...
	    	 (min(exons(2,i),exons(2,j))-max(exons(1,i),exons(1,j)) >= MIN_OVERLAP ) )
	      assert(~((exons(1,i) == exons(1,j)) && (exons(2,i) == exons(2,j))));
	      assert(exons(1,i) ~= exons(1,j));
	    % check that the right (exon) sites are different but left splice is same.
	    %if ( (exons(2,i) ~= exons(2,j)) && (exons(1,i) == exons(1,j)) )
	      % add new events to the list
	      if sum(find(rightsites(1:end)==exons(1,i)))==0
                rightsites=[rightsites, exons(1,i)] ;
                rightidx=[rightidx, which_exons(i)] ;
              end
              if sum(find(rightsites(1:end)==exons(1,j)))==0
                rightsites=[rightsites, exons(1,j)] ;
                rightidx=[rightidx, which_exons(j)] ;
              end
            end  
          end
        end     
      end
     
      % construct the output
      if length(rightsites)>=2,
        if strand=='+',
          exon_alt_3prime(end+1).fiveprimesite=exon_idx ;
          exon_alt_3prime(end).threeprimesites=rightidx ;
          idx_alt_3prime = [idx_alt_3prime, ix] ;
        end
        if strand=='-',
          exon_alt_5prime(end+1).threeprimesite=exon_idx ;
          exon_alt_5prime(end).fiveprimesites=rightidx ;
          idx_alt_5prime = [idx_alt_5prime, ix] ;
        end
      end  % construct output
    end
  end % for exon_idx

  
  % Find alternative sites on the left of the intron,
  % same site on the right of the intron.
  % Store the splice site at the right of the intron in 'starts'.
  starts=[] ;
  for exon_idx = 3:num_exons
    if sum(find(starts==vertices(1,exon_idx)))==0,
      starts=[starts, vertices(1,exon_idx)] ;
      nr_exons=sum(edges(1:exon_idx-1,exon_idx)) ;
      leftsites=[] ;
      leftidx=[] ;
      if nr_exons>=2,
	which_exons = find([edges(1:exon_idx-1,exon_idx)', zeros(1,num_exons-exon_idx+1)]);
        exons=[vertices(:, which_exons)] ;
        for i=1:nr_exons-1
          for j=i+1:nr_exons
	    % check that the 5prime sites are different
	    % make sure exons overlap - either:
	    % - right splice site of exon(i) is in exon(j)
	    % - right splice site of exon(j) is in exon(i)
	    % note that the 'overlap' relationship is not transitive
	    if ((exons(2,i) ~= exons(2,j)) &&...
	    	( ( (exons(2,i) <= exons(2,j)) && (exons(2,i) >= exons(1,j)) )...
	    	|| ( (exons(2,j) <= exons(2,i)) && (exons(2,j) >= exons(1,i)) ) ) && ...
	    	(min(exons(2,i),exons(2,j))-max(exons(1,i),exons(1,j)) >= MIN_OVERLAP ) )
	      assert(~((exons(1,i) == exons(1,j)) && (exons(2,i) == exons(2,j))));
	      assert(exons(2,i) ~= exons(2,j));
	      
	    % check that the right (exon) sites are different but left splice is same.
	    %if ( (exons(2,i) ~= exons(2,j)) && (exons(1,i) == exons(1,j)) )
	      % add new events to the list
              if sum(find(leftsites(1:end)==exons(2,i)))==0
                leftsites=[leftsites, exons(2,i)] ;
                leftidx=[leftidx, which_exons(i)] ;
              end
              if sum(find(leftsites(1:end)==exons(2,j)))==0
                leftsites=[leftsites, exons(2,j)] ;
                leftidx=[leftidx, which_exons(j)] ;
              end
            end
          end
        end     
      end

      % construct the output
      if length(leftsites)>=2,
        if strand=='+'
          exon_alt_5prime(end+1).threeprimesite=exon_idx ;
          exon_alt_5prime(end).fiveprimesites=leftidx ;
          idx_alt_5prime = [idx_alt_5prime, ix] ;
        end
        if strand =='-'
          exon_alt_3prime(end+1).fiveprimesite=exon_idx ;
          exon_alt_3prime(end).threeprimesites=leftidx ;
          idx_alt_3prime = [idx_alt_3prime, ix] ;
        end
      end % construct output
    end
  end % for exon_idx
  


end





fprintf(1,'\n\nNumber of alternative 5 prime sites:\t\t\t\t%d\n',...
	length(idx_alt_5prime));

fprintf(1,'Number of alternative 3 prime sites:\t\t\t\t%d\n',...
	length(idx_alt_3prime));




% filter out the alternative sites which overlap with intron retention sites.
[idx_intron_reten,intron_intron_reten] = detect_intronreten(genes,idx_alt);
[idx_exon_skips, exon_exon_skips] = detect_exonskips(genes, idx_alt) ;

[idx_both,ixa,ixb] = intersect(idx_alt_5prime,idx_intron_reten);

overlap = [];
for ix = 1:length(idx_both)
  fiveexons = exon_alt_5prime(ixa(ix)).fiveprimesites;
  if length(intersect(fiveexons,intron_intron_reten(:,ixb(ix))'))>1
    idx_skip = find(idx_exon_skips==idx_both(ix));
    if ~isempty(idx_skip)
      skipped_exons = exon_exon_skips(:,idx_skip);
      for ixe = 1:size(skipped_exons,2)
	for ixf = 1:length(fiveexons)
	  if all(sort(skipped_exons([1 3],ixe)') == sort([exon_alt_5prime(ixa(ix)).threeprimesite fiveexons(ixf)]))
	    overlap(end+1) = ix;
	  end
	end
      end
    end
  end
end

idx_alt_5prime(ixa(overlap)) = [];
exon_alt_5prime(ixa(overlap)) = [];

fprintf(1,'Filter out 5 prime sites due to intron retention and exon skips: %d\n',...
	length(overlap));



[idx_both,ixa,ixb] = intersect(idx_alt_3prime,idx_intron_reten);
overlap = [];
for ix = 1:length(idx_both)
  threeexons = exon_alt_3prime(ixa(ix)).threeprimesites;
  if length(intersect(threeexons,intron_intron_reten(:,ixb(ix))'))>1
    idx_skip = find(idx_exon_skips==idx_both(ix));
    if ~isempty(idx_skip)
      skipped_exons = exon_exon_skips(:,idx_skip);
      for ixe = 1:size(skipped_exons,2)
	for ixf = 1:length(threeexons)
	  if all(sort(skipped_exons([1 3],ixe)') == sort([exon_alt_3prime(ixa(ix)).fiveprimesite threeexons(ixf)]))
	    overlap(end+1) = ix;
	  end
	end
      end
    end
  end
end

idx_alt_3prime(ixa(overlap)) = [];
exon_alt_3prime(ixa(overlap)) = [];


fprintf(1,'Filter out 3 prime sites due to intron retention and exon skips: %d\n',...
	length(overlap));

