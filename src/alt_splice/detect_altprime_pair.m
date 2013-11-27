function [idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] ...
    = detect_altprime_pair(genes, idx_alt);
% function [idx_alt_5prime,exon_alt_5prime, idx_alt_3prime,exon_alt_3prime] ...
%    = detect_altprime_pair(genes, idx_alt);
%
% Detect the alternative 5 and 3 prime ends of the intron. Note that 5 prime refers to the left
% and 3 prime to the right for a positive strand
% The difference between this function and detect_altprime.m is that the other function collects
% all the alternative sites corresponding to an acceptor or donor.
% note that the 'overlap' relationship is not transitive

MIN_OVERLAP = 11; % two exons that are alternative should overlap by at least the BLAT seed

idx_alt_5prime = [];
idx_alt_3prime = [];
exon_alt_5prime = {};
exon_alt_3prime = {};

if ~isfield(genes, 'strands'),
  for i=1:length(genes),
    genes(i).strands = repmat(genes(i).strand, 1, length(genes(i).transcripts)) ;
  end ;
end ;

for ix=idx_alt
  if (mod(ix,50)==0)
    fprintf(1,'.');
  end
  num_exons = size(genes(ix).splicegraph{1},2) ;
  vertices = genes(ix).splicegraph{1} ;
  edges = genes(ix).splicegraph{2} ;
  
  % Find alternative sites on the right of the intron,
  % same site on the left of the intron.
  for exon_idx = 1:num_exons-2
    nr_exons=sum(edges(exon_idx,exon_idx+1:num_exons));
    if nr_exons>=2,
      which_exons = find([zeros(1,exon_idx),edges(exon_idx,exon_idx+1:num_exons)]);
      exons=[vertices(:, which_exons)] ;
      for i=1:nr_exons-1
	for j=i+1:nr_exons
	  % check that the left splice site of the exons are different
	  % make sure exons overlap - either:
	  % - left splice site of exon(i) is in exon(j)
	  % - left splice site of exon(j) is in exon(i)
	  
	  if ((exons(1,i) ~= exons(1,j)) && ...
	       (( (exons(1,i) >= exons(1,j)) && (exons(1,i) <= exons(2,j)) ) ...
		|| ( (exons(1,j) >= exons(1,i)) && (exons(1,j) <= exons(2,i)) )) && ...
	       (min(exons(2,i),exons(2,j))-max(exons(1,i),exons(1,j)) >= MIN_OVERLAP ) )
	    assert(~((exons(1,i) == exons(1,j)) && (exons(2,i) == exons(2,j))));
	    assert(exons(1,i) ~= exons(1,j));
	    % construct the output
	    if (genes(ix).strands(1)=='+'),
	      exon_alt_3prime(end+1).fiveprimesite=exon_idx ;
	      if exons(1,i) > exons(1,j)
		exon_alt_3prime(end).threeprimesites=which_exons([i j]);
	      else
		exon_alt_3prime(end).threeprimesites=which_exons([j i]);
	      end
	      idx_alt_3prime(end+1) = ix;
	    end
	    if (genes(ix).strands(1)=='-'),
	      exon_alt_5prime(end+1).threeprimesite=exon_idx ;
	      if exons(1,i) > exons(1,j)
		exon_alt_5prime(end).fiveprimesites=which_exons([i j]);
	      else
		exon_alt_5prime(end).fiveprimesites=which_exons([j i]);
	      end
	      idx_alt_5prime(end+1) = ix;
	    end
	  end % if altsplice
	  
	end % for j=i+1:nr_exons
      end % i=1:nr_exons-1
    end
  end % for exon_idx

  
  % Find alternative sites on the left of the intron,
  % same site on the right of the intron.
  for exon_idx = 3:num_exons
    nr_exons=sum(edges(1:exon_idx-1,exon_idx)) ;
    if nr_exons>=2,
      which_exons = find([edges(1:exon_idx-1,exon_idx)', zeros(1,num_exons-exon_idx+1)]);
      exons=[vertices(:, which_exons)] ;
      for i=1:nr_exons-1
	for j=i+1:nr_exons
	  % check that the 5prime sites are different
	  % make sure exons overlap - either:
	  % - right splice site of exon(i) is in exon(j)
	  % - right splice site of exon(j) is in exon(i)
	  if ((exons(2,i) ~= exons(2,j)) &&...
	      (( (exons(2,i) <= exons(2,j)) && (exons(2,i) >= exons(1,j)) )...
	      || ( (exons(2,j) <= exons(2,i)) && (exons(2,j) >= exons(1,i)) )) && ...
	      (min(exons(2,i),exons(2,j))-max(exons(1,i),exons(1,j)) >= MIN_OVERLAP ) )
	    assert(~((exons(1,i) == exons(1,j)) && (exons(2,i) == exons(2,j))));
	    assert(exons(2,i) ~= exons(2,j));

	    % construct the output
	    if (genes(ix).strands(1)=='+'),
	      exon_alt_5prime(end+1).threeprimesite=exon_idx ;
	      if exons(2,i) < exons(2,j)
		exon_alt_5prime(end).fiveprimesites=which_exons([i j]);
	      else
		exon_alt_5prime(end).fiveprimesites=which_exons([j i]);
	      end
	      idx_alt_5prime(end+1) = ix;
	    end
	    if (genes(ix).strands(1)=='-'),
	      exon_alt_3prime(end+1).fiveprimesite=exon_idx ;
	      if exons(2,i) < exons(2,j)
		exon_alt_3prime(end).threeprimesites=which_exons([i j]);
	      else
		exon_alt_3prime(end).threeprimesites=which_exons([j i]);
	      end
	      idx_alt_3prime(end+1) = ix;
	    end
	  end % if altsplice
	  
	end % for j=i+1:nr_exons
      end % for i=1:nr_exons-1
    end

  end % construct output
end % for exon_idx
  



fprintf(1,'\n\nNumber of alternative 5 prime sites:\t\t\t\t%d\n',...
	length(idx_alt_5prime));

fprintf(1,'Number of alternative 3 prime sites:\t\t\t\t%d\n',...
	length(idx_alt_3prime));





