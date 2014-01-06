function [genes, inserted] = insert_intron_edges(genes, CFG)
% genes = insert_intron_edges(genes, CFG);

if ~isfield(CFG, 'debug'),
    CFG.debug = 0;
end;

print_intermediates = 0;

strands = '+-';
P=[] ; both_missing = [0 0] ; one_missing=[0 0] ; multi=0 ; next=0; prev=0 ;
exon_vicinity_cnt1 = [0 0] ;
exon_vicinity_cnt2 = [0 0] ;
merge_idx = [] ;
intron_tol = 1 ;

inserted.intron_in_exon = 0 ;
inserted.alt_53_prime = 0 ;
inserted.exon_skip = 0 ;
inserted.gene_merge = 0 ;
inserted.new_terminal_exon = 0 ;

num_unused_introns = zeros(1,length(genes)) ;

for i = 1:length(genes)
	if CFG.verbose && mod(i,1000)==0, 
        fprintf(CFG.fd_log, '%i of %i genes\n', i, length(genes)); 
    end ;

	s = find(genes(i).strand == '+-');

	if isempty(genes(i).introns{s})
		continue;
	end

	unused_introns = [] ;
    if CFG.debug,
        fprintf(CFG.fd_log, 'processing gene %i; with %i introns; time since last tic\n', i, length(genes(i).introns{s}));
        toc
        tic
    end;

	for j = 1:size(genes(i).introns{s}, 2),
		intron_used = 0 ;

        if ( j > 1 && size(genes(i).splicegraph{1},2) > 1) 
            genes(i) = uniquify_splicegraph(genes(i));
        end;

		%%% find exons within same gene whose end coincides with intron start
		idx1 = find(abs(genes(i).splicegraph{1}(2,:) - genes(i).introns{s}(1,j) + 1) <= intron_tol) ;
		%%% find exons within same gene whose start coincides with intron end
		idx2 = find(abs(genes(i).splicegraph{1}(1,:) - genes(i).introns{s}(2,j) - 1) <= intron_tol) ;

		%%% intron boundaries do not coincide with any exon boundaries
		if isempty(idx1) && isempty(idx2),
			both_missing(s) = both_missing(s) + 1 ;

			if CFG.intron_edges.insert_intron_retention,
				%%% find all exons that completely include added introns 
				idx1__ = find(genes(i).introns{s}(1,j) > genes(i).splicegraph{1}(1,:) &	genes(i).introns{s}(2,j) < genes(i).splicegraph{1}(2,:)) ;
				for idx1_ = idx1__,

					genes(i).splicegraph{1}(:,end+1) = genes(i).splicegraph{1}(:, idx1_) ;
					genes(i).splicegraph{1}(2,end) = genes(i).introns{s}(1,j)-1 ;
							
					genes(i).splicegraph{1}(:,end+1) = genes(i).splicegraph{1}(:, idx1_) ;
					genes(i).splicegraph{1}(1,end) = genes(i).introns{s}(2,j)+1 ;
							
					genes(i).splicegraph{2}(end+1,end+1)=0 ;
					adj_matrix = triu(genes(i).splicegraph{2}) ;
					genes(i).splicegraph{2}(:,end) = adj_matrix(:,idx1_) ; % incoming edges of idx1_
					genes(i).splicegraph{2}(end,:) = adj_matrix(:,idx1_)' ;
							
					genes(i).splicegraph{2}(end+1,end+1)=0 ;
					adj_matrix = triu(genes(i).splicegraph{2}) ;
					genes(i).splicegraph{2}(:,end) = adj_matrix(idx1_,:)' ; % outgoing edges of idx1_
					genes(i).splicegraph{2}(end,:) = adj_matrix(idx1_,:) ;
							
					genes(i).splicegraph{2}(end-1,end) = 1 ;
					genes(i).splicegraph{2}(end,end-1) = 1 ;
							
					genes(i).splicegraph{3}(:,end+1) = genes(i).splicegraph{3}(:,idx1_) ; 
					genes(i).splicegraph{3}(2,end) = 0 ; % cannot be an end
							
					genes(i).splicegraph{3}(:,end+1) = genes(i).splicegraph{3}(:,idx1_) ; 
					genes(i).splicegraph{3}(1,end) = 0 ; % cannot be a start
							
					inserted.intron_in_exon = inserted.intron_in_exon + 1 ;
					assert(all(genes(i).splicegraph{1}(1,:) <= genes(i).splicegraph{1}(2,:)));

                    if CFG.debug,
                        fprintf(CFG.fd_log, '%s\tintron_retention_exon\t%c\t%i\t%i\t%i\t%i\n', genes(i).chr, genes(i).strand, genes(i).splicegraph{1}(1,end-1), genes(i).splicegraph{1}(2,end-1), genes(i).splicegraph{1}(1,end), genes(i).splicegraph{1}(2,end)) ;
                    end;
					intron_used = 1 ;
				end ;
			end ;
			if ~intron_used, unused_introns(end + 1) = j ; end ;
			continue ; % with next intron
		end ;

		% did not find exons in same gene sharing boundaries with intron start
		% find first end in previous gene on same strand
		if isempty(idx1) && i>1 && genes(i-1).chr_num==genes(i).chr_num && genes(i-1).strand==genes(i).strand,
			%%% find all exon ends in previuos gene that coincide with intron start j
			idx1_ = find(abs(genes(i-1).splicegraph{1}(2,:) - genes(i).introns{s}(1,j) + 1) <= intron_tol) ;
			if ~isempty(idx1_),
				prev = prev + 1 ;
				% mark the two genes for merging
				if CFG.intron_edges.gene_merges, 
					merge_idx(end+1,:) = [i-1 i] ;
					intron_used = 1 ;
				end ;
				if ~intron_used, unused_introns(end+1) = j ; end ;
				continue ; % with next intron
			end ;
		end ;

		% did not find exons in same gene sharing boundaries with intron end
		% find second end in next gene on same strand
		if isempty(idx2) && i<length(genes) && genes(i+1).chr_num == genes(i).chr_num && genes(i+1).strand == genes(i).strand,
			%%% find all exon starts in following gene that coincide with intron end j
			idx2_ = find(abs(genes(i+1).splicegraph{1}(1,:) - genes(i).introns{s}(2,j) - 1) <= intron_tol) ;
			if ~isempty(idx2_),
				next = next + 1 ;
				% mark the two genes for merging
				if CFG.intron_edges.gene_merges, 
					merge_idx(end+1,:) = [i i+1] ;
					intron_used = 1 ;
				end ;
				if ~intron_used, unused_introns(end+1) = j ; end ;
				continue ;
			end ;
		end ;

		% did not find exons in same gene sharing boundaries with intron start
		% check whether the intron starts in the vicinity of an exon
		if isempty(idx1) 
			%%% find all exons that overlap intron-start j +/- CFG.intron_edges.vicinity_region
			idx1__ = find(genes(i).splicegraph{1}(1,:) - CFG.intron_edges.vicinity_region <= genes(i).introns{s}(1,j) & ...
									genes(i).splicegraph{1}(2,:) + CFG.intron_edges.vicinity_region >= genes(i).introns{s}(1,j)) ;

            %%% check, if we can find an exon after the current intron and there is continuous coverage between intron end and exon
            if isempty(idx1__),
                idx1__ = find(genes(i).splicegraph{1}(1, :) > genes(i).introns{s}(2, j), 1, 'first');
                if ~isempty(idx1__),
                    gg = genes(i);
                    gg.strand = strands(s);
                    rm_strands = ~isfield(gg, 'strands');
                    gg.strands = strands(s);
                    gg.start = genes(i).introns{s}(2, j) + 1; %%% start of presumable exon
                    gg.stop = genes(i).splicegraph{1}(2, idx1__); %%% stop of next exon
                    maxval = inf; 
                    gg = add_reads_from_bam(gg, CFG.bam_fnames, 'exon_track', '', maxval, CFG.read_filter, CFG.var_aware);
                    if gg.strand == '-',
                        gg.tracks = gg.tracks(:, end:-1:1) ;
                    end ;
                    %%% TODO: make the following a configurable
                    if mean(sum(gg.tracks, 1) > 10) < 0.9,
                        idx1__ = [];
                    end;
                end ;
            end ;

			% only take the case closest to an existing splice site
			diff1 = abs(genes(i).splicegraph{1}(1,idx1__) - genes(i).introns{s}(1,j)) ;
			diff2 = abs(genes(i).splicegraph{1}(2,idx1__) - genes(i).introns{s}(1,j)) ;
			diff = min(diff1, diff2) ;
			[tmp, idx1__min] = min(diff) ;
			idx1__ = idx1__(idx1__min) ;
			for idx1_ = idx1__,
				if genes(i).introns{s}(1,j) - 1 - genes(i).splicegraph{1}(1, idx1_) >= CFG.intron_edges.min_exon_len,
					exon_vicinity_cnt1(s) = exon_vicinity_cnt1(s)+1 ;
					genes(i).splicegraph{1}(:,end+1) = genes(i).splicegraph{1}(:, idx1_) ;
					genes(i).splicegraph{1}(2,end) = genes(i).introns{s}(1,j)-1 ; % set exon end to intron start - 1
					genes(i).splicegraph{2}(end+1,end+1) = 0 ;
					genes(i).splicegraph{2}(:,end) = genes(i).splicegraph{2}(:,idx1_) ;
					genes(i).splicegraph{2}(end,:) = genes(i).splicegraph{2}(idx1_,:) ;
					genes(i).splicegraph{3}(:,end+1) = genes(i).splicegraph{3}(:,idx1_) ; % copy from original exon
					genes(i).splicegraph{3}(2,end) = 0 ; % cannot be an end
								
                    assert(all(genes(i).splicegraph{1}(1,:)<=genes(i).splicegraph{1}(2,:)));

					%for jj=1:size(genes(i).splicegraph{1},2)
					%	assert(genes(i).splicegraph{1}(2,jj)-genes(i).splicegraph{1}(1,jj) >= CFG.intron_edges.min_exon_len_remove) ;
					%end 
								
					% check exons whose start coincides with intron end
					genes(i).splicegraph = add_intron(genes(i).splicegraph, size(genes(i).splicegraph{2},1), 0, idx2, 1) ;
								
					inserted.alt_53_prime = inserted.alt_53_prime + 1 ;
                    if CFG.debug,
                        for idx2_ = idx2,
                            fprintf(CFG.fd_log, '%s\talternative_53_prime1\t%c\t%i\t%i\t%i\n', genes(i).chr, genes(i).strand, genes(i).splicegraph{1}(2,idx1_), genes(i).splicegraph{1}(2,end), genes(i).splicegraph{1}(1,idx2_)) ;
                        end; 
                    end;
					intron_used = 1 ;

				end ;
			end ;

			%%% if no proximal exon was found, insert new terminal exon, if wished
			if ~intron_used && CFG.intron_edges.append_new_terminal_exons,
				inserted.new_terminal_exon = inserted.new_terminal_exon + 1 ;

				iregion = [genes(i).introns{s}(1,j) - CFG.intron_edges.append_new_terminal_exons_len; genes(i).introns{s}(1,j) - 1] ; 
				idx_iregion = find(genes(i).introns{s}(2, :) >= iregion(1) & genes(i).introns{s}(2, :) < iregion(2) - 1) ;
				if ~isempty(idx_iregion),
					if ~(length(idx_iregion) == 1), 
						[tmp,idx_iregion_] = max(genes(i).introns{s}(2, idx_iregion)) ;
						idx_iregion = idx_iregion(idx_iregion_) ;
					end ;
					iregion(1) = genes(i).introns{s}(2, idx_iregion) + 1 ;
					assert(iregion(1) < iregion(2)) ;
				end ;

				genes(i).splicegraph{1}(:, end + 1) = iregion ;
				genes(i).splicegraph{2}(end + 1, end + 1) = 0 ;
				genes(i).splicegraph{3}(1, end + 1) = 1 ; % can be a start
				genes(i).splicegraph{3}(2, end) = 0 ; % cannot be an end

                for tmp_idx = idx2,
                    if genes(i).splicegraph{3}(1, tmp_idx) && genes(i).introns{s}(2,j) + 1 <= genes(i).splicegraph{1}(2, tmp_idx),
                        genes(i).splicegraph{1}(1, tmp_idx) = genes(i).introns{s}(2,j) + 1 ;
                        %genes(i).splicegraph{3}(1, tmp_idx) = 0;
                    end ;
				end ;
                assert(all(genes(i).splicegraph{1}(2, :) >= genes(i).splicegraph{1}(1, :)));

				genes(i).splicegraph = add_intron(genes(i).splicegraph, idx2, 1, size(genes(i).splicegraph{2},1), 0) ;
				intron_used = 1 ;
			end ;

			if intron_used,
				continue; 
			end ;
		end ;

		% did not find exons in same gene sharing boundaries with intron end
		% check whether the intron ends in the vicinity of an exon
		if isempty(idx2) 
			idx2__ = find(genes(i).splicegraph{1}(1,:) - CFG.intron_edges.vicinity_region <= genes(i).introns{s}(2,j) & ...
						genes(i).splicegraph{1}(2,:) + CFG.intron_edges.vicinity_region >= genes(i).introns{s}(2,j)) ;

            %%% check, if we can find an exon after the current intron and there is continuous coverage between intron end and exon
            if isempty(idx2__),
                idx2__ = find(genes(i).splicegraph{1}(1, :) > genes(i).introns{s}(2, j), 1, 'first');
                if ~isempty(idx2__),
                    gg = genes(i);
                    gg.strand = strands(s);
                    rm_strands = ~isfield(gg, 'strands');
                    gg.strands = strands(s);
                    gg.start = genes(i).introns{s}(2, j) + 1; %%% start of presumable exon
                    gg.stop = genes(i).splicegraph{1}(2, idx2__); %%% stop of next exon
                    maxval = inf; 
                    gg = add_reads_from_bam(gg, CFG.bam_fnames, 'exon_track', '', maxval, CFG.read_filter, CFG.var_aware);
                    if gg.strand == '-',
                        gg.tracks = gg.tracks(:, end:-1:1) ;
                    end ;
                    %%% TODO: make configurable
                    if mean(sum(gg.tracks, 1) > 10) < 0.9,
                        idx2__ = [];
                    end;
                end ;
            end ;

			% only take the case closest to an existing splice site
			diff1 = abs(genes(i).splicegraph{1}(1,idx2__) - genes(i).introns{s}(2,j)) ;
			diff2 = abs(genes(i).splicegraph{1}(2,idx2__) - genes(i).introns{s}(2,j)) ;
			diff = min(diff1, diff2) ;
			[tmp, idx2__min] = min(diff) ;
			idx2__ = idx2__(idx2__min) ;
			for idx2_ = idx2__,
				if genes(i).splicegraph{1}(2, idx2_) - genes(i).introns{s}(2,j) >= CFG.intron_edges.min_exon_len,
					exon_vicinity_cnt2(s) = exon_vicinity_cnt2(s) + 1 ;
					genes(i).splicegraph{1}(:, end + 1) = genes(i).splicegraph{1}(:, idx2_) ;
					genes(i).splicegraph{1}(1, end) = genes(i).introns{s}(2, j) + 1 ;
					genes(i).splicegraph{2}(end + 1, end + 1) = 0 ;
					genes(i).splicegraph{2}(:, end) = genes(i).splicegraph{2}(:, idx2_) ;
					genes(i).splicegraph{2}(end, :) = genes(i).splicegraph{2}(idx2_, :) ;
					genes(i).splicegraph{3}(:, end + 1) = genes(i).splicegraph{3}(:, idx2_) ; % copy from original exon
					genes(i).splicegraph{3}(1, end) = 0 ; % cannot be a start
						
                    assert(all(genes(i).splicegraph{1}(1,:)<=genes(i).splicegraph{1}(2,:)));

					%for jj=1:size(genes(i).splicegraph{1},2)
					%	assert(genes(i).splicegraph{1}(2,jj)-genes(i).splicegraph{1}(1,jj) >= CFG.intron_edges.min_exon_len_remove) ;
					%end 
					
					genes(i).splicegraph = add_intron(genes(i).splicegraph, idx1, 1, size(genes(i).splicegraph{2},1), 0) ;
					
					inserted.alt_53_prime = inserted.alt_53_prime + 1 ;
                    if CFG.debug,
                        for idx1_ = idx1,
                            fprintf(CFG.fd_log, '%s\talternative_53_prime2\t%c\t%i\t%i\t%i\n', genes(i).chr, genes(i).strand, genes(i).splicegraph{1}(2,idx1_), genes(i).splicegraph{1}(1,end), genes(i).splicegraph{1}(1,idx2_)) ;
                        end ;
                    end;
					intron_used = 1 ;
							
				end ;
			end ;

			%%% if no proximal exon was found, insert new terminal exon, if wished
			if ~intron_used && CFG.intron_edges.append_new_terminal_exons,

                %%% define range of new exon
				iregion = [genes(i).introns{s}(2,j)+1; genes(i).introns{s}(2,j) + CFG.intron_edges.append_new_terminal_exons_len] ;
                %%% find introns starting within new exon
				idx_iregion = find(genes(i).introns{s}(1, :) > iregion(1) + 1 & genes(i).introns{s}(1, :) <= iregion(2)) ;
				if ~isempty(idx_iregion),
					if ~(length(idx_iregion) == 1), 
						[tmp, idx_iregion_] = min(genes(i).introns{s}(1, idx_iregion)) ;
						idx_iregion = idx_iregion(idx_iregion_) ;
					end ;
                    %%% let new exon end at position before next intron starts
					iregion(2) = genes(i).introns{s}(1, idx_iregion) - 1 ; 
					assert(iregion(1) < iregion(2)) ;
				end ;
				inserted.new_terminal_exon = inserted.new_terminal_exon + 1 ;
				genes(i).splicegraph{1}(:, end + 1) = iregion ;
				genes(i).splicegraph{2}(end + 1, end + 1) = 0 ;
				genes(i).splicegraph{3}(1, end + 1) = 0 ; % cannot be a start
				genes(i).splicegraph{3}(2, end) = 1 ; % can be an end

                %%% adapt terminal exon ends if the new intron starts within them
                for tmp_idx = idx1,
                    if genes(i).splicegraph{3}(2, tmp_idx) && genes(i).introns{s}(1,j) - 1 >= genes(i).splicegraph{1}(1, tmp_idx),
                        genes(i).splicegraph{1}(2, tmp_idx) = genes(i).introns{s}(1,j) - 1 ;
                    end;
				end ;
                assert(all(genes(i).splicegraph{1}(2, :) >= genes(i).splicegraph{1}(1, :)));
						
				genes(i).splicegraph = add_intron(genes(i).splicegraph, idx1, 1, size(genes(i).splicegraph{2}, 1), 0) ;
				intron_used = 1 ;
			end ;

			if intron_used,
				continue, 
			end ;
		end ;
		
		if isempty(idx1) || isempty(idx2),
			one_missing(s) = one_missing(s) + 1 ;
			if ~intron_used, unused_introns(end+1) = j ; end ;
			continue ;
		end ;

        %%% TODO: hard coded limit
		if length(idx1) > 20 || length(idx2) > 20,
			multi = multi + 1 ;
			if ~intron_used, 
                unused_introns(end + 1) = j ; 
            end ;
			continue ;
		end ;
		
		%%% both idx1 and idx2 are not empty and are both shorter than 4
		%%% insert exon skips
		for idx1_ = idx1,
			for idx2_ = idx2,
				if genes(i).splicegraph{2}(idx1_, idx2_) == 0,
					inserted.exon_skip = inserted.exon_skip + 1;
								
					%adj_mat = triu(genes(i).splicegraph{2}) ;
					%id1 = find(adj_mat(idx1_,:)) ;
					%if length(id1)==1 && adj_mat(id1, idx2_),
					%		fprintf(CFG.fd_log, '%s\texon_skip\t%c\t%i\t%i\t%i\t%i\t%i\t%i\n', genes(i).chr, genes(i).strand, genes(i).splicegraph{1}(1,idx1_), genes(i).splicegraph{1}(2,idx1_), genes(i).splicegraph{1}(1,id1),	...
					%						genes(i).splicegraph{1}(2,id1), genes(i).splicegraph{1}(1,idx2_), genes(i).splicegraph{1}(2,idx2_)) ;
					%end ;
				end ;
			end 
		end ;

		genes(i).splicegraph = add_intron(genes(i).splicegraph, idx1, 1, idx2, 1) ;
		used_intron=1 ;

		%for i1=idx1,
		%	for i2=idx2,
		%		P(end+1,:)=[i1,i2] ;
		%		genes(i).splicegraph{2}(i1,i2)=1 ;
		%		genes(i).splicegraph{2}(i2,i1)=1 ;
		%	end ;
		%end ;
	end ;

	idx_unused = find(genes(i).introns{s}(2, unused_introns)>=genes(i).start & genes(i).introns{s}(1,unused_introns)<=genes(i).stop) ;
	unused_introns = unused_introns(idx_unused) ;
	if ~isempty(unused_introns),
        % TODO: check this maybe
	    %i
		%unused_introns
		%keyboard
	end ;
	num_unused_introns(i) = num_unused_introns(i)+length(unused_introns) ;
end ;

if print_intermediates,
    one_missing
    multi
    sum(num_unused_introns)
end;

merge_idx = unique(merge_idx, 'rows') ;
rm_map = zeros(1, length(genes)) ;
for i = 1:size(merge_idx,1)
	
	while (rm_map(merge_idx(i,1))==1)
		merge_idx(i,1) = merge_idx(i,1)-1 ;
	end ;

	% merge transcripts
	for j = 1:length(genes(merge_idx(i, 2)).exons)
		genes(merge_idx(i, 1)).transcripts{end + 1} = genes(merge_idx(i, 2)).transcripts{j} ;
		genes(merge_idx(i, 1)).exons{end + 1} = genes(merge_idx(i, 2)).exons{j} ;
	end ;
	% merge intron lists
	for k = 1:length(genes(merge_idx(i, 2)).introns),
		for j = 1:size(genes(merge_idx(i, 2)).introns{k}, 2)
			genes(merge_idx(i, 1)).introns{k}(:, end + 1) = genes(merge_idx(i, 2)).introns{k}(:, j) ;
		end ;
	end ;
	% merge splice graphs
	genes(merge_idx(i, 1)).splicegraph{1} = [genes(merge_idx(i,1)).splicegraph{1} genes(merge_idx(i,2)).splicegraph{1}] ;
	m = size(genes(merge_idx(i, 1)).splicegraph{2}, 1) ;
	n = size(genes(merge_idx(i, 2)).splicegraph{2}, 1) ;
	genes(merge_idx(i, 1)).splicegraph{2}(m + 1 : n + m, m + 1 : n + m) = genes(merge_idx(i, 2)).splicegraph{2} ;
	genes(merge_idx(i, 1)).splicegraph{3} = [genes(merge_idx(i, 1)).splicegraph{3} genes(merge_idx(i, 2)).splicegraph{3}] ;

	% extend start/stop
	genes(merge_idx(i, 1)).start = min(genes(merge_idx(i, 1)).start, genes(merge_idx(i, 2)).start) ;
	genes(merge_idx(i, 1)).stop = max(genes(merge_idx(i, 1)).stop, genes(merge_idx(i, 2)).stop) ;

	%genes(merge_idx(i,1))=build_splice_graph_caller(genes(merge_idx(i,1))) ;
	%genes(merge_idx(i,1))=infer_splice_graph_caller(genes(merge_idx(i,1))) ;

	rm_map(merge_idx(i, 2)) = 1 ;

	inserted.gene_merge = inserted.gene_merge + 1 ;
end ;
genes(rm_map == 1) = [] ;

%size(P,1)/(size(P,1)+both_missing+one_missing)
%both_missing/(size(P,1)+both_missing+one_missing)
%one_missing/(size(P,1)+both_missing+one_missing)

if ( size(genes(i).splicegraph{1},2) > 1) 
    genes(i) = uniquify_splicegraph(genes(i));
end ;

for i = 1:length(genes),
    assert(all(genes(i).splicegraph{1}(1,:) <= genes(i).splicegraph{1}(2, :))) ;
end ;

for ix = 1:length(genes)
	[dummy,exon_order] = sort(genes(ix).splicegraph{1}(1,:),2,'ascend');
	genes(ix).splicegraph{1} = genes(ix).splicegraph{1}(:, exon_order);
	genes(ix).splicegraph{2} = genes(ix).splicegraph{2}(exon_order, exon_order);
	genes(ix).splicegraph{3} = genes(ix).splicegraph{3}(:, exon_order);
end ;

