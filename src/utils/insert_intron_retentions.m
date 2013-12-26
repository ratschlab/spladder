% written by Georg Zeller & Gunnar Raetsch, Mpi Tuebingen, Germany, 2009
function [genes, inserted] = insert_intron_retentions(genes, CFG)
% [genes, inserted] = insert_intron_retentions(genes, CFG)

inserted.intron_retention = 0 ;

%%% form chunks for quick sorting
chunks = [[genes.chr_num]', cast([genes.strand]', 'int32'), [genes.start]', [genes.stop]'];
[chunks, chunk_idx] = sortrows(chunks) ;
assert(issorted(chunks, 'rows'));

strands = '+-';

for i = 1:length(genes),
	genes(i).introns{1} = [] ;
	genes(i).introns{2} = [] ;
end ;

%%% form all possible combinations of contigs and strands --> regions
regions = init_regions(CFG.bam_fnames);
%%% keep only chromosomes found in genes
keepidx = find(ismember([regions.chr_num], unique([genes.chr_num])));
regions = regions(keepidx);

c = 1; cov = 0 ;
%tic
num_introns_added = 0 ;
num_introns = 0 ; 
for j = 1:length(regions)
	chr = regions(j).chr_num;
	s = find(regions(j).strand == strands);
	
	% fill the chunks on the corresponding chromosome
	while c <= size(chunks,1),
		if chunks(c,1) > chr || chunks(c,2)>strands(s),
			break;
		end;
		if ~(chunks(c,1) == chr),
			error('c logic seems wrong\n') ;
		end ;

		if mod(c, 100)==0
			fprintf('\r %i(%i) genes done (found %i new retentions in %i tested introns, %2.1f%%)', c, ...
			size(chunks,1), num_introns_added, num_introns, 100*num_introns_added/num_introns);
		end	;		

		gg = genes(chunk_idx(c));
		gg.strand = strands(s);
		rm_strands = ~isfield(gg, 'strands');
		gg.strands = strands(s);
		maxval = inf; 
		if ~iscell(CFG.bam_fnames)
			gg = add_reads_from_bam(gg, CFG.bam_fnames, 'exon_track', '', maxval, CFG.read_filter);
		else
			for f = 1:length(CFG.bam_fnames),
				gg = add_reads_from_bam(gg, CFG.bam_fnames{f}, 'exon_track', '', maxval, CFG.read_filter);
			end ;
			%%%% sum of mapped_exon_tracks (odd) and spliced exon tracks (even)
			%gg.tracks = [sum(gg.tracks(1:2:end, :), 1); sum(gg.tracks(2:2:end, :), 1)] ;
		end ;
		if gg.strand == '-',
			gg.tracks = gg.tracks(:, end:-1:1) ;
		end ;
		for k = 1:size(gg.splicegraph{1}, 2),
			idx = [gg.splicegraph{1}(1, k):gg.splicegraph{1}(2, k)] - gg.start + 1 ;
			gg.exon_coverage(k,:) = median(sum(gg.tracks(:, idx), 1), 2) ; % median coverage for exon k
			cov = cov + gg.exon_coverage(k, :) ;
		end ;

		%%% check for all vertex-pairs, if respective intron can be retained
		new_retention = zeros(size(gg.splicegraph{2})) ;
		for k = 1:size(gg.splicegraph{1},2) - 1,
			for l = k + 1:size(gg.splicegraph{1}, 2),
				if gg.splicegraph{2}(k,l) > 0,
					num_introns = num_introns + 1 ;
					idx = [gg.splicegraph{1}(2, k) + 1:gg.splicegraph{1}(1, l) - 1] - gg.start + 1 ;
					icov = sum(gg.tracks(:, idx), 1) ; %%% old: intron coverage icov(1,:) --> mapped, icov(2,:) --> spliced; now: sum of both
				%	if median(icov(1, :) + icov(2, :)) > CFG.intron_retention.min_retention_cov && ...
				%		mean(icov(1, :) + icov(2, :) > 0.5 * mean(icov(1, :) + icov(2, :))) > CFG.intron_retention.min_retention_region && ... % fraction of covered positions
				%		max(gg.exon_coverage(k),gg.exon_coverage(l)) / (1e-6 + min(gg.exon_coverage(k), gg.exon_coverage(l))) <= CFG.intron_retention.min_retention_max_exon_fold_diff && ...
				%		mean(icov(1, :) + icov(2, :)) >= CFG.intron_retention.min_retention_rel_cov * (gg.exon_coverage(k) + gg.exon_coverage(l)) / 2 && ...
				%		mean(icov(1, :) + icov(2, :)) <= CFG.intron_retention.max_retention_rel_cov * (gg.exon_coverage(k) + gg.exon_coverage(l)) / 2
					if median(icov) > CFG.intron_retention.min_retention_cov && ...
						mean(icov > (0.5 * mean(icov))) > CFG.intron_retention.min_retention_region && ... % fraction of covered positions
						max(gg.exon_coverage(k),gg.exon_coverage(l)) / (1e-6 + min(gg.exon_coverage(k), gg.exon_coverage(l))) <= CFG.intron_retention.min_retention_max_exon_fold_diff && ...
						mean(icov) >= CFG.intron_retention.min_retention_rel_cov * (gg.exon_coverage(k) + gg.exon_coverage(l)) / 2 && ...
						mean(icov) <= CFG.intron_retention.max_retention_rel_cov * (gg.exon_coverage(k) + gg.exon_coverage(l)) / 2

						new_retention(k,l) = 1 ;
					%	fprintf(log_fd, '%s\tintron_retention\t%c\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%2.1f\n', gg.chr, gg.strand, gg.splicegraph{1}(1,k), gg.splicegraph{1}(2,k), gg.splicegraph{1}(1,l), gg.splicegraph{1}(2,l), ...
					%			floor(median(icov(1,:)+icov(2,:))), floor(gg.exon_coverage(k)), floor(gg.exon_coverage(l)), 100*mean(icov(1,:)+icov(2,:)>0)) ;
						inserted.intron_retention = inserted.intron_retention+1 ;
					end ;
				end ;
			end ;
		end ;
		any_added = 0 ;
		if any(new_retention(:))
			new_retention = expm(new_retention) ;
			while (1),
				any_added=0 ;
				for k = 1:size(new_retention, 2),
					for l = k + 1:size(new_retention, 2),
						if new_retention(k,l) > 0,
							gg.splicegraph = add_intron_retention(gg.splicegraph, k, l) ;
							new_retention(end + 1, end + 1) = 0 ;
							new_retention(k, l) = 0 ;
							any_added = 1 ;
							num_introns_added = num_introns_added + 1 ;
							%fprintf(log_fd, '%s\tintron_retention\t%i\t%i\t%i\t%i\t%i\t%2.1f\n', gg.chr, gg.splicegraph{1}(2,k), gg.splicegraph{1}(1,l), floor(median(icov(1,:)+icov(2,:))), gg.exon_coverage(k), gg.exon_coverage(l), 100*mean(icov(1,:)+icov(2,:)>0)) ;
							break ;
						end ;
					end ;
					if any_added, break ; end ;
				end ;
				[dummy,exon_order] = sort(gg.splicegraph{1}(1, :), 2, 'ascend');
				gg.splicegraph{1} = gg.splicegraph{1}(:, exon_order);
				gg.splicegraph{2} = gg.splicegraph{2}(exon_order, exon_order) ;
				if length(gg.splicegraph) > 2
					gg.splicegraph{3} = gg.splicegraph{3}(:, exon_order) ;
				end
				new_retention = new_retention(exon_order, exon_order) ;
				if ~any_added, break ; end ;
			end ;
		end ;
		if any_added,
			[dummy,exon_order] = sort(gg.splicegraph{1}(1, :), 2, 'ascend');
			gg.splicegraph{1} = gg.splicegraph{1}(:, exon_order);
			gg.splicegraph{2} = gg.splicegraph{2}(exon_order, exon_order) ;
			if length(gg.splicegraph) > 2,
				gg.splicegraph{3} = gg.splicegraph{3}(:, exon_order) ;
			end ;
		end ;
		if rm_strands
			gg = rmfield(gg, 'strands') ;
		end
		gg = rmfield(gg, 'tracks') ;
		gg = rmfield(gg, 'exon_coverage') ;
		genes(chunk_idx(c)) = gg ;
		c = c + 1;
	end;
end;

%[num_introns num_introns_added]
% eof
