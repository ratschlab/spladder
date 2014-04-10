% written by Andre Kahles, Mpi Tuebingen, Germany, 2012

function [genes, inserted] = insert_cassette_exons(genes, CFG)
% [genes, inserted] = insert_cassette_exons(genes, CFG)

inserted.cassette_exon = 0 ;

%%% form chunks for quick sorting
chunks = [[genes.chr_num]', cast([genes.strand]', 'int32'), [genes.start]', [genes.stop]'];
[chunks, chunk_idx] = sortrows(chunks) ;
assert(issorted(chunks, 'rows'));

strands = '+-';

%%% form all possible combinations of contigs and strands --> regions
regions = init_regions(CFG.bam_fnames);
%%% keep only chromosomes found in genes
keepidx = find(ismember([regions.chr_num], unique([genes.chr_num])));
regions = regions(keepidx);

c = 1; cov = 0 ;
num_exons_added = 0 ;
num_exons = 0 ; 

for j = 1:length(regions),
	chr = regions(j).chr_num;
	s = find(regions(j).strand == strands);
	
	% fill the chunks on the corresponding chromosome
	while c <= size(chunks,1),
		if chunks(c,1) > chr || chunks(c,2)>strands(s),
			break;
		end ;
		if ~(chunks(c,1) == chr),
			error('c logic seems wrong\n') ;
		end ;

		if CFG.verbose && mod(c, 100) == 0,
			fprintf('\r %i(%i) genes done (found %i new cassette exons in %i tested intron pairs, %2.1f%%)', c, ...
			size(chunks,1), num_exons_added, num_exons, 100*num_exons_added/num_exons);
		end	;

		gg = genes(chunk_idx(c));
		gg.strand = strands(s);
		rm_strands = ~isfield(gg, 'strands');
		gg.strands = strands(s);
		maxval = inf; 
		if ~iscell(CFG.bam_fnames)
			gg = add_reads_from_bam(gg, CFG.bam_fnames, 'exon_track', '', maxval, CFG.read_filter, CFG.var_aware, CFG.only_primary);
		else
			for f = 1:length(CFG.bam_fnames),
				gg = add_reads_from_bam(gg, CFG.bam_fnames{f}, 'exon_track', '', maxval, CFG.read_filter, CFG.var_aware, CFG.only_primary);
			end ;
			%%%% sum of mapped_exon_tracks (odd) and spliced exon tracks (even)
			%gg.tracks = [sum(gg.tracks(1:2:end,:),1);sum(gg.tracks(2:2:end,:),1)] ;
		end ;
		if gg.strand == '-',
			gg.tracks = gg.tracks(:, end:-1:1) ;
		end ;

        %%% add introns implied by splicegraph to the list
        all_introns = gg.introns{s};
        for k = 1:size(gg.splicegraph{2}, 1),
            for l = k + 1:size(gg.splicegraph{2}, 1),
                if gg.splicegraph{2}(k, l) == 1,
                    all_introns(:, end + 1) = [gg.splicegraph{1}(2, k) + 1; gg.splicegraph{1}(1, l) - 1];
                end;
            end;
        end;
        all_introns = unique(all_introns', 'rows')';
   
        %%% use only relevant introns (inside gene boundaries)
        if ~isempty(all_introns),
            rm_idx = all_introns(2, :) <= gg.start | all_introns(1, :) >= gg.stop;
            all_introns(:, rm_idx) = [];
            clear rm_idx;
        end;

        segment_bounds = sort(unique(all_introns(:)));

		%%% check for all intron-pairs, if exon could exist between them
		new_cassette = zeros(length(all_introns)) ;
        for k = 1:size(all_introns, 2),
            for l = k + 1:size(all_introns, 2),
                if all_introns(2, k) >= all_introns(1, l),
                    continue;
                end;
                %%% only take intron pair, if outer ends are supported by current exons
                if isempty(find(gg.splicegraph{1}(2, :) == all_introns(1, k) - 1)) || isempty(find(gg.splicegraph{1}(1, :) == all_introns(2, l) + 1)),
                    continue;
                end;
                curr_exon = [all_introns(2, k) + 1, all_introns(1, l) - 1];
                %%% do not allow curr_exon to overlap existing exon
                if ~isempty(find(gg.splicegraph{1}(1, :) < curr_exon(2) & gg.splicegraph{1}(2, :) > curr_exon(1)))
                    continue;
                end;

                num_exons = num_exons + 1;

                if ~ismember(curr_exon, gg.splicegraph{1}', 'rows'),
                    idx = [curr_exon(1):curr_exon(2)] - gg.start + 1;
                    exon_cov = sum(gg.tracks(:, idx), 1);
                    pre_segment_end = find(segment_bounds < curr_exon(1) - 1, 1, 'last');
                    if ~isempty(pre_segment_end),
                        pre_segment_cov = sum(gg.tracks(:, [segment_bounds(pre_segment_end) : curr_exon(1) - 1] - gg.start + 1), 1);
                    else
                        pre_segment_cov = sum(gg.tracks(:, 1 : curr_exon(1) - gg.start), 1);
                    end;
                    min_len_pre = min(length(pre_segment_cov), length(exon_cov));

                    aft_segment_start = find(segment_bounds > curr_exon(2) + 1, 1, 'first');
                    if ~isempty(aft_segment_start),
                        aft_segment_cov = sum(gg.tracks(:, [curr_exon(2) + 1 : segment_bounds(aft_segment_start) ] - gg.start + 1), 1);
                    else
                        aft_segment_cov = sum(gg.tracks(:, curr_exon(2) + 1 - gg.start + 1 : end), 1);
                    end;
                    min_len_aft = min(length(aft_segment_cov), length(exon_cov));
                    if mean(exon_cov > 0.2*mean(exon_cov)) > CFG.cassette_exon.min_cassette_region && ...
                       median(exon_cov) > CFG.cassette_exon.min_cassette_cov && ...
                       max(median(exon_cov(end-min_len_aft+1:end)), median(aft_segment_cov(1:min_len_aft))) / min(median(exon_cov(end-min_len_aft+1:end)), median(aft_segment_cov(1:min_len_aft))) - 1 >= CFG.cassette_exon.min_cassette_rel_diff && ...
                       max(median(exon_cov(1:min_len_pre)), median(pre_segment_cov(end-min_len_pre+1:end))) / min(median(exon_cov(1:min_len_pre)), median(pre_segment_cov(end-min_len_pre+1:end))) - 1 >= CFG.cassette_exon.min_cassette_rel_diff,
                        new_cassette(k, l) = 1;
                        inserted.cassette_exon = inserted.cassette_exon + 1;
                    end;
                end;
            end;
        end;
        any_added = 0 ;
		if any(new_cassette(:))
            curr_sg = gg.splicegraph{1};
            for k = 1:size(new_cassette, 2),
                for l = k + 1:size(new_cassette, 2),
                    if new_cassette(k,l) > 0,
                        exons_pre = find(curr_sg(2, :) == (all_introns(1, k) - 1));
                        exons_aft = find(curr_sg(1, :) == (all_introns(2, l) + 1));

                        gg.splicegraph = add_cassette_exon(gg.splicegraph, [all_introns(2, k) + 1, all_introns(1, l) - 1], exons_pre, exons_aft) ;
                        new_cassette(k,l) = 0 ;
                        any_added = 1 ;
                        num_exons_added = num_exons_added + 1 ;
                    end ;
                end ;
            end ;
            [dummy,exon_order] = sort(gg.splicegraph{1}(1,:),2,'ascend');
            gg.splicegraph{1} = gg.splicegraph{1}(:,exon_order);
            gg.splicegraph{2} = gg.splicegraph{2}(exon_order,exon_order) ;
            if length(gg.splicegraph)>2
                gg.splicegraph{3} = gg.splicegraph{3}(:, exon_order) ;
            end;
            if ~any_added, break ; end ;
		end ;
		if any_added,
			[dummy,exon_order] = sort(gg.splicegraph{1}(1,:),2,'ascend');
			gg.splicegraph{1} = gg.splicegraph{1}(:,exon_order);
			gg.splicegraph{2} = gg.splicegraph{2}(exon_order,exon_order) ;
			if length(gg.splicegraph)>2
				gg.splicegraph{3} = gg.splicegraph{3}(:, exon_order) ;
			end;
		end ;
		if rm_strands
			gg = rmfield(gg, 'strands') ;
		end
        %%% clean up gene structure
		gg = rmfield(gg, 'tracks') ;
		genes(chunk_idx(c)) = gg ;
		c = c + 1;
	end;
end;
%[num_exons num_exons_added]
