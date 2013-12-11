% (author) Georg Zeller & Gunnar Raetsch, Mpi Tuebingen, Germany, 2009
% (author) Andre Kahles, MSKCC NYC, USA, 2013 

function introns = get_intron_list(genes, CFG)
% introns = get_intron_list(genes, CFG)

%%% form chunks for quick sorting
chunks = [[genes.chr_num]', cast([genes.strand]', 'int32'), [genes.start]', [genes.stop]'];
[chunks, chunk_idx] = sortrows(chunks) ;
assert(issorted(chunks, 'rows'));

strands = '+-';

introns = cell(length(genes), 2);;

%%% collect all possible combinations of contigs and strands
regions = init_regions(CFG.bam_fnames);
keepidx = find(ismember([regions.chr_num], unique([genes.chr_num])));
regions = regions(keepidx);

c = 1;
num_introns_filtered = 0 ;

for j = 1:length(regions)
	chr = regions(j).chr_num;
	s = find(regions(j).strand == strands);
	
	% fill the chunks on the corresponding chromosome
	while c <= size(chunks,1),
		if chunks(c,1) > chr || chunks(c,2)>strands(s),
			break;
		end
		if ~(chunks(c,1) == chr),
			error('c logic seems wrong\n') ;
		end ;

		if CFG.verbose == 1 && mod(c, 100) == 0,
			fprintf('%i (%i) genes done (%i introns taken)\n', c, size(chunks,1), num_introns_filtered);
		end ;			

		gg = genes(chunk_idx(c));
		gg.strand = strands(s);
        gg.start = max(gg.start - 5000, 1) ;
        gg.stop = gg.stop + 5000 ;

		maxval = inf; 
        if ~iscell(CFG.bam_fnames),
            gg = add_reads_from_bam(gg, CFG.bam_fnames, 'intron_list', '', maxval, CFG.read_filter);
        else
            % merge intron lists of several bam files
            segments = [] ;
            for f = 1:length(CFG.bam_fnames),
                gg = add_reads_from_bam(gg, CFG.bam_fnames{f}, 'intron_list', '', maxval, CFG.read_filter);
                segments = [segments; gg.segment_lists{end} gg.segment_scores{end}] ;
            end ;
            segments = sortrows(segments) ;

            rm_idx = [] ;
            for i = 1:size(segments,1) - 1,
                if segments(i,1) == segments(i + 1, 1) && segments(i, 2) == segments(i + 1, 2),
                    rm_idx(end + 1) = i ;
                    segments(i + 1, 3) = segments(i, 3) + segments(i + 1, 3) ;
                end ;
            end ;
            segments(rm_idx, :) = [] ;

            if isempty(segments),
                introns{chunk_idx(c), s} = [];
                c = c + 1;
                continue;
            else
                idx = find(segments(:, 3) > CFG.read_filter.mincount) ;
                gg.segment_lists = {segments(idx, 1:2)} ;
                gg.segment_scores = {segments(idx, 3)} ;
            end;
        end ;
        num_introns_filtered = num_introns_filtered + size(gg.segment_lists{1}, 1) ;

        %% gg.segment_lists{1} => intron list, ONE based, half open
		if strands(s)=='+',
            introns{chunk_idx(c), s} = double([gg.segment_lists{1}(:, 1)'; gg.segment_lists{1}(:, 2)'-1]+gg.start-1) ; % intron list is one based, closed, plus strand relative counting 
		else
			introns{chunk_idx(c), s} = double(gg.stop-[gg.segment_lists{1}(:, 2)'-1; gg.segment_lists{1}(:, 1)']+1) ; 
		end;
		c = c + 1;
	end;
end;
