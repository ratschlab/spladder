% (author) Georg Zeller & Gunnar Raetsch, Mpi Tuebingen, Germany, 2009
% (author) Andre Kahles, MSKCC NYC, USA, 2013 

function genes = append_intron_list(genes, fn_bam, CFG)
% genes = append_intron_list(genes, fn_bam, CFG)

if ~isfield(CFG.intron_filter, 'intron'),
    CFG.intron_filetr.intron = 20000;
end 
if ~isfield(CFG.intron_filter, 'exon_len'),
    CFG.intron_filter.exon_len = 20; % relatively specific
end
if ~isfield(CFG.intron_filter, 'mismatch'),
    CFG.intron_filter.mismatch = 2;
end ;
if ~isfield(CFG.intron_filter, 'mincount'),
    CFG.intron_filter.mincount = 5 ;
end ;

%%% form chunks for quick sorting
chunks = [[genes.chr_num]', [genes.strand]', [genes.start]', [genes.stop]'];
[chunks, chunk_idx] = sortrows(chunks) ;
assert(issorted(chunks, 'rows'));

strands = '+-';

for i = 1:length(genes),
	genes(i).introns{1} = [] ;
	genes(i).introns{2} = [] ;
end ;

%%% collect all possible combinations of contigs and strands
regions = init_regions(fn_bam);
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

		if CFG.verbose && mod(c, 50) == 0,
			fprintf('%i(%i) genes done (%i introns taken)\n', c, size(chunks,1), num_introns_filtered);
		end ;			

		gg = genes(chunk_idx(c));
		gg.strand = strands(s);
        gg.start = gg.start - 5000 ;
        gg.stop = gg.stop + 5000 ;

		maxval = inf; 
        if ~iscell(fn_bam),
            gg = add_reads_from_bam(gg, fn_bam, 'intron_list', '', maxval, CFG.intron_filter);
        else
            % merge intron lists of several bam files
            segments = [] ;
            for f = 1:length(fn_bam),
                gg = add_reads_from_bam(gg, fn_bam{f}, 'intron_list', '', maxval, CFG.intron_filter);
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
            segments(rm_idx,: ) = [] ;
            idx = find(segments(:, 3) > CFG.intron_filter.mincount) ;
            gg.segment_lists = {segments(idx, 1:2)} ;
            gg.segment_scores = {segments(idx, 3)} ;
        end ;
        num_introns_filtered = num_introns_filtered + size(gg.segment_lists{1}, 1) ;

        %% gg.segment_lists{1} => intron list, ONE based, half open
		if strands(s)=='+',
			genes(chunk_idx(c)).introns{s} = double([gg.segment_lists{1}(:, 1)'; gg.segment_lists{1}(:, 2)'-1]+gg.start-1) ; % intron list is one based, closed, plus strand relative counting 
		else
			genes(chunk_idx(c)).introns{s} = double(gg.stop-[gg.segment_lists{1}(:, 2)'-1; gg.segment_lists{1}(:, 1)']+1) ; 
		end;
		c = c + 1;
	end;
end;
