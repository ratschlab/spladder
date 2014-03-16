function blocks = add_reads_from_bam(blocks, base_dir, which_data, tmp, maxval, filter, var_aware)
% blocks = add_reads_from_bam(blocks, base_dir, which_data, tmp, maxval, filter, var_aware)

% blocks coordinates are assumed to be in closed intervals

if nargin<5
	maxval = 50;
end
if nargin<6
	filter.intron = 20000;
	filter.exon_len = 8;
	filter.mismatch = 1;
end
if nargin < 7,
    var_aware = 0;
end;

if isstruct(which_data)
	%% assume this is the genome info struct
	which_data = tmp;
end

types = strsplit(which_data, ',');

if isempty(types{1})
	fprintf('add_reads_from_bam: nothing to do\n');
	return
end

verbose = 0;
pair = 0;
if any(strcmp(types, 'pair_coverage'))
    pair = 1;
end

clipped = 0;

%reads from both strands are in one sam file
if ~iscell(base_dir),
    filenames = strsplit(base_dir, ',');
else
    filenames = base_dir;
end;

for b = 1:length(blocks),
    clear introns;

	if verbose == 1 && mod(b, 10) == 0,
		fprintf('\radd_exon_track_from_bam: %i(%i)', b, length(blocks));
	end;
	block_len = blocks(b).stop - blocks(b).start + 1;

	%% get data from bam
	if any(strcmp(types, 'exon_track'))
		mapped = 1;
		spliced = 1;
		[introns, coverage, pair_cov] = get_all_data(blocks(b), mapped, spliced, filenames, filter, pair, clipped, var_aware);
	end
	if any(strcmp(types, 'mapped_exon_track'))
		mapped = 1;
		spliced = 0;
		[introns, mapped_coverage, pair_cov] = get_all_data(blocks(b), mapped, spliced, filenames, filter, pair, clipped, var_aware);
	end
	if any(strcmp(types, 'spliced_exon_track'))
		mapped = 0;
		spliced = 1;
		[introns, spliced_coverage, pair_cov] = get_all_data(blocks(b), mapped, spliced, filenames, filter, pair, clipped, var_aware);
	end
	if any(strcmp(types, 'polya_signal_track'))
		mapped = 1;
		spliced = 1;
        clipped = 1;
		[introns, polya_signals, pair_cov] = get_all_data_uncollapsed(blocks(b), mapped, spliced, filenames, filter, clipped, var_aware);
	end
	if any(strcmp(types, 'end_signal_track'))
		mapped = 1;
		spliced = 1;
        clipped = 0;
		[introns, read_end_signals, pair_cov] = get_all_data_uncollapsed(blocks(b), mapped, spliced, filenames, filter, clipped, var_aware);
	end
	if ~exist('introns', 'var')
		% no exon coverage needed at all
		mapped = 0;
		spliced = 1;
        clipped = 0;
		[introns, spliced_coverage, pair_cov] = get_all_data(blocks(b), mapped, spliced, filenames, filter, pair, clipped, var_aware);
	end

    introns = sortrows(introns);
	%% process introns
	if blocks(b).strand == '+'
		introns = introns - blocks(b).start + 1;  %%% ONE based, introns are in closed intervals!
		introns(:, 2) = introns(:, 2) + 1;      %%% augment introns to half open
	elseif ~isempty(introns)
        introns = blocks(b).stop - introns(:, 2:-1:1) + 1;
        introns(:, 2) = introns(:, 2) + 1;
	end;
 
	% add requested data to block
	for j = 1:length(types),
		if strcmp(types{j}, 'pair_coverage'),
			if ~isfield(blocks, 'tracks'),
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end+1,:) = zeros(1, block_len);
			end;
			pair_cov(pair_cov > maxval) = maxval;
			if blocks(b).strand == '+' 
				blocks(b).tracks(end, :) = pair_cov;
			else
				blocks(b).tracks(end, :) = pair_cov(end:-1:1);
			end;	
		elseif strcmp(types{j}, 'exon_track')
			%% add exon track to block
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if ~isfield(blocks, 'tracks')
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end + 1,:) = zeros(1, block_len);
			end;
			coverage(coverage>maxval) = maxval;
			if blocks(b).strand == '+' 
				blocks(b).tracks(end, :) = coverage;
			else
				blocks(b).tracks(end, :) = coverage(end:-1:1);
			end;
		elseif strcmp(types{j}, 'mapped_exon_track')
			%% add mapped exon track to block
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if ~isfield(blocks, 'tracks'),
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end + 1, :) = zeros(1, block_len);
			end;
			mapped_coverage(mapped_coverage>maxval) = maxval;
			if blocks(b).strand == '+' 
				blocks(b).tracks(end, :) = mapped_coverage;
			else
				blocks(b).tracks(end, :) = mapped_coverage(end:-1:1);
			end;
		elseif strcmp(types{j}, 'spliced_exon_track')
			%% add spliced exon track to block
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if ~isfield(blocks, 'tracks'),
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end + 1, :) = zeros(1, block_len);
			end;
			spliced_coverage(spliced_coverage>maxval) = maxval;
			if blocks(b).strand == '+' 
				blocks(b).tracks(end, :) = spliced_coverage;
			else
				blocks(b).tracks(end, :) = spliced_coverage(end:-1:1);
			end;
		elseif strcmp(types{j}, 'intron_track')
			%% add intron track to file
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			intron_list = introns;
			if ~isfield(blocks, 'tracks'),
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end + 1, :) = zeros(1, block_len);
			end;
		
			intron_coverage = zeros(1, block_len);
			if ~isempty(intron_list)
				for k = 1:size(intron_list, 1)
					from_pos = max(1, intron_list(k, 1));
					to_pos = min(block_len, intron_list(k, 2));
					intron_coverage(from_pos:to_pos) = intron_coverage(from_pos:to_pos) + 1;
				end;
			end;
			intron_coverage(intron_coverage>maxval) = maxval;
			blocks(b).tracks(end, :) = intron_coverage;
		elseif strcmp(types{j}, 'intron_list')
			%% add intron list to file
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			intron_list = introns;
		
			if ~isempty(intron_list)
				rm_idx = find(intron_list(:, 1) < 1);
				rm_idx = [rm_idx; find(intron_list(:, 2) > block_len)];
				intron_list(rm_idx, :) = [];
			end;

			% calc number of occurences
			[tmp fidx] = unique(intron_list, 'rows', 'first');
			[intron_list lidx] = unique(intron_list, 'rows', 'last');
			intron_quality = lidx - fidx + 1;
		
			if ~isempty(intron_list),
				[tmp sortidx] = sort(intron_list(:, 2));
				intron_list = intron_list(sortidx, :);
				intron_quality = double(intron_quality(sortidx));
			end;
            
            if isfield(filter, 'mincount'),
                take_map = zeros(1,length(intron_list)) ;
                for kk = 1:length(intron_quality),
                    take_map(kk) = intron_quality(kk) >= filter.mincount ;
                end ;
                %fprintf('dropped %i introns, keeping %i introns\n', sum(take_map==0), sum(take_map~=0)) ;
                intron_list = intron_list(find(take_map ~= 0 ), :) ;
                intron_quality = intron_quality(find(take_map ~= 0 ), :) ;
            end ;
		
			if ~isfield(blocks, 'segment_lists'),
				blocks(b).segment_lists = {intron_list};
				blocks(b).segment_scores = {intron_quality};
			else
				blocks(b).segment_lists{end+1} = intron_list;
				blocks(b).segment_scores{end+1} = intron_quality;
			end;
		elseif strcmp(types{j}, 'polya_signal_track')
			%% add polyA signal track to block
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if ~isfield(blocks, 'tracks'),
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end + 1, :) = zeros(1, block_len);
			end;

            tmp_signals = sparse(zeros(size(polya_signals)));
            %%% get only end positions of reads
            for p_idx = 1:size(polya_signals, 1),
                polya_signals(p_idx, find(polya_signals(p_idx, :), 1, 'last')) = 1;
            end;
            %%% collapse reads
            polya_signals = tmp_signals;
            clear tmp_signals
            polya_signals = sum(polya_signals); 

			if blocks(b).strand == '+',
				blocks(b).tracks(end, :) = polya_signals;
			else
				blocks(b).tracks(end, :) = polya_signals(end:-1:1);
			end;
		elseif strcmp(types{j}, 'end_signal_track')
			%% add read end signal track to block
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if ~isfield(blocks, 'tracks'),
				blocks(b).tracks = zeros(1, block_len);
			else
				blocks(b).tracks(end+1,:) = zeros(1, block_len);
			end;

            tmp_signals = sparse(zeros(size(read_end_signals)));
            %%% get only end positions of reads
            for p_idx = 1:size(read_end_signals, 1),
                tmp_signals(p_idx, find(read_end_signals(p_idx, :), 1, 'last')) = 1;
            end;
            %%% collapse reads
            read_end_signals = tmp_signals;
            clear tmp_signals
            read_end_signals = sum(read_end_signals); 

			if blocks(b).strand=='+',
				blocks(b).tracks(end, :) = read_end_signals;
			else
				blocks(b).tracks(end, :) = read_end_signals(end:-1:1);
			end;
		else 
			fprintf('unknown type of data: %s', types{j})
		end;
	end;
end;

return

function [introns, coverage, pair_cov] = get_all_data(block, mapped, spliced, filenames, filter, pair, clipped, var_aware) 

	block_len = block.stop - block.start + 1;
	% get all data from bam file
	coverage = zeros(1, block_len);
	pair_cov = zeros(1, block_len);
	introns = [];
	for j = 1:length(filenames)
		fname = filenames{j};
		if ~(exist(fname, 'file') == 2)
			fprintf('add_reads_from_bam: did not find file %s\n', fname);
			continue;
		end

		subsample = 1000;%% no subsampling
		if isfield(block, 'subsample_reads'),
			subsample = floor(block.subsample_reads * 1000);
		end;

		contig_name = block.chr;
		strand = block.strand;
		maxminlen = 0;
        %%% check for filter maps -> requires uncollapsed reads
        if isfield(filter, 'masks')
            collapse = 0;
        else
            collapse = 1;
        end;
        %%% get reads from bam file
		if pair
			[coverage_tmp, introns_cell, pair_cov_tmp] = get_reads(fname, contig_name, block.start, block.stop, strand, collapse, subsample, filter.intron, filter.exon_len, filter.mismatch, mapped, spliced, maxminlen, pair, clipped, 0, -1, var_aware);
			pair_cov = pair_cov+double(pair_cov_tmp);
		else
			[coverage_tmp, introns_cell] = get_reads(fname, contig_name, block.start, block.stop, strand, collapse, subsample, filter.intron, filter.exon_len, filter.mismatch, mapped, spliced, maxminlen, pair, clipped, 0, -1, var_aware);
		end;

        %%% compute total coverages
        if isfield(filter, 'maps'),
            %%% reconstruct sparse matrix from uncollapsed read map
            if isempty(coverage_tmp)
              coverage_tmp = zeros(0, block.stop - block.start+ 1 );
            else
              coverage_tmp = sparse(coverage_tmp(1,:)', coverage_tmp(2,:)', ones(size(coverage_tmp,2),1), max(coverage_tmp(1,:)), block.stop - block.start + 1);
            end;

            %%% apply filters
            if isfield(filter.maps, 'repeat_map')
                curr_idx = filter.maps.repeat_map{block.chr_num}(block.start:block.stop);
                rm_idx = sum(coverage_tmp(:, ~curr_idx)) == 0;
                coverage_tmp(rm_idx) = [];
            end;
            if isfield(filter.maps, 'indel_map')
                curr_idx = filter.maps.indel_map{block.chr_num}(block.start:block.stop);
                rm_idx = curr_idx(max(coverage_tmp, [], 2)' == 0)' == 1;
                coverage_tmp(rm_idx) = [];
            end;
            if isfield(filter.maps, 'gene_overlap_map')
                curr_idx = filter.maps.gene_overlap_map{block.chr_num}(block.start:block.stop);
                rm_idx = sum(coverage_tmp(:, curr_idx)) > 0;
                coverage_tmp(rm_idx) = [];
            end;
            coverage = coverage + sum(double(coverage_tmp), 1);
        else
            coverage = coverage + double(coverage_tmp);
        end;

		if iscell(introns_cell),
			introns = [introns; double([introns_cell{:}]')];
		else
			introns = [introns; double(introns_cell)'];
		end;
	end
return

function [introns, coverage, pair_cov] = get_all_data_uncollapsed(block, mapped, spliced, filenames, filter, var_aware) 

	block_len = block.stop - block.start + 1;
	% get all data from bam file
	coverage = [];
	pair_cov = zeros(1, block_len);
	introns = [];
	for j = 1:length(filenames)
		fname = filenames{j};
		if ~(exist(fname, 'file')==2)
			fprintf('add_reads_from_bam: did not find file %s\n', fname);
			continue;
		end;

		subsample = 1000;%% no subsampling
		if isfield(block, 'subsample_reads'),
			subsample = floor(block.subsample_reads*1000);
		end

		contig_name = block.chr;
		strand = block.strand;
		collapse = 0;
		maxminlen = 0;
        pair = 0;
        coverage_tmp = get_reads(fname, contig_name, block.start, block.stop, strand, collapse, ...
                                                 subsample, filter.intron, filter.exon_len, filter.mismatch, ...
                                                 mapped, spliced, maxminlen, pair, clipped, 0, -1, var_aware);
		coverage = [coverage; sparse(coverage_tmp(1,:), coverage_tmp(2, :), 1)];
	end
return

