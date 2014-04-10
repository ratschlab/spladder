function detect_truncations(CFG)
% detect_truncations(CFG)

if ~isfield(CFG, 'genes'),
    load(CFG.anno_fname, 'genes');
else
    genes = CFG.genes ;
end;

inserted.truncation = 0 ;


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

c = 1;
num_truncations_added = 0 ;
num_introns = 0 ; 
truncations_left = {};
truncations_right = {};

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

		if mod(c, 50)==0
			fprintf('\r %i(%i) genes done (found %i truncations in %i tested introns, %2.1f%%)', c, ...
			size(chunks,1), num_truncations_added, num_introns, 100*num_truncations_added/num_introns);
		end	;		

		gg = genes(chunk_idx(c));
        %fprintf(1, '%i\r', chunk_idx(c));
		gg.strand = strands(s);
		gg.strands = strands(s);
		maxval = inf; 
		if ~iscell(CFG.bam_fnames)
			gg = add_reads_from_bam(gg, CFG.bam_fnames, 'exon_track', '', maxval, CFG.read_filter, CFG.var_aware, CFG.only_primary);
		else
			for f = 1:length(CFG.bam_fnames),
				gg = add_reads_from_bam(gg, CFG.bam_fnames{f}, 'exon_track', '', maxval, CFG.read_filter, CFG.var_aware, CFG.only_primary);
			end ;
		end ;
		if gg.strand == '-',
			gg.tracks = gg.tracks(:, end:-1:1) ;
		end ;
		for k = 1:size(gg.splicegraph{1}, 2),
			idx = [gg.splicegraph{1}(1, k):gg.splicegraph{1}(2, k)] - gg.start + 1 ;
			gg.exon_coverage(k,:) = mean(sum(gg.tracks(:, idx), 1), 2) ; % mean coverage for exon k (use median for retention ...)
		end ;

		%%% check for all vertex-pairs, if respective intron can be retained
		new_truncation_left = zeros(size(gg.splicegraph{2})) ;
		new_truncation_right = zeros(size(gg.splicegraph{2})) ;
		for k = 1:size(gg.splicegraph{1},2) - 1,
			for l = k + 1:size(gg.splicegraph{1}, 2),
				if gg.splicegraph{2}(k,l) > 0,
					num_introns = num_introns + 1 ;
					idx = [gg.splicegraph{1}(2, k) + 1:gg.splicegraph{1}(1, l) - 1] - gg.start + 1 ;
					icov = sum(gg.tracks(:, idx), 1) ;
                    %%% check if left exon can be extended into the intron 
                    if size(idx, 2) > 10 && mean(gg.exon_coverage(k)) > 1 && mean(icov(1:10)) > 0.75 * mean(gg.exon_coverage(k)),
                        %%% find end point
                        tmp = conv(icov, ones(1, 10)/10);
                        idx2 = find(tmp(5:end-5) < 0.1 * gg.exon_coverage(k), 1, 'first');
                        lens = gg.splicegraph{1}(2, :) - gg.splicegraph{1}(1, :) + 1;
                        if ~isempty(idx2) && mean(icov(idx2:end)) < 0.1 * gg.exon_coverage(k) && mean(icov(1:idx2 - 1)) > 0.5 * gg.exon_coverage(k) && sum(gg.exon_coverage(l:end)' .* lens(l:end)) / sum(lens(l:end)) < 0.02 * gg.exon_coverage(k),
                            val = idx(1) + idx2 + gg.start - 2;
                            %%% check if there exists an end-terminal exon that includes the newly found truncation
                            idx3 = find(gg.splicegraph{3}(2, :) & abs(gg.splicegraph{1}(1, :) - gg.splicegraph{1}(1, k)) <= 5);
                            idx3 = idx3(idx3 ~= k);
                            if isempty(idx3) || val > max(gg.splicegraph{1}(2, idx3)) + 10, 
                                new_truncation_left(k, l) = idx(1) + idx2 + gg.start - 2;
                                inserted.truncation = inserted.truncation + 1 ;
                                fprintf(1, 'adding truncation-left in gene %s at position %i\n', gg.name, new_truncation_left(k, l));
                                truncations_left{end + 1, 1} = gg.name;
                                truncations_left{end, 2} = new_truncation_left(k, l);
                                %keyboard;
                            end;
                        end;
                    end;
                    %%% check if right exon can be extended into the intron 
                    if size(idx, 2) > 10 && mean(gg.exon_coverage(l)) > 1 && mean(icov(end-10:end)) > 0.75 * mean(gg.exon_coverage(l)),
                        %%% find end point
                        tmp = conv(icov, ones(1, 10)/10);
                        idx2 = find(tmp(5:end-5) < 0.1 * mean(gg.exon_coverage(l)), 1, 'last');
                        lens = gg.splicegraph{1}(2, :) - gg.splicegraph{1}(1, :) + 1;
                        if ~isempty(idx2) && mean(icov(1:idx2)) < 0.1 * gg.exon_coverage(l) && mean(icov(idx2:end)) > 0.5 * gg.exon_coverage(l) && sum(gg.exon_coverage(1:k)' .* lens(1:k)) / sum(lens(1:k)) < 0.02 * gg.exon_coverage(l),
                            val = idx(1) + idx2 + gg.start - 2;
                            %%% check if there exists a start-terminal exon that includes the newly found truncation
                            idx3 = find(gg.splicegraph{3}(1, :) & abs(gg.splicegraph{1}(2, :) - gg.splicegraph{1}(2, l)) <= 5);
                            idx3 = idx3(idx3 ~= l);
                            if isempty(idx3) || val < min(gg.splicegraph{1}(1, idx3)) - 10, 
                                new_truncation_right(k, l) = idx(1) + idx2 + gg.start - 2;
                                inserted.truncation = inserted.truncation + 1 ;
                                fprintf(1, 'adding truncation-right in gene %s at position %i\n', gg.name, new_truncation_right(k, l));
                                truncations_right{end + 1, 1} = gg.name;
                                truncations_right{end, 2} = new_truncation_right(k, l);
                                %keyboard;
                            end;
                        end;
                    end;
				end ;
			end ;
		end ;
		%if any(new_truncation_left(:))
        %    t = find(new_truncation_left);
        %    for tt = 1:length(t),
        %        [k, l] = ind2sub(size(new_truncation_left), t(tt));
        %        gg.splicegraph{1}(:, end + 1) = [gg.splicegraph{1}(1, k); new_truncation_left(k,l) - 1];
        %        gg.splicegraph{2}(end + 1, end + 1) = 0;
        %        gg.splicegraph{2}(end + 1, :) = gg.splicegraph{2}(k, :);
        %        gg.splicegraph{2}(:, end + 1) = gg.splicegraph{2}(:, k);
        %        num_truncations_added = num_truncations_added + 1 ;
        %        gg.splicegraph{3}(:, end + 1) = [0; 1];
        %    end;
        %    [dummy,exon_order] = sort(gg.splicegraph{1}(1, :), 2, 'ascend');
        %    gg.splicegraph{1} = gg.splicegraph{1}(:, exon_order);
        %    gg.splicegraph{2} = gg.splicegraph{2}(exon_order, exon_order) ;
        %    if length(gg.splicegraph) > 2
        %        gg.splicegraph{3} = gg.splicegraph{3}(:, exon_order) ;
        %    end
		%end ;
		%if any(new_truncation_right(:))
        %    t = find(new_truncation_right);
        %    for tt = 1:length(t),
        %        [k, l] = ind2sub(size(new_truncation_right), t(tt));
        %        gg.splicegraph{1}(:, end + 1) = [new_truncation_right(k,l) + 1; gg.splicegraph{1}(2, l)];
        %        gg.splicegraph{2}(end + 1, end + 1) = 0;
        %        gg.splicegraph{2}(end + 1, :) = gg.splicegraph{2}(k, :);
        %        gg.splicegraph{2}(:, end + 1) = gg.splicegraph{2}(:, k);
        %        num_truncations_added = num_truncations_added + 1 ;
        %        gg.splicegraph{3}(:, end + 1) = [0; 1];
        %    end;
        %    [dummy,exon_order] = sort(gg.splicegraph{1}(1, :), 2, 'ascend');
        %    gg.splicegraph{1} = gg.splicegraph{1}(:, exon_order);
        %    gg.splicegraph{2} = gg.splicegraph{2}(exon_order, exon_order) ;
        %    if length(gg.splicegraph) > 2
        %        gg.splicegraph{3} = gg.splicegraph{3}(:, exon_order) ;
        %    end
		%end ;
		c = c + 1;
	end;
end;

fd_out = fopen(CFG.out_fname, 'w');
for i = 1:size(truncations_left, 1),
    fprintf(fd_out, '%s\t%i\tL\n', truncations_left{i, :});
end;
for i = 1:size(truncations_right, 1),
    fprintf(fd_out, '%s\t%i\tR\n', truncations_right{i, :});
end;
fclose(fd_out);
