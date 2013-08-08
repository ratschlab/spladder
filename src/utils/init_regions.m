function regions=init_regions(fn_bams, organism)
% regions=init_regions(fn_bams, organism)

regions(1).chr = [];
regions(1).chr_num = [];
regions(1).strand = [];
regions(1).start = [];
regions(1).stop = [];
regions(1).offset =[];
regions(1).id = [];
regions(1).organism = [];

for i = 1:length(fn_bams),
    if exist(fn_bams{i}, 'file') == 0,
        continue;
    else
      regions = regions([]) ;
      header_info = get_header(fn_bams{i});
      strands = '+-';
      for c = 1:size(header_info, 1),
        %+strand
        for s = 1:length(strands),
          region = [];
          region.chr = header_info{c, 1};
          region.chr_num = c;
          region.strand = strands(s);
          region.start = 1;
          region.stop = header_info{c, 2}(1) ;
          switch strands(s)
           case '+',
            region.offset      = region.start - 1;
           case '-',
            region.offset      = region.stop + 1;
          end
          region.id = (c-1)*2+s;
          region.organism = [];
          if exist('organism', 'var')==1,
            region.organism = organism;
          end;
          regions(region.id) = region;
        end;
      end;
    end;
    break;
end;
