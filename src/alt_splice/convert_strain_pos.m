function converted_pos = convert_strain_pos(chr_num, pos, source_strain, target_strain)
% converted_pos = convert_strain_pos(chr_num, pos, source_strain, target_strain)
%
% global variables mapping_tables mapping_strains

% handle the case where the strain name also contains the sample name
if any(source_strain==':')
    items=separate(source_strain, ':') ;
    assert(length(items)==2) ;
    source_strain=items{1};
    source_sample=items{2} ;
end ;

if any(target_strain==':')
    items=separate(target_strain, ':') ;
    assert(length(items)==2) ;
    target_strain=items{1};
    target_sample=items{2} ;
end ;

source_strain = lower(source_strain) ;
target_strain = lower(target_strain) ;

if isequal(source_strain, target_strain),
    converted_pos=pos ;
    return ;
end ;

%mapfile = '/agbs/cluster/stegle/data/arabcall/out/hdf5/maps.hdf5';
%mapfile = '/agbs/agkbshare/projects/20strains/arabcall/out/hdf5/maps.hdf5';
%mapfile = '/kyb/agbs/stegle/work/20ArabStrains/genome_coordinates/out/hdf5/maps.gz.hdf5' ;

%mapfile = '/agbs/agkbshare/projects/20strains/arabcall/out/hdf5_v7c/maps.gz.hdf5' ;
%if fexist('/tmp/global2/A_thaliana_magic/hdf5_v7c/maps.gz.hdf5'),
%    mapfile = '/tmp/global2/A_thaliana_magic/hdf5_v7c/maps.gz.hdf5' ;
%end ;
%mapfile = '/agbs/agkbshare/projects/20strains/arabcall/out/hdf5_v7c/maps.gz.hdf5'; 
mapfile= '/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/results/maps/hdf5_v7c/maps.gz.hdf5' ;
%if fexist('/tmp/global2/A_thaliana_magic/hdf5_v7c/maps_all.gz.hdf5'),
%    mapfile = '/tmp/global2/A_thaliana_magic/hdf5_v7c/maps_all.gz.hdf5' ;
%end ;


if isequal(source_strain, target_strain) && isequal(source_strain, 'col_0'),
    % nothing to do, nothing to gain ... 
    % for non-col_0 strains we do the transformation via col_0 anyway to obtain mappable positions
    converted_pos=pos ;
    return ;
end ;
if ~isequal(source_strain, 'col_0') && ~isequal(target_strain, 'col_0'),
    converted_pos1 = convert_strain_pos(chr_num, pos, source_strain, 'col_0') ;
    converted_pos = convert_strain_pos(chr_num, converted_pos1, 'col_0', target_strain) ;
    return ;
end ;
if isequal(source_strain, 'col_0') 
    strain=target_strain ;
else    
    assert(isequal(target_strain, 'col_0')) 
    strain=source_strain ;
end ;

max_chr_num = 5 ;
assert(chr_num<=max_chr_num) ;

global mapping_tables mapping_strains 

idx=strmatch([source_strain ':' target_strain], mapping_strains , 'exact') ;
assert(length(idx)<=1) ;
if isempty(idx),
    mapping_strains{end+1} = [source_strain ':' target_strain];
    mapping_tables{end+1, max_chr_num} = [] ;
    idx = size(mapping_tables,1) ;
    fprintf('created conversion table entry %s\n', mapping_strains{end});
end ;
if isempty(mapping_tables{idx, chr_num})
    fprintf('Filling conversion table entry %s chromosome %i\n', mapping_strains{idx}, chr_num);

    info = hdf5info(mapfile);
	strain_idx = find(strcmp(['/' lower(strain)], lower({info.GroupHierarchy.Groups.Name})));
	map = hdf5read(info.GroupHierarchy.Groups(strain_idx).Datasets(chr_num));
	% map(1, :) contains the chromosome number
	% map(2, :) contains the col_0 positions 
	% map(3, :) contains the strain positions

    if isequal(target_strain, 'col_0'),
        map = map(2:3, :);
        inverse_map = -ones(1, max(map(2,:))); 
        map = map(:, map(1, :)~=-1);
        for j = 1:size(map, 2)
            if map(2, j)>0,
                inverse_map(map(2, j)) = map(1, j);
            end
        end	
        mapping_tables{idx, chr_num} = inverse_map ;
    else
        map = map(2:3, :);
        inverse_map = -ones(1, max(map(1,:))); 
        map = map(:, map(2, :)~=-1);
        for j = 1:size(map, 2)
            if map(1, j)>0,
                inverse_map(map(1, j)) = map(2, j);
            end
        end	
        mapping_tables{idx, chr_num} = inverse_map ;

        assert(isequal(source_strain, 'col_0')) ;
    end ;
end ;

mapping_table = mapping_tables{idx, chr_num} ;
assert(~isempty(mapping_table)) ;

converted_pos = map_col0(pos, mapping_table);

return

function array = map_col0(array, map)
% array = map_col0(array, map)

	for j = 1:size(array, 1)
		for k = 1:size(array, 2)
			if j>0 && k>0&& (array(j,k)>0) 
                if array(j,k)<=length(map),
                    array(j,k) = map(array(j,k));
                else
                    array(j,k) = -1 ;
                end ;
			end
		end
	end

return

