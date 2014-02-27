function [converted_intervals, rm_idx] = convert_strain_pos_intervals(chr_num, intervals, source_strain, target_strain) ;
% converted_intervals = convert_strain_pos_intervals(chr_num, intervals, source_strain, target_strain) ;

converted_intervals = convert_strain_pos(chr_num, intervals, source_strain, target_strain) ;

rm_idx=[] ;

if all(converted_intervals(:)>0),
    return ;
end ;

for i=1:size(converted_intervals,1),
    if any(converted_intervals(i,:)<=0),
        map = convert_strain_pos(chr_num, intervals(i,1):intervals(i,2), source_strain, target_strain) ; 
        if all(map<=0),
            rm_idx(end+1) = i ; 
        else
            converted_intervals(i,1) = min(map(map>0)) ;
            converted_intervals(i,2) = max(map(map>0)) ;
        end ;
    end ;
end ;
converted_intervals(rm_idx,:) = [] ;
