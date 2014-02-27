function path=cwd ;
% path=cwd ;
% 
% removes the automounter prefixes

path = pwd ;

dirs{1,1} = '/.amd_mnt/huangho/export/kwaid0' ;
dirs{1,2} = '/fml/ag-raetsch' ;

dirs{2,1} = '/.amd_mnt/yangtse/export/altai1' ;
dirs{2,2} = '/agbs/cluster' ;

dirs{3,1} = '/.amd_mnt/steinlach/export/home/kyb/agbs' ;
dirs{3,2} = '/kyb/agbs' ;

for i=1:size(dirs,1)
  if length(path)>=length(dirs{i,1})
    if isequal(dirs{i,1},path(1:length(dirs{i,1}))),
      path(1:length(dirs{i,1}))= [] ;
      path = [dirs{i,2} path] ;
    end ;
  end ;
end ;
