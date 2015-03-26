function genes_cell2struct(anno_fname)
% GENES_CELL2STRUCT   Converts genes stored as a cell to struct.
%
%   genes_cell2struct(anno_fname)
%
%   -- input --
%   anno_fname:   name of file where genes as cell are stored
%
%   -- output --
%   genes as a struct
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2011 Max Planck Society
%

load(anno_fname, 'genes');
if iscell(genes)
  genes_cell = genes;
  clear genes;
  for g = 1:length(genes_cell), 
    gene = genes_cell{g};
    for e = 1:length(gene.exons)
      gene.exons{e} = double(gene.exons{e});
    end    
    gene.exons = reshape(gene.exons, 1, length(gene.exons));
    gene.id = double(gene.id);
    gene.start = double(gene.start);
    gene.stop = double(gene.stop);
    genes(g) = gene;
  end
  save(anno_fname, 'genes');
end
