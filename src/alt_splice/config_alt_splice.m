

%%% add paths to the local workspace
addpath ~/git/tools/rproc.torque
addpath ~/git/tools/utils_matlab
addpath ~/git/tools/ngs/

addpath ~/git/projects/2012/mGene_core/utils/
addpath ~/git/tools/genomes/splicegraphs
addpath ~/git/tools/genomes/splicegraphs/detect_altsplice/
addpath ~/git/software/spladder/src
addpath ~/git/software/spladder/src/utils
addpath ~/git/software/spladder/mex
addpath ~/git/software/spladder/bin

addpath ./writer
addpath ~/git/projects/2013/THCA/alternative_splicing/writer

if ~exist('merge_strategy', 'var'),
    merge_strategy = 'merge_graphs';
end;


CFG.global.base_dir = '/cbio/grlab/projects/TCGA/Fagin';

CFG.samples = {};
CFG.gene_list = {};
CFG.list_config = {};
CFG.list_bam = {};
CFG.strains = {};

%%% define list of genes of interest (can be annotation or any other gene list)
global_gene_list = '/cbio/grlab/projects/TCGA/PanCancer/annotation/gencodeV14.v7.protein_coding.no_cds.mat' ;

is_half_open = 0;

%%% reference genome used for the analysis
reference_strain = 'hg19' ;
reference_gio = '/cbio/grlab/projects/TCGA/genome/hg19.fa.gio/genome.config';

experiments = {'THCA'};

for e_idx = 1:length(experiments),

    %%% define a sample list based on the respective experiment
    %%% experiments can be different conditions, cancer types, ...
    [ret_status, files] = system(sprintf('ls -1 %s/../PanCancer/rna/bam-08-19-2013/%s/*.v4.star.bam', CFG.base_dir, experiments{e_idx}));
    %files = strsplit(files, sprintf('\n'), true); %%% Octave
    files = separate(files, sprintf('\n'), 1);

    %%% set samples / BAM files
    for i = 1:length(files),
        if exist(strrep(files{i}, 'star.bam', 'done'), 'file'),
            samples{end + 1} = [experiments{e_idx} '_' regexprep(regexprep(files{i}, '.*/', ''), '.bam', '')]; 
            list_bam{end + 1} = files{i};
            gene_list{end + 1} = global_gene_list;
            list_config{end + 1} = reference_gio;
            %%% set genome and sample paths
            if strcmp(merge_strategy, 'merge_graphs')
                strains{end + 1} = [reference_strain ':' samples{end}] ; 
            elseif strcmp(merge_strategy, 'merge_bams')
                list_config = reference_gio;
                strains{end + 1} = [reference_strain ':' samples{end}] ; 
            else
                error(sprintf('Unknown merge strategy: %s\n', merge_strategy));
            end;
        else
            fprintf(1, '%s not finished yet\n', files{i});
        end;
    end;
end ;

%%% set result dir
CFG.global.result_dir = sprintf('%s/alternative_splicing/THCA.coding/', CFG.global.base_dir);
unix(sprintf('mkdir -p %s', CFG.global.result_dir));

%confidence_levels = [3 2 1 0] ;
confidence_levels = [3];

transform_coordinates = 0;

use_rproc = 0;
replicate_idxs = 1;
same_genestruct_for_all_samples = 1;
curate_alt_prime_events = 1;

%%%%%%%%% configurations for spladder %%%%%%%%%
if ~isfield(CFG, 'spladder'),
    CFG.spladder = struct();
end;
CFG.spladder.validate_splicegraphs = 1;
CFG.spladder.gen_graph_do_prune = 0;
CFG.spladder.gen_graph_do_gen_isoforms = 0;

