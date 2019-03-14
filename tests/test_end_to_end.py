import gzip
import h5py
import os
import scipy as sp
import glob

import pytest
import pickle


data_dir = os.path.join(os.path.dirname(__file__), 'data')
sample_dir = os.path.join(os.path.dirname(__file__), 'testcase1')
from spladder import spladder


def _compare_hdf5(expected, actual):
    for k in expected:
        if isinstance(expected[k], h5py._hl.group.Group):
            for l in expected[k]:
                assert sp.all(expected[k][l][:] == actual[k][l][:])
        else:
            if k in ['psi']:
                assert sp.all(expected[k][:].astype('str') == actual[k][:].astype('str'))
            else:
                assert sp.all(expected[k][:] == actual[k][:])

def _compare_ps(expected, actual):

    e_str = []
    for line in expected:
        if line.startswith('%'):
            continue
        e_str.append(line)
    a_str = []
    for line in actual:
        if line.startswith('%'):
            continue
        a_str.append(line)
    assert ''.join(e_str) == ''.join(a_str)

def _assert_files_equal_testing(e, a):

    da = sp.loadtxt(a, dtype='str', delimiter='\t')
    de = sp.loadtxt(e, dtype='str', delimiter='\t')

    ### check header
    assert sp.all(da[0, :] == de[0, :])
    da = da[1:, ]
    de = de[1:, ]

    ### check text cols
    assert sp.all(da[:, [0, 1]] == de[:, [0, 1]])
    da = da[:, 2:]
    de = de[:, 2:]

    ### check p-values (up to certain precision)
    da = sp.around(da.astype('float'), decimals=6)
    de = sp.around(de.astype('float'), decimals=6)
    assert sp.all(da == de)

def _assert_files_equal(expected_path, actual_path):
    def o(f):
        if f.endswith('.gz'):
            return gzip.open(f, 'rb')
        elif f.endswith('.hdf5'):
            return h5py.File(f, 'r')
        elif f.endswith('.ps'):
            return open(f, 'r')
        elif f.endswith('.pickle'):
            return open(f, 'rb')
        else:
            return open(f, 'rb')

    with o(expected_path) as e:
        with o(actual_path) as a:
            if expected_path.endswith('.hdf5'):
                _compare_hdf5(e, a)
            elif expected_path.endswith('.ps'):
                _compare_ps(e, a)
            elif expected_path.endswith('.pickle'):
                ta = pickle.load(a, encoding='latin1')
                te = pickle.load(e, encoding='latin1')
                if len(ta) == 0 or isinstance(ta[0], sp.int64):
                    assert sp.all(ta == te)
                else:
                    assert sp.all([_compare_gene(_[0], _[1]) for _ in zip(ta, te)])
            elif os.path.basename(expected_path).startswith('test_results'):
                _assert_files_equal_testing(e, a)
            else:
                assert e.read() == a.read(), 'actual and expected content differ!\nactual path: %s\nexpected path: %s\n' % (expected_path, actual_path)

def _codeUTF8(s):
    return s.view(sp.chararray).encode('utf-8')

def _compare_gene(a, b):
    if sp.issubdtype(a.strain.dtype, sp.str_):
        _astrain = _codeUTF8(a.strain)
    else:
        _astrain = a.strain
    if sp.issubdtype(b.strain.dtype, sp.str_):
        _bstrain = _codeUTF8(b.strain)
    else:
        _bstrain = b.strain

    return ((a.chr == b.chr) &
            (a.strand == b.strand) &
            (sp.all(a.exons1 == b.exons1)) &
            (sp.all(a.exons2 == b.exons2)) &
            (sp.all(_astrain == _bstrain)) &
            (a.event_type == b.event_type) &
            (a.gene_idx == b.gene_idx) &
            (a.num_detected == b.num_detected))

def _check_files_spladder(result_dir, out_dir, prefix):
    files = []
    for event_type in ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime']:
        files.append('{}_{}_C3.confirmed.pickle'.format(prefix, event_type))
        files.append('{}_{}_C3.counts.hdf5'.format(prefix, event_type))
        files.append('{}_{}_C3.pickle'.format(prefix, event_type))
        if os.path.exists(os.path.join(result_dir, '{}_{}_C3.confirmed.gff3'.format(prefix, event_type))):
            files.append('{}_{}_C3.confirmed.gff3'.format(prefix, event_type))
        if os.path.exists(os.path.join(result_dir, '{}_{}_C3.confirmed.txt.gz'.format(prefix, event_type))):
            files.append('{}_{}_C3.confirmed.txt.gz'.format(prefix, event_type))

    for fname in files:
        _assert_files_equal(
            os.path.join(result_dir, fname),
            os.path.join(out_dir, fname)
        )

def _check_files_testing(result_dir, out_dir, suffixes):
    files = []
    for p in suffixes:
        files.extend([os.path.basename(x) for x in glob.glob(os.path.join(result_dir, p))])

    for fname in files:
        _assert_files_equal(
            os.path.join(result_dir, fname),
            os.path.join(out_dir, fname)
        )


@pytest.mark.parametrize("test_id,case", [
    ['basic', 'pos'],
    ['basic', 'neg']
])
def test_end_to_end_merge(test_id, case, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'data')
    result_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'results_merged_{}'.format(case))

    out_dir = str(tmpdir)

    my_args = ['spladder',
               'build',
               '-a', os.path.join(data_dir, 'annotation_{}.gtf'.format(case)),
               '-o', out_dir,
               '-b', ','.join([os.path.join(data_dir, 'align', '{}_{}.bam'.format(case, i+1)) for i in range(5)]),
               '--merge-strat', 'merge_graphs',
               '--extract-as',
               '-n', '15',
               '--output-conf-icgc',
               '--output-txt',
               '--output-txt-conf',
               '--output-gff3',
               '--output-struc',
               '--output-struc-conf',
               '--output-bed',
               '--output-conf-bed',
               '--output-conf-tcga',
               '-v']

    spladder.main(my_args)

    ### check that files are identical
    _check_files_spladder(result_dir, out_dir, prefix='merge_graphs')

    #### test visualization
    #my_args = ['spladder',
    #           'viz',
    #           '-o', out_dir,
    #           '-b', ':'.join([','.join([os.path.join(data_dir, 'align', '{}_{}.bam'.format(case, i+1)) for i in range(3)]),
    #                           ','.join([os.path.join(data_dir, 'align', '{}_{}.bam'.format(case, i+4)) for i in range(2)])]),
    #           '-L', 'group1,group2',
    #           '-f', 'ps',
    #           '-v']

    #spladder.main(my_args)

    #### check that files are identical
    #_assert_files_equal(os.path.join(result_dir, 'plots', 'gene_overview_C3_GENE1.ps'), os.path.join(out_dir, 'plots', 'gene_overview_C3_GENE1.ps'))

@pytest.mark.parametrize("test_id,case", [
    ['basic', 'pos'],
     ['basic', 'neg']
])
def test_end_to_end_single(test_id, case, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'data')
    result_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'results_single_{}'.format(case))

    out_dir = str(tmpdir)

    my_args = ['spladder',
               'build',
               '-a', os.path.join(data_dir, 'annotation_{}.gtf'.format(case)),
               '-o', out_dir,
               '-b', os.path.join(data_dir, 'align', '{}_1.bam'.format(case)),
               '--merge-strat', 'single',
               '--extract-ase',
               '-n', '15',
               '--output-conf-icgc',
               '-v']

    spladder.main(my_args)

    my_args = ['spladder',
               'build',
               '-a', os.path.join(data_dir, 'annotation_{}.gtf'.format(case)),
               '-o', out_dir,
               '-b', os.path.join(data_dir, 'align', '{}_1.bam'.format(case)),
               '--merge-strat', 'single',
               '--extract-ase',
               '-n', '15',
               '--output-conf-icgc',
               '--sparse-bam',
               '-v']

    spladder.main(my_args)


    ### check that files are identical
    _check_files_spladder(result_dir, out_dir, prefix='{}_1'.format(case))

@pytest.mark.parametrize("test_id", [
    'events'
])
def test_end_to_end_testing(test_id, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'data')
    result_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'results_merged')

    out_dir = str(tmpdir)

    my_args = ['spladder',
               'build',
               '-a', os.path.join(data_dir, 'testcase_{}_spladder.gtf'.format(test_id)),
               '-o', out_dir,
               '-b', ','.join([os.path.join(data_dir, 'align', 'testcase_{}_1_sample{}.bam'.format(test_id, idx+1)) for idx in range(20)]),
               '--event-types', 'exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip',
               '--merge-strat', 'merge_graphs',
               '--extract-ase',
               '--readlen', '50',
               '--output-conf-icgc',
               '--output-txt',
               '--output-txt-conf',
               '--output-gff3',
               '--output-struc',
               '--output-struc-conf',
               '--output-bed',
               '--output-conf-bed',
               '--output-conf-tcga',
               '--output-conf-icgc',
               '-v']

    spladder.main(my_args)

    ### check that event files are identical
    _check_files_spladder(result_dir, out_dir, prefix='merge_graphs')

    bamsA = os.path.join(out_dir, 'bamlistA.txt')
    bamsB = os.path.join(out_dir, 'bamlistB.txt')
    with open(bamsA, 'w') as fh:
        fh.write('\n'.join([os.path.join(data_dir, 'align', 'testcase_{}_1_sample{}.bam'.format(test_id, idx+1)) for idx in range(10)]) + '\n')
    with open(bamsB, 'w') as fh:
        fh.write('\n'.join([os.path.join(data_dir, 'align', 'testcase_{}_1_sample{}.bam'.format(test_id, idx+11)) for idx in range(10)]) + '\n')

    my_args = ['spladder',
               'test',
               '-o', out_dir,
               '-v',
               '--diagnose-plots',
               '-f', 'ps',
               '--readlen', '50',
               '--merge-strat', 'merge_graphs',
               '--event-types', 'exon_skip',
               '-a', bamsA,
               '-b', bamsB]

    spladder.main(my_args)

    ### check that files are identical
    _check_files_testing(os.path.join(result_dir, 'testing'), os.path.join(out_dir, 'testing'), suffixes=['*.tsv'])
    #_check_files_testing(os.path.join(result_dir, 'testing', 'plots'), os.path.join(out_dir, 'testing', 'plots'), suffixes=['*.ps'])


