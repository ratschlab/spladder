import gzip
import h5py
import os
import scipy as sp

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
            assert sp.all(expected[k][:] == actual[k][:])


def _assert_files_equal(expected_path, actual_path):
    def o(f):
        if f.endswith('.gz'):
            return gzip.open(f, 'rb')
        elif f.endswith('.hdf5'):
            return h5py.File(f, 'r')
        elif f.endswith('.pickle'):
            return open(f, 'rb')
        else:
            return open(f, 'r')

    with o(expected_path) as e:
        with o(actual_path) as a:
            if expected_path.endswith('.hdf5'):
                _compare_hdf5(e, a)
            elif expected_path.endswith('.pickle'):
                ta = pickle.load(a, encoding='latin1')
                te = pickle.load(e, encoding='latin1')
                if len(ta) == 0 or isinstance(ta[0], sp.int64):
                    assert sp.all(ta == te)
                else:
                    assert sp.all([_compare_gene(_[0], _[1]) for _ in zip(ta, te)])
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
            (a.id == b.id) &
            (a.num_detected == b.num_detected))

def _check_files(result_dir, out_dir, prefix):
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

@pytest.mark.parametrize("test_id,case", [
    ['basic', 'pos'],
    ['basic', 'neg']
])
def test_end_to_end_merge(test_id, case, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'data')
    result_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'results_merged_{}'.format(case))

    out_dir = str(tmpdir)

    my_args = ['spladder',
               '-a', os.path.join(data_dir, 'annotation_{}.gtf'.format(case)),
               '-o', out_dir,
               '-b', ','.join([os.path.join(data_dir, 'align', '{}_{}.bam'.format(case, i+1)) for i in range(5)]),
               '--merge-strat', 'merge_graphs',
               '--extract-as',
               '-n', '15']

    spladder.main(my_args)

    ### check that files are identical
    _check_files(result_dir, out_dir, prefix='merge_graphs') 

@pytest.mark.parametrize("test_id,case", [
    ['basic', 'pos'],
    ['basic', 'neg']
])
def test_end_to_end_single(test_id, case, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'data')
    result_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'results_single_{}'.format(case))

    out_dir = str(tmpdir)

    my_args = ['spladder',
               '-a', os.path.join(data_dir, 'annotation_{}.gtf'.format(case)),
               '-o', out_dir,
               '-b', os.path.join(data_dir, 'align', '{}_1.bam'.format(case)),
               '--merge-strat', 'single',
               '--extract-as',
               '-n', '15']
               #'-b', ','.join([os.path.join(data_dir, 'align', '{}_{}.bam'.format(case, i+1)) for i in range(5)]),

    spladder.main(my_args)

    ### check that files are identical
    _check_files(result_dir, out_dir, prefix='{}_1'.format(case)) 

