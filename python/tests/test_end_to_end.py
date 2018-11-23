import gzip
import h5py
import os
import scipy as sp

import pytest


data_dir = os.path.join(os.path.dirname(__file__), 'data')
sample_dir = os.path.join(os.path.dirname(__file__), 'testcase1')
from spladder import spladder


def _compare_hdf5(expected, actual):
    for k in expected.keys():
        if isinstance(expected[k], h5py._hl.group.Group):
            for l in expected[k].keys():
                assert sp.all(expected[k][l][:] == actual[k][l][:])
        else:
            try:
                assert sp.all(expected[k][:] == actual[k][:])
            except:
                import pdb
                pdb.set_trace()


def _assert_files_equal(expected_path, actual_path):
    def o(f):
        if f.endswith('.gz'):
            return gzip.open(f, 'rb')
        elif f.endswith('.hdf5'):
            return h5py.File(f, 'r')
        else:
            return open(f, 'r')

    with o(expected_path) as e:
        with o(actual_path) as a:
            if expected_path.endswith('.hdf5'):
                _compare_hdf5(e, a)
            else:
                assert e.read() == a.read()


def _check_files(event_type, result_dir, out_dir):
    files = []
    for event_type in ['exon_skip', 'intron_retention', 'alt_3prime', 'alt_5prime']:
        files.append('merge_graphs_{}_C3.confirmed.pickle'.format(event_type))
        files.append('merge_graphs_{}_C3.counts.hdf5'.format(event_type))
        files.append('merge_graphs_{}_C3.pickle'.format(event_type))
        if os.path.exists(os.path.join(result_dir, 'merge_graphs_{}_C3.confirmed.gff3'.format(event_type))):
            files.append('merge_graphs_{}_C3.confirmed.gff3'.format(event_type))
        if os.path.exists(os.path.join(result_dir, 'merge_graphs_{}_C3.confirmed.txt.gz'.format(event_type))):
            files.append('merge_graphs_{}_C3.confirmed.txt.gz'.format(event_type))

    for fname in files:
        _assert_files_equal(
            os.path.join(result_dir, fname),
            os.path.join(out_dir, fname)
        )

@pytest.mark.parametrize("test_id,case", [
    ['new', 'pos'],
    ['new', 'neg']
])


def test_end_to_end_merge(test_id, case, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'data')
    result_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'results_merged_{}'.format(case))

    out_dir = str(tmpdir)

    my_args = ['spladder',
               '-a', os.path.join(data_dir, 'annotation_{}.gtf'.format(case)),
               '-o', out_dir,
               '-b', ','.join([os.path.join(data_dir, 'align', '{}_{}.bam'.format(case, i+1)) for i in range(5)]),
               '--merge_strat', 'merge_graphs',
               '-T', 'y',
               '-n', '15']

    spladder.main(my_args)

    ### check that files are identical
    _check_files(event_type, result_dir, out_dir) 

def test_end_to_end_single(test_id, case, tmpdir):
    data_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'data')
    result_dir = os.path.join(os.path.dirname(__file__), 'testcase_{}'.format(test_id), 'results_single_{}'.format(case))

    out_dir = str(tmpdir)

    my_args = ['spladder',
               '-a', os.path.join(data_dir, 'annotation_{}.gtf'.format(case)),
               '-o', out_dir,
               '-b', ','.join([os.path.join(data_dir, 'align', '{}_{}.bam'.format(case, i+1)) for i in range(5)]),
               '--merge_strat', 'merge_graphs',
               '-T', 'y',
               '-n', '15']

    spladder.main(my_args)

    ### check that files are identical
    _check_files(event_type, result_dir, out_dir) 

