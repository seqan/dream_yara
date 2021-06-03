#!/usr/bin/env python2
"""Execute the tests for Yara.

The golden test outputs are generated by the script generate_outputs.sh.

You have to give the root paths to the source and the binaries as arguments to
the program.  These are the paths to the directory that contains the 'projects'
directory.

Usage:  run_tests.py SOURCE_ROOT_PATH BINARY_ROOT_PATH
"""
import logging
import os.path
import sys
import glob

# Automagically add util/py_lib to PYTHONPATH environment variable.
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..',
                                    'include', 'seqan', 'util', 'py_lib'))
sys.path.insert(0, path)

import seqan.app_tests as app_tests

log_transforms = [
	app_tests.RegexpReplaceTransform("[0-9\.\-e]+ sec", "0.0 sec"),
	app_tests.RegexpReplaceTransform("Free [0-9]+ of [0-9]+ MB", "Free 0 of 0 MB")
]

sam_transforms = [
	app_tests.RegexpReplaceTransform("@PG.*", "@PG")
]

def main(source_base, binary_base):
    """Main entry point of the script."""

    # gold standard binary files created on little endian
    if sys.byteorder != 'little':
        print 'Skipping tests for Yara on big endian'
        print '====================================='
        return 0

    print 'Executing test for Yara'
    print '=============================='
    print

    ph = app_tests.TestPathHelper(
        source_base, binary_base,
        'tests')  # tests dir

    # ============================================================
    # Auto-detect the binary path.
    # ============================================================

    path_to_indexer = app_tests.autolocateBinary(
      binary_base, 'bin', 'yara_indexer')

    path_to_mapper = app_tests.autolocateBinary(
      binary_base, 'bin', 'yara_mapper')

    # ============================================================
    # Built TestConf list.
    # ============================================================

    # Build list with TestConf objects, analoguely to how the output
    # was generated in generate_outputs.sh.
    conf_list = []

    # ============================================================
    # Run Indexer Tests
    # ============================================================

    for organism in ['adeno']:

        # Get file extensions for the fm index files
        exts = [os.path.basename(fname).split('.', 1)[-1]
                for fname in glob.glob(ph.inFile('gold/%s-genome.*' % organism))]

        conf = app_tests.TestConf(
            program=path_to_indexer,
            args=[ph.inFile('input/%s-genome.fa' % organism),
                  '-o', ph.outFile('%s-genome' % organism)],
            to_diff=[(ph.inFile('gold/%s-genome.%s' % (organism, ext)),
                     ph.outFile('%s-genome.%s' % (organism, ext)), 'md5') for ext in exts])
        conf_list.append(conf)

    # ============================================================
    # Run Single-End Mapper Tests
    # ============================================================

#    mapper_args = [
#                    ['--threads', '1' ]],
#                    ['--threads', '8' ]
#                  ]
#    mapper_suffix = ['t1', 't8']

    mapper_args = [['--threads', '1']]
    mapper_suffix = ['t1']

    for organism in ['adeno']:
        for i in range(0, len(mapper_args)):

            conf = app_tests.TestConf(
                program=path_to_mapper,
                args=[ph.inFile('gold/%s-genome' % organism),
                      ph.inFile('input/%s-reads_1.fa' % organism),
                      '-o', ph.outFile('%s-reads_1.%s.sam' % (organism, mapper_suffix[i]))] +
                      mapper_args[i],
                to_diff=[(ph.inFile('gold/%s-reads_1.%s.sam' % (organism, mapper_suffix[i])),
                          ph.outFile('%s-reads_1.%s.sam' % (organism, mapper_suffix[i])),
                          sam_transforms)])
            conf_list.append(conf)

    # ============================================================
    # Execute the tests
    # ============================================================

    failures = 0
    for conf in conf_list:
        res = app_tests.runTest(conf)
        # Output to the user.
        print ' '.join([conf.program] + conf.args),
        if res:
             print 'OK'
        else:
            failures += 1
            print 'FAILED'

    # Cleanup.
    ph.deleteTempDir()

    print '=============================='
    print '     total tests: %d' % len(conf_list)
    print '    failed tests: %d' % failures
    print 'successful tests: %d' % (len(conf_list) - failures)
    print '=============================='
    # Compute and return return code.
    return failures != 0


if __name__ == '__main__':
    sys.exit(app_tests.main(main))
