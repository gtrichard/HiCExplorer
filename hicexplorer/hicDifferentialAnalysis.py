from __future__ import division
import argparse
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString, convertInfsToZeros_ArrayFloat
from hicmatrix.HiCMatrix import check_cooler
import numpy as np
import logging
log = logging.getLogger(__name__)

from scipy.stats import mannwhitneyu

def parse_arguments():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
Prints information about a matrix or matrices including matrix size,
number of elements, sum of elements, etc.
An example usage is:
$ hicInfo -m matrix1.h5 matrix2.h5 matrix3.h5
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrices', '-m',
                                help='The two matrices to detect differential interaction patterns on. '
                                'HiCExplorer supports the following file formats: h5 (native HiCExplorer format) '
                                'and cool.',
                                nargs=2,
                                required=True)
    parserRequired.add_argument('--outputFileName', '-o',
                                help='File names for the result of the differential analysis.',
                                required=True)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--windowSize', '-w',
                           help='Side length of the square we use to slide over the values.',
                           default=10,
                           type=int,
                           required=False)
    parserOpt.add_argument('--pvalue', '-p',
                           help='P-value acceptance level.',
                           default=0.01,
                           type=float,
                           required=False)
    parserOpt.add_argument('--qvalue', '-q',
                           help='Q-value (FDR) acceptance level.',
                           default=0.01,
                           type=float,
                           required=False)
    parserOpt.add_argument('--help', '-h', action='help', help='show this help message and exit')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def distance_count(pMatrix):
    expected_interactions = np.zeros(pMatrix.shape[0])
    row, col = pMatrix.nonzero()
    distance = np.absolute(row - col)

    for i, distance_ in enumerate(distance):
        expected_interactions[distance_] += pMatrix.data[i]
    
    return expected_interactions

def normalize_genomic_distance(pMatrix):
    genomic_distance_count = distance_count(pMatrix)
    row, col = pMatrix.nonzero()

    distance = np.ceil(np.absolute(row - col)).astype(np.int32)

    data_type = type(pMatrix.data[0])
        
    expected = genomic_distance_count[distance]
    pMatrix.data = pMatrix.data.astype(np.float64)
    pMatrix.data /= expected
    pMatrix.data = convertInfsToZeros_ArrayFloat(pMatrix.data).astype(data_type)

    return pMatrix.data

def main(args=None):

    args = parse_arguments().parse_args(args)
    
    hic_matrix_one = hm.hiCMatrix(args.matrices[0])
    hic_matrix_two = hm.hiCMatrix(args.matrices[1])

    if hic_matrix_one.getBinSize() != hic_matrix_two.getBinSize():
        log.error('Matrices do not have the same resolution')
    
    log.debug('Normalizating matrices.')

    hic_matrix_one.matrix.data = normalize_genomic_distance(hic_matrix_one.matrix)
    hic_matrix_two.matrix.data = normalize_genomic_distance(hic_matrix_two.matrix)
    log.debug('Search for differentiated regions')

    differential_regions = []
    pvalues_list = []
    i = 0
    j = 0
    while i <  (hic_matrix_one.matrix.shape[0] - args.windowSize) and i <  (hic_matrix_two.matrix.shape[0] - args.windowSize) :
        while j <  (hic_matrix_one.matrix.shape[1] - args.windowSize) and j <  (hic_matrix_two.matrix.shape[1] - args.windowSize):

            neighborhood_one = hic_matrix_one.matrix[i:i+args.windowSize, j:j+args.windowSize].toarray().flatten()
            neighborhood_two = hic_matrix_two.matrix[i:i+args.windowSize, j:j+args.windowSize].toarray().flatten()

            if np.sum(neighborhood_one) == np.sum(neighborhood_two):
                log.debug('sum identical!')
                log.debug('neighbor 1: {}'.format(neighborhood_one))
                log.debug('neighbor 2: {}'.format(neighborhood_two))

                j += args.windowSize
                continue
            test_result = mannwhitneyu(neighborhood_one, neighborhood_two)

            pvalue = test_result[1]
            if pvalue < args.pvalue:
                differential_regions.append([i, j, pvalue])
                pvalues_list.append(pvalue)
            j += args.windowSize
        i += args.windowSize
        j = 0
    log.debug('Window approach computed, correcting with FDR.')
    if len(differential_regions) == 0:
        log.info('No differential regions found!')
        exit(0)

    pvalues_list = np.sort(pvalues_list)
    largest_p_i = -np.inf
    for i, p in enumerate(pvalues_list):
        if p <= (args.qvalue * (i + 1) / len(pvalues_list)):
            if p >= largest_p_i:
                largest_p_i = p

    differential_regions_fdr_accepted = []
    for region in differential_regions:
        if region[0] < largest_p_i:
            differential_regions_fdr_accepted.append(region)
    

    # map regions to genomic positions
    mapped_regions = []
    for region in differential_regions_fdr_accepted:
        chr_x, start_x, end_x, _ = hic_matrix_one.getBinPos(region[0])
        chr_y, start_y, end_y, _ = hic_matrix_one.getBinPos(region[1])
        window_size = hic_matrix_one.getBinPos(region[0])[1] - hic_matrix_one.getBinPos(region[0]+args.windowSize)[1]
        mapped_regions.append((chr_x, start_x, end_x, chr_y, start_y, end_y, window_size))
    
    with open(args.outputFileName, 'w') as fh:
       for region in mapped_regions:
                fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % region)

    hic_matrix_one.save('normalized.'+args.matrices[0])
    hic_matrix_two.save('normalized.'+args.matrices[1])
