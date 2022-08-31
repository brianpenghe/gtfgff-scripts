##################################
#                                #
# Last modified 04/05/2014       #
#                                #
# Georgi Marinov                 #
#                                #
##################################
from __future__ import print_function

from argparse import ArgumentParser
import numpy
import logging
import os

logger = logging.getLogger('gene_coverage_wig_gtf')

try:
    import pyBigWig
except ImportError:
    logger.warn('pyBigWig is not available, no direct bigwig support')
    pyBigWig = None

try:
    import pysam
except ImportError:
    logger.warn('pysam is not available, no direct BAM support')
    pysam = None

version = '2.0'

NORMALIZATIONS = {
    'sum': numpy.mean,
    'max': numpy.max,
    'mean': numpy.mean,
    'median': numpy.median,
}


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.quiet:
        logging.basicConfig(level=logging.WARN)
    else:
        logging.basicConfig(level=logging.INFO)

    if args.version:
        parser.exit(message="version {}\n".format(version))

    if args.source_type is not None:
        logger.info('will only consider genes from source %s',
                    args.source_type)

    if args.max_gene_length is not None:
        logger.info(
            'will only consider genes longer than %s and shorter than %s',
            args.min_gene_length,
            args.max_gene_length)
    else:
        logger.info('will only consider genes longer than %s',
                    args.min_gene_length)

    if args.normalize:
        logger.info('will normalize scores')

    if args.gene_normalization in NORMALIZATIONS:
        logger.info('Will normalize gene bins by %s',
                    args.gene_normalization)

    if args.all_gene_models:
        logger.info('will only all genes')
    else:
        logger.info('will only use genes with one isoform')

    if args.gene_type:
        logger.info('will only consider genes of type %s', args.gene_type)

    geneDict = loadAnnotation(args.gtf, args.source_type, args.gene_type)
    logger.info('genes passed type filters %s', len(geneDict))

    gene_coverage = loadGeneCoverage(args.filename, geneDict,
                                     args.all_gene_models)

    if args.print_list:
        geneListFilename = args.output + '.geneList'
    else:
        geneListFilename = None

    outputArray = createCoveragePercentiles(
        geneDict, gene_coverage,
        args.min_gene_length, args.max_gene_length,
        geneListFilename,
        args.normalize,
        args.gene_normalization,
        args.expression_threshold,
    )

    with open(args.output, 'wt') as outfile:
        for i in range(100):
            outfile.write(str(i) + '\t' + str(outputArray[i])+'\n')


def make_parser():
    parser = ArgumentParser()
    parser.add_argument('--gtf', required=True, help='GTF file name')
    parser.add_argument(
        'filename',
        help='file of type bedgraph or bigwig to compute coverage of')
    parser.add_argument('-o', '--output', required=True,
                        help='output file name')
    parser.add_argument('--source-type', help='source type name')
    parser.add_argument('--max-gene-length', default=None, type=int,
                        help='maximum gene length to consider')
    parser.add_argument('--min-gene-length', type=int, default=1000,
                        help='minimum gene length to consider')
    parser.add_argument(
        '--all-gene-models', default=False, action='store_true',
        help='use all gene models, not just the ones with a single model')
    parser.add_argument('--gene-type', help='limit to specified gene type')
    parser.add_argument('--print-list', default=False, action='store_true',
                        help='write gene ids considered to <outfile>.geneList')
    parser.add_argument('--normalize', action='store_true', default=False,
                        help='normalize scores by number of genes')
    parser.add_argument('--gene-normalization',
                        choices=['none'] + sorted(NORMALIZATIONS.keys()),
                        default='none',
                        help='Per gene model normalization. sum is total '
                             'number of reads assigned to gene. '
                             'max is the maximum gene bin size.')
    parser.add_argument('--expression-threshold', default=0.0, type=float,
                        help='at least one bin must be >= this threshold to be included')

    parser.add_argument('-q', '--quiet', action='store_true',
                        help='only report errors')
    parser.add_argument('--version', action='store_true',
                        help='report version number')
    return parser


def loadAnnotation(filename, source_type, gene_type_filter):
    with open(filename, 'rt') as gtfStream:
        return parseAnnotation(
            gtfStream,
            source_type,
            gene_type_filter)


def parseAnnotation(stream, source_type, gene_type_filter):
    geneDict = {}
    for line in stream:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        if source_type is not None and fields[1] != source_type:
            continue
        chromosome = fields[0]
        strand = fields[6]
        left = int(fields[3])
        right = int(fields[4])
        if gene_type_filter is not None:
            try:
                geneType = getGFFAttributeValueByKey(fields[8], 'gene_type')
            except ValueError:
                continue
            if geneType != gene_type_filter:
                continue
        geneID = getGFFAttributeValueByKey(fields[8], 'gene_id')
        transcriptID = getGFFAttributeValueByKey(fields[8], 'transcript_id')
        if transcriptID in geneDict.setdefault(geneID, {}):
            if chromosome != geneDict[geneID][transcriptID][0][0]:
                continue
        else:
            geneDict[geneID][transcriptID] = []
        geneDict[geneID][transcriptID].append(
            (chromosome, left, right, strand))

    logger.info('finished inputting annotation %s', len(geneDict.keys()))
    return geneDict


def getGFFAttributeValueByKey(field, name):
    start = field.find(name)
    if start == -1:
        return None
    start += len(name)
    start = field.find('"', start)
    if start == -1:
        return None
    start += 1
    end = field.index('"', start)
    return field[start:end]


def loadGeneCoverage(filename, geneDict, all_gene_models):
    """Read coverage data for genes of interest

    :Parameters:
      - filename: source file to read
      - geneDict: GTF annotations
      - all_gene_models: flag if we should use all gene models instead
        of just single gene models
    """
    _, ext = os.path.splitext(filename)
    if ext in ('.bam', '.sam'):
        gene_coverage = readSAM(filename, geneDict, all_gene_models)
        return gene_coverage

    instream = guessFileOpen(filename)
    if instream is None:
        logger.error('Unable to open %s. Not a supported file-type',
                     os.path.abspath(filename))
        raise RuntimeError('Unsupported file type')

    if pyBigWig and \
       isinstance(instream, pyBigWig.pyBigWig) and \
       instream.isBigWig():
        gene_coverage = readBigwig(instream, geneDict, all_gene_models)
    else:
        gene_coverage = readWiggle(instream, geneDict, all_gene_models)

    return gene_coverage


def guessFileOpen(filename):
    """Try opening the file with supported file readers

    returns the opened object
    """
    stream = open(filename, 'rt')
    header = stream.read(100)

    is_text = True
    for c in header:
        if not (c.isprintable() or c.isspace()):
            is_text = False
            break
    stream.seek(0)
    if is_text:
        return stream

    if pysam:
        try:
            stream = pysam.AlignmentFile(filename, "r")
            return stream
        except ValueError:
            pass

    if pyBigWig:
        try:
            stream = pyBigWig.open(filename)
            return stream
        except RuntimeError as e:
            logger.debug('%s is not a bigwig file. %s', filename, str(e))

    if not os.path.exists(filename):
        logger.error('%s is not a local file', filename)

    return None


def readSAM(filename, geneDict, all_gene_models):
    coverageDict = initializeCoverageDict(geneDict, all_gene_models)
    to_pm = calculate_bam_normalization(filename)
    all_seen = 0
    reads_seen = 0
    with pysam.AlignmentFile(filename) as aligned:
        for chromosome in coverageDict:
            if chromosome.startswith("chr"):
                for read in aligned.fetch(chromosome):
                    all_seen += 1
                    if not read.is_unmapped and read.get_tag('NH') == 1:
                        reads_seen += 1
                        for start, stop in read.get_blocks():
                            for pos in range(start, stop):
                                if pos in coverageDict[chromosome]:
                                    coverageDict[chromosome][pos] += 1

    #dump_wiggle(coverageDict)
    print("Reads seen: {} / {}. ".format(reads_seen, all_seen))
    to_pm = 1e6 / reads_seen
    for chromosome in coverageDict:
        for pos in sorted(coverageDict[chromosome]):
            coverageDict[chromosome][pos] = numpy.round(coverageDict[chromosome][pos] * to_pm, 5)

    logger.info('finished inputting bigwig')
    logger.info('genes passed type filters %s', len(geneDict))
    logger.info('Normalizing {} reads to per million factor {}'.format(reads_seen, to_pm))
    return coverageDict


def dump_wiggle(coverageDict):
    print('dw', coverageDict.keys())
    with open("temp.wig", "wt") as output:
        for chromosome in coverageDict:
            last_score = None
            last_start = None
            last_end = None
            for pos in sorted(coverageDict[chromosome]):
                score = coverageDict[chromosome][pos]
                if last_score is None:
                    last_score = score
                    last_start = pos
                elif last_end < pos or score != last_score:
                    output.write("{}\t{}\t{}\t{}\n".format(
                        chromosome,
                        last_start,
                        last_end,
                        numpy.round(last_score, 5)))
                    last_score = score
                    last_start = pos
                last_end = pos + 1
            output.write("{}\t{}\t{}\t{}\n".format(
                chromosome, last_start, last_end, numpy.round(score, 5)))


def calculate_bam_normalization(filename, prefix="chr"):
    with pysam.AlignmentFile(filename) as aligned:
        uniq = 0
        for read in aligned:
            if not read.is_unmapped and read.reference_name.startswith(prefix):
                if read.get_tag('NH') == 1:
                    uniq += 1
        to_pm = 1e6 / uniq
        logger.info('Normalizing {} uniq reads to per million factor {}'.format(uniq, to_pm))
    return to_pm


def readBigwig(bigwig, geneDict, all_gene_models):
    coverageDict = initializeCoverageDict(geneDict, all_gene_models)
    chromosomes = set(coverageDict).intersection(set(bigwig.chroms()))
    for chromosome in chromosomes:
        for left, right, score in bigwig.intervals(chromosome):
            for j in range(left, right):
                if j in coverageDict[chromosome]:
                    coverageDict[chromosome][j] = score

    logger.info('finished inputting bigwig')
    logger.info('genes passed type filters %s', len(geneDict))
    return coverageDict


def readWiggle(wiggle, geneDict, all_gene_models):
    coverageDict = initializeCoverageDict(geneDict, all_gene_models)
    for line in wiggle:
        if line.startswith('track'):
            continue
        if line.startswith('#'):
            continue
        fields = line.replace(' ', '\t').strip().split('\t')
        chromosome = fields[0]
        left = int(fields[1])
        right = int(fields[2])
        score = float(fields[3])
        if chromosome not in coverageDict:
            continue
        for j in range(left, right):
            if j in coverageDict[chromosome]:
                coverageDict[chromosome][j] = score

    logger.info('finished inputting wiggle')
    return coverageDict


def initializeCoverageDict(GeneDict, all_gene_models):
    CoverageDict = {}
    genesToRemove = set()
    for geneID in GeneDict:
        if all_gene_models is False and len(GeneDict[geneID]) > 1:
            genesToRemove.add(geneID)
            continue
        for transcriptID in GeneDict[geneID]:
            transcript = GeneDict[geneID][transcriptID]
            for (chromosome, left, right, strand) in transcript:
                for j in range(left, right):
                    CoverageDict.setdefault(chromosome, {})[j] = 0
    for geneID in genesToRemove:
        del GeneDict[geneID]
    logger.info('Removed %s multi-model genes', len(genesToRemove))
    return CoverageDict


def createCoveragePercentiles(
        GeneDict, coverageDict,
        minGeneLength, maxGeneLength=None,
        geneListFilename=None,
        doNormalize=False,
        gene_normalization=None,
        expression_threshold=0.0):
    outputArray = numpy.zeros(shape=100)

    geneListStream = open(geneListFilename, 'wt') if geneListFilename else None

    geneNumber = 0
    for geneID in GeneDict:
        NucleotideList, chromosome = buildNucleotideList(GeneDict[geneID])
        geneLength = len(NucleotideList)
        if geneLength < minGeneLength:
            continue
        if maxGeneLength is not None and geneLength > maxGeneLength:
            continue
        final_vector = numpy.zeros(shape=100)
        bins = numpy.linspace(0, geneLength, num=101, dtype=int)
        start = bins[0]
        for i, end in enumerate(bins[1:]):
            region = NucleotideList[start:end]
            counts = [coverageDict[chromosome][pos] for pos in region]
            final_vector[i] = numpy.mean(counts)
            start = end
        assert len(final_vector) == 100
        final_vector_sum = numpy.sum(final_vector)
        if final_vector_sum > 0 and numpy.max(final_vector) > expression_threshold:
            if geneListStream:
                geneListStream.write(geneID)
                geneListStream.write('\t')
                geneListStream.write('\t'.join([str(x) for x in final_vector]))
                geneListStream.write('\n')

            if gene_normalization in NORMALIZATIONS:
                normalization_func = NORMALIZATIONS[gene_normalization]
                final_vector /= normalization_func(final_vector)

            outputArray += final_vector
            geneNumber += 1

    logger.info('%s genes considered', geneNumber)

    if doNormalize:
        outputArray /= float(geneNumber)

    if geneListStream:
        geneListStream.close()

    return outputArray


def buildNucleotideList(geneModel):
    NucleotideList = []
    for transcriptID in geneModel:
        for (chromosome, left, right, strand) in geneModel[transcriptID]:
            for i in range(left, right):
                NucleotideList.append(i)
    NucleotideList = sorted(set(NucleotideList))
    if strand == '-' or strand == 'R':
        NucleotideList.reverse()
    return NucleotideList, chromosome


if __name__ == '__main__':
    main()
