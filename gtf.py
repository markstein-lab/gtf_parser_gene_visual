import gene_draw
import numpy as np
"""
Read GFF/GTF files. Works with gzip compressed files and pandas.
    http://useast.ensembl.org/info/website/upload/gff.html
"""


from collections import defaultdict
import gzip
import pandas as pd
import re


GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')


def dataframe(filename):
    """Open an optionally gzipped GTF file and return a pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)
    # breakpoint()
    for i, line in enumerate(lines(filename)):
        # breakpoint()
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                # breakpoint()
                result[key] = [None] * i
        # breakpoint()
        # Ensure this row has some value for each column.
        for key in result.keys():
            # breakpoint()
            result[key].append(line.get(key, None))
    send = pd.DataFrame(result)
    # breakpoint()
    return send


def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open
    # breakpoint()
    with fn_open(filename) as fh:
        for line in fh:
            # breakpoint()
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')
    # breakpoint()
    for i, col in enumerate(GTF_HEADER):
        # breakpoint()
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]
    '''
    infos[7] = 'transcript_name "ABCB7-RA"' -> line 2 in abcb7.gtf
    '''
    # breakpoint()
    for i, info in enumerate(infos, 1):
        # breakpoint()
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)

    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value

def get_exon_ra_rb(transcript_name):
    ra_or_rb = transcript_name.split('-')[1]
    return ra_or_rb

if __name__ == '__main__':
    result = dataframe("abcb7.gtf")
    # breakpoint()
    # print(result)
    '''
    result.transcript_name
    '''
    gene_start = result.start[0]
    gene_end = result.end[0]
    # for fields in result.iteritems():
    ra_exon = []
    rb_exon = []
    exon_pos =[]
    for i, (f,t_n) in enumerate(zip(result.feature, result.transcript_name)):
        # breakpoint()
        if f == 'exon':
            to_add = [int(result.start[i]),int(result.end[i])]
            exon_pos.append(to_add)
            # -----------------
            if get_exon_ra_rb(t_n) == 'RA':
                ra_exon.append((t_n, to_add))
            else:
                rb_exon.append((t_n, to_add))
            # ---------------
            
    # breakpoint()
    # ra_exon = tuple(ra_exon)
    # rb_exon = tuple(rb_exon)
    exon_pos = tuple(exon_pos)
    
    exon_plus_ra_rb = [exon_pos, ra_exon, rb_exon]

    geneLength = [int(gene_start),int(gene_end)]

    # np.savetxt("parsed_data.txt", result.values, fmt='%s')

    gene = gene_draw.GeneImage(geneLength,exon_plus_ra_rb)
    gene.show()