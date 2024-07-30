"""Parses and aligns sequences of human and swine H1 HAs descended from 1918.

Requires the files listing anomalous sequences: 'phylogenetics/JVI_82_8947_Anomalies.txt' and 'phylogenetics/JDB_Anomalies.txt'.
  
These functions utilize the `EMBOSS needle`_ program for the alignments,
and requires the ``mapmuts`` program. Originally written by Jesse Bloom
and further modified by Mike Doud, February 2016.
"""
import sys
import os
import re
import random
import datetime
import mapmuts.align
import mapmuts.sequtils

def DateToOrdinal(datestring, refyear=1968):
    """Converts a date string to an ordinal date.

    *datestring* is a date given by a string such as '2007/2/13' (for
    Feb-13-2007), or '2007/2//' if no day is specified, or
    '2007//' if no day or month is specified. The '/' characters can
    also be '-'.

    *refdate* is an integer year from the approximate timeframe we are examining
    which is used to anchor the datestring date on the assumption
    that each year has 365.25 days.

    The returned value is a number (decimal) giving the date. If no
    day is specified, the 15th (halfway through the month) is chosen.
    If no month or day is specified, July 1 (halfway through the
    year) is chosen.

    >>> print "%.2f" % DateToOrdinal('2007/4/27', 1968)
    2007.32

    >>> print "%.2f" % DateToOrdinal('2007/4/', 1968)
    2007.29

    >>> print "%.2f" % DateToOrdinal('2007//', 1968)
    2007.50

    >>> print "%.2f" % DateToOrdinal('2007-4-27', 1968)
    2007.32

    """
    if not isinstance(refyear, int):
        raise ValueError('refyear is not an integer')
    refdate = datetime.date(refyear, 1, 1).toordinal()
    try:
        if '/' in datestring:
            (year, month, day) = datestring.split('/')
        else:
            (year, month, day) = datestring.split('-')
    except ValueError:
        raise ValueError("Invalid datestring of: %s" % str(datestring))
    if year and month and day:
        (year, month, day) = (int(year), int(month), int(day))
        date = datetime.date(year, month, day)
    elif year and month:
        (year, month) = (int(year), int(month))
        date = datetime.date(year, month, 15)
    elif year:
        year = int(year)
        date = datetime.date(year, 7, 1)
    else:
        raise ValueError("Invalid datestring of: %s" % str(datestring))
    return (date.toordinal() - refdate) / 365.25 + refyear


def ParseHeader(head):
    """Parses an influenza sequence header for year, subtype, and host.

    *head* is a string specifying a sequence header, in a format such as
    these examples::

        cds:CAA24268 Human A/Puerto Rico/8/1934 1934// H1N1 Puerto Rico NP
        cds:BAN57584 Avian A/muscovy duck/Vietnam/LBM330/2013 2013/01/07 H5N1 Viet Nam NP
        cds:AGN51218 Swine A/swine/Illinois/A01158960/2013 2013/02/08 H1N1 USA NP

    This function parses the header and returns a 4-tuple 
    *(strain, host, date, subtype)*. The date is extracted using
    *DateToOrdinal*.

    For example, the three example header above would return (if dates
    are rounded to one decimal for display):

        * *('A Puerto Rico/8/1934', 'Human', 1934.5, 'H1N1')*
        
        * *('A/muscovy duck/Vietnam/LBM330/2013', 'Avian', 2013.0, 'H5N1')*
        
        * *('A/swine/Illinois/A01158960/2013', 'Swine', 2013.1, 'H1N1')*

    If the header cannot be parsed, returns *None*. This happens when:

        * The header contains the string *mixed*.

        * No 4-digit year is specified (i.e. the date is less than 1800)

        * The subtype is not fully known

        * No host is specified

    """
    if 'mixed' in head:
        return
    headmatch = re.compile('^\S+ (?P<host>[\w ]*) (?P<strain>(A|swine)(\/[\w\. \-\'\:\(\)\?\+]*){1,6}) +(?P<date>(\d*|[Uu]nknown)\/\d*\/\d*) (?P<subtype>\w+)')
    m = headmatch.search(head)
    if not m:
        raise ValueError("Failed to match header:\n%s" % head)
    host = m.group('host')
    if not host:
        return
    date = m.group('date')
    if 'unknown' in date.lower() or date == '//':
        return
    date = DateToOrdinal(date)
    if date < 1800:
        return
    subtypematch = re.compile('^H\d+N\d+$')
    subtype = m.group('subtype')
    if not subtypematch.search(subtype):
        return
    return (m.group('strain'), host, date, subtype)


def GetUnique(seqs):
    """Gets unique sequences from a set.

    *seqs* is a list of sequences as *(header, sequence)* 2-tuples.

    This method returns a new list *uniqueseqs* which is like *seqs*
    except that non-unique sequences have been removed. Specifically,
    all sequences that are identical to another sequence are removed.
    The first sequence is retained.
    """
    seq_d = dict([(seq, True) for (head, seq) in seqs])
    uniqueseqs = []
    for (head, seq) in seqs:
        if seq_d[seq]:
            uniqueseqs.append((head, seq))
            seq_d[seq] = False
    return uniqueseqs


def NeedleCDSandProtAlignments(refseq, seqs, needlecmd, tempfile='_alignments.temp'):
    """Uses EMBOSS needle to align sequences in *seqs* to *refseq*.

    The sequences in *seqs* and *refseq* should be coding DNA sequences.
    They must be translatable, but ambiguous nucleotides and truncation
    of incomplete sequences are allowed.

    *refseq* is a string giving the reference sequence.

    *seqs* is a list of *(header, sequence)* 2-tuples.

    *needlecmd* is the path to the EMBOSS needle program.

    *tempfile* is the name of a file temporarily created to store
    the alignments, and then deleted.

    Returns the 2-tuple *(prot_alignment, cds_alignment)*.

    *prot_alignment* contains each of the proteins encoded in 
    *seqs* as a *(header, protsequence)* 2-tuple, with the
    protein sequence aligned to that in *refseq* and with all
    gaps relative to *resfeq* stripped away.

    *cds_alignment* is an alignment of the coding DNA sequences
    in *prot_alignment*, with the nucleotide alignments done
    according to the protein alignments.
    """
    prots = []
    heads = []
    refprot = mapmuts.sequtils.Translate([('head', refseq)])[0][1]
    for (head, seq) in seqs:
        try:
            (head, prot) = mapmuts.sequtils.Translate([(head, seq)], readthrough_n=True, truncate_incomplete=True)[0]
            heads.append(head)
            prots.append(prot)
        except:
            sys.stderr.write("PROBLEM translating sequence %s" % head)
            raise
    try:
        mapmuts.align.Needle(refprot, prots, needlecmd, 'protein', tempfile)
        alignments = mapmuts.sequtils.ReadFASTA(tempfile)
    finally:
        if os.path.isfile(tempfile):
            os.remove(tempfile)
    assert len(alignments) == 2 * len(prots) == 2 * len(heads) == 2 * len(seqs)
    prot_alignment = []
    cds_alignment = []
    for i in range(len(prots)):
        prot = prots[i]
        head = heads[i]
        seq = seqs[i][1]
        assert seqs[i][0] == head
        (refa, prota) = (alignments[2 * i][1], alignments[2 * i + 1][1])
        assert len(refa) == len(prota)
        iref = iprot = 0
        alignedprot = []
        alignedcds = []
        for (aa_ref, aa_prot) in zip(refa, prota):
            assert (aa_ref == '-' or aa_ref == refprot[iref])
            assert (aa_prot == '-' or aa_prot == prot[iprot])
            if aa_ref == '-' and aa_prot != '-':
                iprot += 1
            elif aa_prot == '-' and aa_ref != '-':
                alignedprot.append(aa_prot)
                alignedcds.append('---')
                iref += 1
            elif aa_ref != '-' and aa_prot != '-':
                alignedprot.append(aa_prot)
                alignedcds.append(seq[3 * iprot : 3 * iprot + 3])
                iref += 1
                iprot += 1
            else:
                raise ValueError("Both prots in alignment have gap")
        alignedprot = ''.join(alignedprot)
        alignedcds = ''.join(alignedcds)
        assert alignedprot == mapmuts.sequtils.Translate([(head, alignedcds)], readthrough_n=True, truncate_incomplete=True, translate_gaps=True)[0][1]
        assert len(alignedprot) == len(refprot)
        prot_alignment.append((head, alignedprot))
        cds_alignment.append((head, alignedcds))
    assert len(prot_alignment) == len(cds_alignment)
    return (prot_alignment, cds_alignment)


def MakeAlignment(infile, nperyearhostsubtype, refseqfile, outfile, base_directory):
    random.seed(1)

    # input / output files and script parameters
    # For further downsampling:
    only_keep_every_N_years = False
    every_N_yr = 2 # only meaningful if only keeping every N years...
    if not only_keep_every_N_years:
        every_N_yr = 1

    minlength = 1695 # only keep sequences of at least this length
    purgeambiguous = True # remove sequences with ambiguous nucleotides
    needlecmd = '/app/emboss/6.6.0/bin/needle' # path to EMBOSS needle    
    
    anomalies = frozenset([line.strip() for line in open('%s/phylogenetics/JVI_82_8947_Anomalies.txt' % base_directory).readlines() + open('%s/phylogenetics/JDB_Anomalies.txt' % base_directory).readlines() if not line.isspace()]) # strain names of anomalous sequences
    subtype = 'H1N\d{1,2}' # only keep sequences with subtypes matching this regular expression
    hosts_subtypes_years = {
        # For each host, valus are dictionaries keyed by 
        # matching subtypes. Values are dictionaries keyed by years to
        # retain, with keys the number of years to subtract off that
        # year for date stamping. 
        # 24 years are subtracted of human seasonal H1N1 to conform to clock.
        'Human':{
                'H1N1':dict([(year, 0) for year in range(1918, 1958)] + [(year, 24) for year in range(1977, 2009)] + [(year, 0) for year in range(2009, 2015)]),
                },
        'Swine':{
                'H1N1':dict([(year, 0) for year in range(1918, 2015)]),
                },
            }

    # get reference sequence
    if not os.path.isfile(refseqfile):
        raise IOError("Could not find refseqfile of %s" % refseqfile)
    refseq = mapmuts.sequtils.ReadFASTA(refseqfile)
    if len(refseq) != 1:
        raise IOError("Failed to read exactly one sequence from refseqfile %s" % refseqfile)
    refprot = mapmuts.sequtils.Translate(refseq)[0][1]
    (refhead, refseq) = refseq[0]
    print "\nUsing a reference sequence of %s.\nThis sequence is %d nucleotides in length, and encodes a protein of %d residues." % (refhead, len(refseq), len(refprot))

    # make alignment
    if not os.path.isfile(infile):
        raise IOError("Could not find infile %s" % infile)
    seqs = mapmuts.sequtils.ReadFASTA(infile)
    random.shuffle(seqs) # so they are not in order in case the download ordered them in some way. This is important since some of the purging of redundant sequences is done in order of appearance.
    print "\nRead %d sequences from %s." % (len(seqs), infile)
    print "\nPurging pre-defined anomalous sequences..."
    cleanseqs = []
    for (head, seq) in seqs:
        for anomaly in anomalies:
            if anomaly in head:
                break
        else:
            cleanseqs.append((head, seq))
    seqs = cleanseqs
    print "Retained %d sequences after removing anomalous ones." % len(seqs)
    print "\nPurging sequences that are not at least %d nucleotides long..." % minlength
    seqs = [(head, seq) for (head, seq) in seqs if len(seq) >= minlength]
    print "Retained %d sequences that are at least %d nucleotides long." % (len(seqs), minlength)
    if purgeambiguous:
        print "\nPurging any sequences with ambiguous nucleotide identities."
        m = re.compile('^[ATCGatcg]+$')
        cleanseqs = []
        for (head, seq) in seqs:
            if m.search(seq):
                cleanseqs.append((head, seq))
        seqs = cleanseqs
        print "Retained %d sequences after purging those with ambiguous nucleotides." % len(seqs)
    print "\nMaking sure the host / year / subtype can be parsed..."
    seqs = [(ParseHeader(head), seq) for (head, seq) in seqs if ParseHeader(head)]
    print "Retained %d sequences for which this information could be parsed." % len(seqs)
    print "\nOnly keep sequences of %s subtype..." % subtype
    subtype_match = re.compile(subtype)
    seqs = [(head_tup, seq) for (head_tup, seq) in seqs if subtype_match.search(head_tup[3])]
    print "Retained %d sequences of subtype %s." % (len(seqs), subtype)

    print "\nParsing sequences by host and year..."
    seqs_by_host_subtype_year = {}
    for ((strain, host, date, subtype), seq) in seqs:
        year = int(round(date))
        if only_keep_every_N_years:
            if not(year % every_N_yr) == 0:
                continue # only proceed every N years. ie if only keep every_n_yr = 2, 1918 will be added, 1919 will not, ...
        
        if not (host in hosts_subtypes_years and subtype in hosts_subtypes_years[host] and year in hosts_subtypes_years[host][subtype]):
            continue # don't keep as not a valid host / year combination
        head = '%.2f_STRAIN_%s_HOST_%s_SUBTYPE_%s' % (date - hosts_subtypes_years[host][subtype][year], strain, host, subtype)
        chars_to_replace = {' ':'', "'":'', '"':''}
        for (char, replacement) in chars_to_replace.iteritems():
            head = head.replace(char, replacement)
        if host not in seqs_by_host_subtype_year:
            seqs_by_host_subtype_year[host] = {}
        if subtype not in seqs_by_host_subtype_year[host]:
            seqs_by_host_subtype_year[host][subtype] = {}
        if year not in seqs_by_host_subtype_year[host][subtype]:
            seqs_by_host_subtype_year[host][subtype][year] = []
        seqs_by_host_subtype_year[host][subtype][year].append((head, seq))
    nretained = 0
    for host in seqs_by_host_subtype_year.iterkeys():
        print "For host %s, retained following subtypes:" % host
        for subtype in seqs_by_host_subtype_year[host].iterkeys():
            iseqs = []
            for yearseqs in seqs_by_host_subtype_year[host][subtype].itervalues():
                iseqs += yearseqs
            nretained += len(iseqs)
            print "  %s: %d total sequences from %d different years." % (subtype, len(iseqs), len(seqs_by_host_subtype_year[host][subtype]))
    print "Overall, %d of the %d sequences were retained as an appropriate host / subtype / year." % (nretained, len(seqs))
    print "\nNow retaining just %d sequences per host / subtype / year, retaining those with the highest counts of that exact sequence for each year..." % nperyearhostsubtype
    seqs = []
    for host in seqs_by_host_subtype_year.iterkeys():
        for subtype in seqs_by_host_subtype_year[host].iterkeys():
            for year in seqs_by_host_subtype_year[host][subtype].iterkeys():
                iseq_d = {}
                for (head, seq) in seqs_by_host_subtype_year[host][subtype][year]:
                    if seq in iseq_d:
                        iseq_d[seq].append((head, seq))
                    else:
                        iseq_d[seq] = [(head, seq)]
                counts_iseqs = [(len(iseq_d[seq]), iseq_d[seq]) for seq in iseq_d.keys()]
                counts_iseqs.sort()
                counts_iseqs.reverse()
                i = 0
                while i < nperyearhostsubtype and i < len(counts_iseqs):
                    seqs.append(("%s_n%d" % (counts_iseqs[i][1][0][0], i + 1), counts_iseqs[i][1][0][1]))
                    i += 1
    print "Retained %d total sequences." % len(seqs)
    
    print "Removing redundant sequencies..."
    seqs = GetUnique(seqs)
    print "Retained %d sequences after removing any identical sequences." % len(seqs)

    print "\nNow translating and aligning..." 
    (prot_alignments, cds_alignments) = NeedleCDSandProtAlignments(refseq, seqs, needlecmd)
    print "Successfully aligned all sequences."
    print "\nWriting the aligned coding DNA sequences to %s." % outfile
    mapmuts.sequtils.WriteFASTA(cds_alignments, outfile)
