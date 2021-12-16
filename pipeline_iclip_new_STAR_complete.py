#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=====================================
iCLIP processing and mapping pipeline
=====================================

"""
from ruffus import *
import sys
import os
from cgatcore import pipeline as P


PARAMS = P.get_parameters("pipeline_iclip.yml")

    
# Cutadapt
@follows (mkdir("processed"))
@transform ('*.fastq.gz', regex(r'(.*).fastq.gz'), r'processed/\1.trim.fastq.gz') 
def cutadapt (infile, outfile):
    ''' trims 3' adapter and removes low quality reads '''
    statement = ''' cutadapt -q %(cutadapt_minphred)s 
    --minimum-length %(cutadapt_minlength)s 
    -a %(general_adapter)s
    -o %(outfile)s %(infile)s
    '''
    P.run(statement)
  
    
# Remove UMI
@transform (cutadapt, suffix('.trim.fastq.gz'), '.processed.fastq.gz')
def umiprocess (infile, outfile):
    ''' removes UMI from start of read and moves to end of read name '''
    statement = ''' umi_tools extract -I %(infile)s 
    --bc-pattern=%(umiprocess_pattern)s
    --stdout=%(outfile)s
    '''
    P.run(statement)

    
# Demultiplexing ### need to figure out how to redirect to new folder
@split (umiprocess, 'demultiplexed/*.fastq')
def demux (infile, outfiles):    
    ''' demultiplex based on barcode at start of read using fastx-toolkit '''
    startfile = " ".join(infile)
    statement = '''  zcat %(startfile)s |
    fastx_barcode_splitter.pl --bcfile %(demux_barcodes)s 
    --prefix demux --suffix .fastq --bol 
    --mismatches %(demux_mismatches)i 
    --partial %(demux_partial)i 
    '''
    P.run(statement)


# FASTQC
@follows (mkdir("fastqc1"), demux)
@transform ('demux*.fastq', regex(r'demux(.*).fastq'), 'fastqc1/.fastqc')
def fastqc1 (infile,outfile):
    ''' does fastqc on reads after processing '''
    statement = ''' fastqc %(infile)s -o %(general_outputdir)s/fastqc1 > %(outfile)s
    '''
    P.run(statement)
    
    



# STAR mapping 
@follows (mkdir("STARmapped"), demux)
@transform ('demux*.fastq', regex(r'demux(.*).fastq'), r'STARmapped/\1.bam')
def STARmap (infile,outfile):
    ''' maps non-repetitive elements to genome '''
    outprefix = P.snip(outfile, ".bam")
    job_threads = PARAMS["STARmap_threads"]
    statement = '''
    STAR  --runMode alignReads
    --runThreadN %(job_threads)i
    --genomeDir %(STARmap_genome)s
    --readFilesIn %(infile)s
    --outFilterMultimapNmax 10
    --limitBAMsortRAM 10000000000
    --outFileNamePrefix %(outprefix)s
    --outSAMattributes All
    --outStd BAM_SortedByCoordinate
    --outSAMtype BAM SortedByCoordinate
    --outFilterType BySJout
    --outReadsUnmapped Fastx
    --outSAMattrRGline ID:foo
    --alignEndsType Extend5pOfRead1
    --outFilterScoreMinOverLread 0 
    --outFilterMatchNminOverLread 0 
    --outFilterMatchNmin 0
    --outFilterMismatchNoverReadLmax 0.06
    --outSAMmultNmax 1
    > %(outfile)s 
    '''
    P.run(statement)

    
    
# samtools index1
@transform (STARmap, suffix('.bam'), '.bam.bai')
def index1 (infile, outfile):
    statement = ''' samtools index %(infile)s
    '''
    P.run(statement)


# Deduplicate
@follows (index1)
@transform (STARmap, regex(r'(.*).bam'), r'\1.dedup.bam')
def dedup (infile,outfile):
    ''' deduplicate samples based on UMI using umi_tools '''
    statement = '''
    umi_tools dedup -I %(infile)s --output-stats=deduplicated -S %(outfile)s
    '''
    P.run(statement)


# samtools index2
@transform (dedup, suffix('.dedup.bam'), 'dedup.bam.bai')
def index2 (infile, outfile):
    ''' creates index of deduplicated bam file, generates .bai '''
    statement = '''samtools index %(infile)s
    '''
    P.run(statement)



@follows (cutadapt, umiprocess, demux, STARmap, index1, dedup, index2)
def full():
    pass
    
# ---------------------------------------------------
# Generic pipeline tasks


#
#def main(argv=None):
#    if argv is None:
#        argv = sys.argv
#    P.main(argv)
#
#
#if __name__ == "__main__":
#    sys.exit(P.main(sys.argv))
if __name__ == "__main__":
    sys.exit( P.main(sys.argv) )
