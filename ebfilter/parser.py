#! /usr/bin/env python

from .run import *
import argparse

def create_parser():
    parser = argparse.ArgumentParser(prog = "EBFilter")
    
    parser.add_argument("--version", action = "version", version = "EBFilter-0.3.0")

    subparsers = parser.add_subparsers()

    ####################
    # EBFilter main
    main_parser = subparsers.add_parser("main",
                         help = "EBFilter main function",
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    main_parser.add_argument("targetMutationFile", metavar = "target.vcf", type = str,
                              help = "the path to the mutation file")
    
    main_parser.add_argument("targetBamPath", metavar = "target.bam", type = str,
                              help = "the path to the target bam file")
    
    main_parser.add_argument("controlBamPathList", metavar = "controlBam_list.txt", type = str,
                              help = "the list of paths to control bam files")
    
    main_parser.add_argument("outputPath", metavar = "output.vcf", type = str,
                              help = "the path to the output")
    
    main_parser.add_argument('-f', choices=['vcf', 'anno'], default = 'vcf',
                        help = "the format of mutation file vcf or annovar (tsv) format")
    
    main_parser.add_argument('-t', metavar = "thread_num", default='1', type=int,
                        help = "the number of threads")
    
    main_parser.add_argument('-q', metavar = "mapping_qual_thres", default='20', type=int,
                        help = "threshold for mapping quality for calculating base counts")
    
    main_parser.add_argument('-Q', metavar = "base_qual_thres", default='15', type=int,
                        help = "threshold for base quality for calculating base counts")
    
    main_parser.add_argument('--ff', metavar = "filter_flags", default='UNMAP,SECONDARY,QCFAIL,DUP', type=str,
                        help = "skip reads with mask bits set")

    main_parser.add_argument("--loption", help = "use samtools mpileup -l option", action='store_true', default=False)

    main_parser.add_argument("--region", default = '', type = str,
                        help = "restrict the chromosomal region for mutation. active only if loption is on")

    main_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    main_parser.add_argument('--panel', metavar = "panel", required=True, type=str,
                        help = "control panel database(Tabix VCF)")

    main_parser.set_defaults(func = ebfilter_main)


    ####################
    # Create Control Panel Database
    panel_parser = subparsers.add_parser("panel",
                         help = "Create control panel database",
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    panel_parser.add_argument("targetMutationFile", metavar = "target.vcf", type = str,
                              help = "the path to the mutation file")

    panel_parser.add_argument("controlBamPathList", metavar = "controlBam_list.txt", type = str,
                              help = "the list of paths to control bam files")

    panel_parser.add_argument("outputPathPrefix", metavar = "/path/to/output", type = str,
                             help = "the path to the output prefix")

    panel_parser.add_argument('-q', metavar = "mapping_qual_thres", default='20', type=int,
                        help = "threshold for mapping quality for calculating base counts")
    
    panel_parser.add_argument('-Q', metavar = "base_qual_thres", default='15', type=int,
                        help = "threshold for base quality for calculating base counts")
    
    panel_parser.add_argument('--ff', metavar = "filter_flags", default='UNMAP,SECONDARY,QCFAIL,DUP', type=str,
                        help = "skip reads with mask bits set")

    panel_parser.add_argument("--loption", help = "use samtools mpileup -l option", action='store_true', default=False)

    panel_parser.add_argument("--region", default = '', type = str,
                        help = "restrict the chromosomal region for mutation. active only if loption is on")

    panel_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    panel_parser.set_defaults(func = create_control_panel_database)
    
    
    return parser

