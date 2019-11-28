#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage:
    aureme utils --run=ID --compounds_to_sbml --cpd=/path/to/root/txt_file [-v]
    aureme utils --run=ID --gbk_to_faa --gbk=/path/to/root/txt_file [-v]

options:
    --run=ID    Pathname to the aureme directory.
  
"""

import configparser
import csv
import docopt
import eventlet
import mpwt
import os
import os.path
import re
import subprocess
import time

from padmet.utils.connection import sbmlGenerator, gbk_to_faa

from aureme.parse import parse_config_file


def command_help():
    print(docopt.docopt(__doc__))


def utils_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    verbose = args['-v']

    if args['--compounds_to_sbml'] and args['--cpd']:
            txt_file = args['--cpd']
            compounds_to_sbml_utils(run_id, txt_file, verbose)
            return
    if args['--gbk_to_faa'] and args['--gbk']:
        gbk_file = args['--gbk']
        gbk_to_faa_utils(run_id, gbk_file, verbose)

        
def compounds_to_sbml_utils(species, input_file, verbose):
    """
    It reads the input_file at TXT format and it creates a SBML file.
    """
    run_id = "/shared/" + species
    root_file = run_id + "/" + input_file
    txt_file = root_file + ".txt"
    sbml_file = root_file + ".sbml"
    if os.path.exists(txt_file):
        if verbose:
            print("------>RUNNING STEP : creating file %s.sbml" %input_file)
        sbmlGenerator.compound_to_sbml(txt_file, sbml_file, verbose)
    else:
        print(txt_file + " is missing.")
        exit()

        
def gbk_to_faa_utils(species, input_file, verbose):
    """
    It reads the input_file at GBK format and it creates a FASTA file.
    """
    run_id = "/shared/" + species
    root_file = run_id + "/" + input_file
    gbk_file = root_file + ".gbk"
    faa_file = root_file + ".faa"
    if os.path.exists(gbk_file):
        if verbose:
            print("------>RUNNING STEP : gbk to faa")
        gbk_to_faa.gbk_to_faa(gbk_file=gbk_file, output=faa_file,
                              verbose=verbose)
    else:
        print(faa_file + " is missing.")
        exit()
