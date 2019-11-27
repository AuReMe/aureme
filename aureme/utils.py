#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage:
    aureme utils --run=ID --compounds_to_sbml --cpd=/path/to/root/txt_file [-v]

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

from aureme.parse import parse_config_file


def command_help():
    print(docopt.docopt(__doc__))


def utils_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    verbose = args['-v']

    if args['--compounds_to_sbml'] and args['--cpd']:
            txt_file = args['--cpd']
            compounds_to_sbml(run_id, txt_file, verbose)

    
def compounds_to_sbml(species, txt_file, verbose):
    """
    It read the txt_file and it creates a sbml file.
    """
    run_id = "/shared/" + species
    config_data = parse_config_file(run_id)
    padmet_u = config_data['padmet_u']
    root_file = run_id + "/" + txt_file

    if os.path.exists(root_file + ".txt"):
        print('"------>RUNNING STEP : creating file %s.sbml' %txt_file)
        cmds = ["python3.7", padmet_u + "/connection/sbmlGenerator.py",
                "--compound=" + root_file + ".txt",
                "--output=" + root_file + ".sbml"]
        if verbose:
            cmds.append('-v')
            subprocess.call(cmds)
    else:
        print(root_file + ".txt is missing.")
        exit()
