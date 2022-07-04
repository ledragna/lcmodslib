#!/usr/bin/env python3
"""
Example of Json
{"qmdata": {"functional": "B3PW91",
                    "basis": "gen 6D 10F empiricaldispersion=gd3bj",
                    "molchr": 0,
                    "molspn": 1,
                    "addroot": "freq=vcd ",
                    "addline": "@/home/m.fuse/basis/SNSD.gbs"}
                    }
"""

import os
import sys
import argparse
import copy
import json
import warnings

from lcmodslib.base import gio
from lcmodslib.base import gmanip

## Generic functions
def read_optfile(fname):
    """
    Read the option file. expected in json style
    """
    if not os.path.exists(fname):
        raise FileNotFoundError('Missing JSON file')
    with open(fname, 'r', encoding='UTF-8') as fopen:
        data = json.load(fopen)
    return data

def build_parser():
    """Builds options parser.
    Builds the full option parser.
    Returns
    -------
    :obj:`ArgumentParser`
        `ArgumentParser` object
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('fname', type=str, help="The gaussian log file with the equilibrium geometry")
    parser.add_argument('optfile', nargs='?', help='option file in json with QM options')
    parser.add_argument('patom', default=1, type=int,
                        help='Pivot atom on the solute around which look for solvent molecules')
    parser.add_argument('solattype', defoult='OW', type=str,
                        help='atom type of the solvent to search (es. "OW")')
    parser.add_argument('-ns', '--nsolvent', default=0, type=int,
                        help='number of solvent molecule to be included in QM')
    txt = """Reference frame type:
        single: only one frame used usually the cluster centroid
        all: all the frame in the trajectory
        """
    parser.add_argument('-m', '--mode', choices=['single', 'all'],
                        default='single', help=txt)
    parser.add_argument('--frame', default=-1, type=int, help='the frame to be\
                        used as reference in single mode')
    parser.add_argument('-p', '--prefix', default='fromtraj',
                        type=str, help='prefix name for the output files')
    parser.add_argument('--silent', action='store_true', help="Suppres almost all the printing")
    return parser