#!/usr/bin/env python3

# Author: William F Godoy godoywf@ornl.gov
# Documents the list of mfem examples for adios2stream output 
# that have Paraview support

# Usage: Run from mfem-build/examples
# ./run_adios2stream.py example_prefix
# Example:
# ./run_adios2stream.py ex2p

import argparse
import subprocess
import os
import glob
import shutil
from genericpath import exists

cases = {'ex1p': ['square-disc.mesh',
                  'star.mesh',
                  'star-mixed.mesh',
                  'escher.mesh',
                  'fichera.mesh',
                  'fichera-mixed.mesh',
                  'toroid-wedge.mesh',
                  'square-disc-p2.vtk -o 2',
                  'square-disc-p3.mesh -o 3',
                  'square-disc-nurbs.mesh -o -1',
                  'star-mixed-p2.mesh -o 2',
                  'disc-nurbs.mesh -o -1',
                  'pipe-nurbs.mesh -o -1',
                  'ball-nurbs.mesh -o 2',
                  'fichera-mixed-p2.mesh -o 2',
                  'star-surf.mesh',
                  'square-disc-surf.mesh',
                  'inline-segment.mesh',
                  'amr-quad.mesh',
                  'amr-hex.mesh',
                  'mobius-strip.mesh',
                  'mobius-strip.mesh -o -1 -sc'],

         'ex2p': ['beam-tri.mesh',
                  'beam-quad.mesh',
                  'beam-tet.mesh',
                  'beam-hex.mesh',
                  'beam-wedge.mesh',
                  'beam-tri.mesh -o 2 -sys',
                  'beam-quad.mesh -o 3 -elast',
                  'beam-quad.mesh -o 3 -sc',
                  'beam-quad-nurbs.mesh',
                  'beam-hex-nurbs.mesh'
                  ]
         }


def ArgParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("example", type=str, help='example executable')
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = ArgParser()
    example = str(args.example)
    run_command = ''

    if example in cases.keys():
        if example.endswith('p'):
            run_command = 'mpirun '
            run_prefix = '-n 4 ' + example + ' '
        else:
            run_command = './' + example
            run_prefix = ' '

        for case in cases[example]:
            run_args = run_prefix + '-m ../data/' + case
            subprocess.check_call(run_command + run_args, shell=True)
        
        outdir = 'mfem_out_'+example
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        for bp in glob.glob(r'*.bp'):
            shutil.move(bp, outdir+'/'+bp)
