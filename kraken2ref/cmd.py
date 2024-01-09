import argparse
import os
import shutil
import logging


# get the version number from a file that is created by setuptools_scm 
# when the package is installed. 
try:
    from .version import version as __version__
    from .version import version_tuple
except ImportError:
    __version__ = "unknown version"
    version_tuple = (0, 0, "unknown version")

# NOTE about imports and using this script in develpoment and after production.
# The script uses absolute imports that will only work after being installed by pip.
# To run the script during development, run as follows from top level directory of the repo:
# python -m kraken_flu.cmd { OPTIONS }

def args_parser():
    """
    Command line argument parser
    """        
    parser = argparse.ArgumentParser(
        description = "kraken2ref: extract reads to reference sequences from kraken2 outputs") 

    parser.add_argument(
        '-v', '--version', 
        action='version', 
        version='kraken2ref ' + __version__ )

    # TODO: add the parameters
    
    return parser
    
def main():
    
    args = args_parser().parse_args()

    # TODO: run the command

if __name__ == "__main__":
    exit(main())