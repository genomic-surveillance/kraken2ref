import re
import os.path
from cached_property import cached_property
import logging

logging.basicConfig( format='%(asctime)s %(message)s', level=logging.DEBUG )

class KrakenTaxonomyReport():   
    """
    This class handles a kraken2 taxonomy report output file.
    
    Parameters:
        in_file: str, required
            path to the kraken2 taxonomy report file (should be called "report.txt")
            
        min_abs_reads: int, optional, defaults to 5
            minimum absolute number of reads that are directly assigned to a taxon. A taxon with
            fewer reads will be discarded.
        
    """
    
    def __init__( self, in_file: str, min_abs_reads: int = 5 ):
        self.in_file = in_file
        
        if not os.path.isfile( self.in_file ):
            raise ValueError(f'path { in_file} does not exist or is not a file')
        
    @cached_property
    def taxonomy_graph( self ):
        """
        A graph representation of the taxonomy with number of reads that are assigned directly to each node
        """
        # TODO
        pass