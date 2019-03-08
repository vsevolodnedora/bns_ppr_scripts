# contains set of methods to load/dump select/sort/get files

from lists import *

class MakePath:

    def __init__(self):
        pass

    @staticmethod
    def outflow(sim, outflow = '_0'):
        return Paths.ppr_sims + sim + '/' + 'outflow' + outflow + '/'

    @staticmethod
    def collated(sim):
        return Paths.ppr_sims + sim + '/' + 'collated' + '/'

    @staticmethod
    def collated2(sim):
        return Paths.ppr_sims + sim + '/' + 'collated2' + '/'