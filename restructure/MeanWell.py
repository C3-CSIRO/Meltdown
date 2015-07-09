# -*- coding: utf-8 -*-


class MeanWell:
    def __init__(self, tm, tmError, isComplex, replicates, contents):
        #relevant info for a mean well
        self.tm = tm
        self.tmError = tmError
        self.isComplex = isComplex
        self.replicates = replicates
        self.contents = contents
        return
