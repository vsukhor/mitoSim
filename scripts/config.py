

class Config:
    """ Encapsulates configuration parameters.
    """

    def __init__(self):

        self.workingDir_out = None
        self.runs = None
        self.time_total = None
        self.logfreq = None
        self.savefreq = None
        self.edge_length = None
        self.mtmassini = None
        self.segmassini = None
        self.use_fission = None
        self.rate_fission = None
        self.use_11_fusion = None
        self.fusion_rate_11 = None
        self.use_12_fusion = None
        self.fusion_rate_12 = None
        self.use_1L_fusion = None
        self.fusion_rate_1L = None

    def readin(self, f):

        """ Read the configuration data from the log file 'f'
        """

        skip_lines = lambda n: [f.readline() for _ in range(n)]

        skip_lines(2)            # Runstarted, empty
        self.workingDir_out = f.readline().split()[2]
        self.runs = [int(f.readline().split()[2]),
                     int(f.readline().split()[2])]
        skip_lines(2)            # empty, reading...
        self.time_total = float(f.readline().split()[2])
        self.logfreq = int(f.readline().split()[2])
        self.savefreq = int(f.readline().split()[2])
        self.edge_length = float(f.readline().split()[2])
        self.mtmassini = int(f.readline().split()[2])
        self.segmassini = int(f.readline().split()[2])
        self.use_fission = bool(f.readline().split()[2])
        self.rate_fission = float(f.readline().split()[2])
        self.use_11_fusion = bool(f.readline().split()[2])
        self.fusion_rate_11 = float(f.readline().split()[2])
        self.use_12_fusion = bool(f.readline().split()[2])
        self.fusion_rate_12 = float(f.readline().split()[2])
        self.use_1L_fusion = bool(f.readline().split()[2])
        self.fusion_rate_1L = float(f.readline().split()[2])

    def __eq__(self, c):

        """ Equality operator.
        """

        if self.workingDir_out != c.workingDir_out: return False
        #       if self.runs != c.runs:             return False
        if self.time_total != c.time_total:         return False
        if self.logfreq != c.logfreq:               return False
        if self.savefreq != c.savefreq:             return False
        if self.edge_length != c.edge_length:       return False
        if self.mtmassini != c.mtmassini:           return False
        if self.segmassini != c.segmassini:         return False
        if self.use_fission != c.use_fission:       return False
        if self.rate_fission != c.rate_fission:     return False
        if self.use_11_fusion != c.use_11_fusion:   return False
        if self.fusion_rate_11 != c.fusion_rate_11: return False
        if self.use_12_fusion != c.use_12_fusion:   return False
        if self.fusion_rate_12 != c.fusion_rate_12: return False
        if self.use_1L_fusion != c.use_1L_fusion:   return False
        if self.fusion_rate_1L != c.fusion_rate_1L: return False

        return True

    summary_fields = [
        'time_total',
        'logfreq',
        'savefreq',
        'edge_length',
        'mtmassini',
        'segmassini',
        'use_fission',
        'rate_fission',
        'use_11_fusion',
        'fusion_rate_11',
        'use_12_fusion',
        'fusion_rate_12',
        'use_1L_fusion',
        'fusion_rate_1L',
    ]

    def to_summary(self):

        return \
            self.time_total, \
            self.logfreq, \
            self.savefreq, \
            self.edge_length, \
            self.mtmassini, \
            self.segmassini, \
            self.use_fission, \
            self.rate_fission, \
            self.use_11_fusion, \
            self.fusion_rate_11,\
            self.use_12_fusion, \
            self.fusion_rate_12, \
            self.use_1L_fusion, \
            self.fusion_rate_1L
