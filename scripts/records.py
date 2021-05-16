"""Container class for time-dependent variables recorded in the log file.
"""

import numpy as np


class Records:
    """ Container for time-dependent variables recorded in the log file.
    """

    #: Run indexes.
    runs_read_in = []

    time_unit = None
    time_label = None

    def __init__(self):

        #: Run indexes; if len(self.runs) > 1,
        #: all the fields except 'runs',
        self.runs = None

        # 'seed', 'lattice_dims' and 'inds' are average values over the 'runs'.
        #: Rng seeds used to produce the runs.
        self.seed = None

        #: MC iteration indexes.
        self.inds = []
        #: System time.
        self.time = []

        #: Rescaled time.
        self.t = []

        #: Time interval from preceding reaction event.
        self.tau = []

        #: Reaction type.
        self.rt = {'code': [],
                   'str':  []}

        #: Number of nodes.
        self.n = [{'val': [], 'fit': None, 'pars': None},
                  {'val': [], 'fit': None, 'pars': None},
                  {'val': [], 'fit': None, 'pars': None}]

        #: Segment numbers, by type.
        self.m = {'11': [],
                  '22': [],
                  '33': [],
                  '13': []}

        #: Total graph mass (in edges).
        self.mtm = []

        #: Number of segments, total.
        self.mtn = []

        #: Number of clusters.
        self.cln = []

        #: Number of reaction events executed and the curent propensity,
        #: for each reaction type:
        self.score = {'fission':  {'num': [], 'val': []},
                      'fusion11': {'num': [], 'val': []},
                      'fusion12': {'num': [], 'val': []},
                      'fusion1L': {'num': [], 'val': []}}

    def add(self, rec):
        """ Extract information from a string record 'rec'.

        Add the result to the already available data.
        """

        q = rec.split()

        self.inds.append(int(q[0]))
        self.time.append(float(q[2]))
        self.tau.append(float(q[4]))
        self.rt['code'].append(int(q[6]))
        self.rt['str'].append(q[7])
        [self.n[i]['val'].append(int(q[9+i])) for i in range(3)]
        [self.m[k].append(int(q[13+2*i])) for i, k in enumerate(self.m.keys())]
        self.mtm.append(int(q[21]))
        self.mtn.append(int(q[23]))
        self.cln.append(int(q[25]))
        self.score['fission']['num'].append(int(q[28]))
        self.score['fission']['val'].append(float(q[29]))
        self.score['fusion11']['num'].append(int(q[31]))
        self.score['fusion11']['val'].append(float(q[32]))
        self.score['fusion12']['num'].append(int(q[34]))
        self.score['fusion12']['val'].append(float(q[35]))
        self.score['fusion1L']['num'].append(int(q[37]))
        self.score['fusion1L']['val'].append(float(q[38]))

    @staticmethod
    def scale_time_to(recs, unit):
        """ Scale time to the desired unit 'unit'.

        Acceptable units: 'd', 'hours', 'min', 's', 'secs'
        """

        for r in recs:
            if unit == 'd':
                r.t = [t / 3600 / 24 for t in r.time]
            elif unit == 'hours':
                r.t = [t / 3600 for t in r.time]
            elif unit == 'min':
                r.t = [t / 60 for t in r.time]
            elif unit == 's' or unit == 'sec':
                r.t = [t for t in r.time]
            else:
                assert (False, 'Wrong time unit')

        Records.time_unit = unit
        Records.time_label = 'Time (' + unit + ')'

    @staticmethod
    def plot_nodes(recs, pat, with_fit=False, figsize=None):
        """ Plot time evolution of the number of nodes.
        """

        import matplotlib.pyplot as plt

        ls = [':', '--', '-.']
        colors = plt.cm.rainbow(np.linspace(0, 1, len(recs)))

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

        for j, r in enumerate(recs):
            c = colors[j]
            for i in range(3):
                label = 'run ' + str(r.runs[0]) + ' n' + str(i+1)
                ax.plot(r.t, r.n[i]['val'], c=c, ls=ls[i], lw=0.5, label=label)
                if with_fit:
                    ax.plot(r.t, r.n[i]['fit'], c=c, lw=0.5)

        max_x = max([r.t[-1] for r in recs])
        ax.set_xlim([0., max_x * 1.3])  # to accommodate the legend
        ax.set_xlabel(Records.time_label)
        plt.legend()
        plt.grid(True)
        fig.suptitle('number of nodes: runs ' + pat)

        plt.show()

    @staticmethod
    def plot_segments_by_type(recs, pat, figsize=None):

        """ Plot time evolution of the number of segments.
        """

        import matplotlib.pyplot as plt

        ls = ['-', ':', '--', '-.']
        colors = plt.cm.rainbow(np.linspace(0, 1, len(recs)))

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

        for j, r in enumerate(recs):
            c = colors[j]
            for i, (k, v) in enumerate(r.m.items()):
                label = 'run ' + str(r.runs[0]) + ' m' + k
                ax.plot(r.t, v, c=c, ls=ls[i], lw=0.5, label=label)

        max_x = max([r.t[-1] for r in recs])
        ax.set_xlim([0., max_x * 1.3])  # to accommodate the legend
        ax.set_xlabel(Records.time_label)
        plt.legend()
        plt.grid(True)
        fig.suptitle('number of segments: runs ' + pat)

        plt.show()
