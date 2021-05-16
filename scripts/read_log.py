""" Analysis of the simulation logs.

An examplary script for loading and the basic analysis of the simulation
log records. Functionality is similar to the accompanying jupyter
notebook.
"""


def import_log_files(path, runs):

    import sys
    from records import Records

    pat = str(runs[0]) + ' : ' + str(runs[1])

    fn = lambda k: path + 'log_m_' + str(k) + '.txt'
    append_static_recs = \
        lambda rs: [Records.runs_read_in.append(u) for u in rs.runs]

    cf, r = _read_log(fn(runs[0]))
    recs = [r]
    append_static_recs(r)
    for i in range(runs[0]+1, runs[1]+1):
        c, r = _read_log(fn(i))
        if c == cf:
            recs.append(r)
            append_static_recs(r)
        else:
            print(f'Error: Runs between {runs[0]} and {runs[1]} do use '
                  f'differing configurations at run {i}. \nExiting!')
            sys.exit(-1)
#    ravg = Records.average(recs)
    Records.scale_time_to(recs, 's')

    return recs, pat


def _read_log(fname):
    """ Read a log file to produce instance of Config and Records.

    :param fname: Name of the log file
    """

    from config import Config
    from records import Records

    print('Reading data from: ' + fname)

    cf = Config()
    rs = Records()

    with open(fname, 'r') as f:
        cf.readin(f)
        while True:
            h = f.readline().split()
            if h[0] == 'RUN':
                break
        rs.runs = [int(h[2])]
        rs.seed = [int(f.readline().split()[2])]
        for _ in range(2):
            f.readline()
        while True:
            r = f.readline()
            if r == '' or r[0] == '\n':
                break
            rs.add(r)

    print(f'read in: {len(rs.inds)} records')

    return cf, rs


def plot_timedata(name, x, y, fit=None, n=1,
                  figsize=None, labels=None):
    """ Plot data 'y' and eventually 'fit', both specified at 'x'.
    """

    import matplotlib.pyplot as plt
    from records import Records

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # plot the data:
    for j in range(n):
        ax.plot(x[j], y[j], ls='--', lw=0.5, label=labels[j])
    # plot the fit:
    if fit is not None:
        for j in range(n):
            ax.plot(x[j], fit[j], ls='-', lw=.5)

    max_x = max([xx[-1] for xx in x])
    ax.set_xlim([0., max_x*1.2])    # to accommodate the legend

    ax.legend()
    ax.set_xlabel(Records.time_label)
    plt.grid(True)
    plt.title(name)
    plt.show()


def main():

    from records import Records

    # Set the directory to the log files and the min, max
    # Monte Carlo run indexes:
    path = '../tests/'
    run_first = 28
    run_last = 29

    # Import data from the files:
    recs, pat = import_log_files(path, [run_first, run_last])

    # The log file is a record of time-dependent parameters.
    # Thsese evolve in real time measured in seconds.
    # The correct time values are ensured by application
    # of the Gillespie algorithm.
    # Acceptable units: 'd', 'hours', 'min', 's', 'secs'.
    Records.scale_time_to(recs, 's')

    # Prepare some vatiables for common use:
    # extract time for convenience:
    t = [r.t for r in recs]
    # set line lables to reflect run indexes:
    labels = ['run ' + str(i) for i in range(len(recs))]

    # System free energy is represented by the reaction scores:
    vv = [[sc['val'] for sc in r.score.values()] for r in recs]
    scores_total = [[sum(u) for u in zip(*v)] for v in vv]

    plot_timedata('total scores', t, scores_total, n=len(recs), labels=labels)

    # Plot evolution of the the number of nodes, by node degree (1 to 3):
    Records.plot_nodes(recs, pat)

    # ... and the number of segments, by segment type. The type is
    # specified by degrees of the nodes:
    # the reaction permit four segment tpes: 11, 13, 33 and 22 (the
    # latter designetes a disconnected cycle)
    Records.plot_segments_by_type(recs, pat)

    # Here is the total number of segments:
    plot_timedata('total number of segments',
                  t, [r.mtn for r in recs],
                  n=len(recs), labels=labels)

    # and the number of segment clusters (disconnected graph components):
    plot_timedata('number of clusters',
                  t, [r.cln for r in recs],
                  n=len(recs), labels=labels)


if __name__ == '__main__':

    main()
