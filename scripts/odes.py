import numpy as np
import scipy.integrate as spi


# ======================================================================================================================
def node_numbers_equil(c1, c2, L):

    """ Equilibrium solution of the ODE system.
    """

    from numpy.polynomial import Polynomial as pol

    # Eq. for the free ends is a 3rd degree polynomial. Its roots are:
    n1 = pol([
            -2 * L,
            1 - c1,
            (1 - c2) * c1,
            c1 * c2
        ]).roots()

    # Because the solution represents the number of nodes, we are interested in positive real roots.
    # Of the three roots only one should correspond to these constraints.
    n1 = n1[n1.imag == 0].real
    n1 = n1[n1 >= 0]
    if len(n1) != 1:
        print(f'Warning: not unique solution for c1 = {c1}, c2 = {c2}, L = {L}: {n1}')
    else:
        n1 = n1[0]

    n2 = n1 * (n1 - 1) * c1 / 2     # number of bulk nodes
    n3 = 2 * n1 * n2 * c2 / 3       # number of three-way junctions

    return n1, n2, n3


# ======================================================================================================================
def x2(x13, L):

    """ Number of bulk nodes given the total graph size 'L'
        and the number of free ends x13[0] and the brnchings x13[1].
    """

    try:
        iter(x13[0])

    except TypeError:
        # x1, x3 are not iterable:
        return L - (x13[0] + 3 * x13[1]) / 2

    else:
        # x1, x3 are iterable:
        return np.array([L - (x1 + 3 * x3) / 2 for x1, x3 in zip(*x13)])


# ======================================================================================================================
def eqs(t, x, b, a1, a2, L):

    """ Right-hand side of the ODE system.
    """

    x1, x3 = x

    return np.array([
            (a2/2 - a1) * x1*x1 +
            (a1 - a2 * L - b) * x1 +
            1.5 * a2 * x1 * x3 -
            1.5 * b * x3 +
            2. * b * L,

            a2 * L * x1 -
            0.5 * a2 * x1*x1 -
            1.5 * a2 * x1 * x3 -
            1.5 * b * x3
    ])


# ======================================================================================================================
def _jacobian(x1, x3, a1, a2, b, L):

    """ Jacobian matrix for the ODE system.
    """

    return np.array([[a1 - b - a2 * L + (a2 - 2 * a1) * x1 + 1.5 * a2 * x3,   1.5 * (a2 * x1 - b)],
                     [a2 * (L - x1) - 1.5 * a2 * x3,                         -1.5 * (a2 * x1 + b)]])


# ======================================================================================================================
def plot_node_numbers(x, y, L, z3, figsize=None):

    """ Plot the node numbers 'z3' as a function of reduced parameters 'x' and 'y'.
    """

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=figsize)
    fig.suptitle(f'Number of nodes for L={L}')

    colors = ('r', 'g', 'b')
    xx, yy = np.meshgrid(np.log10(x), np.log10(y))
    zlim = max(max(max(z3)))

    for i in range(3):
        z = np.array([[oo[i] for oo in o] for o in z3])
        ax = fig.add_subplot(1, 3, i+1, projection='3d')
        ax.plot_surface(xx, yy, z, color=colors[i], linewidth=0.2, edgecolor='k', antialiased=False)
        ax.set_zlim(top=zlim)
        ax.set_xlabel('log(c1)')
        ax.set_ylabel('log(c2)')
        ax.set_zlabel('number of nodes')
        ax.set_title(f'degree {i+1}')

    fig.tight_layout()
    plt.show()


# ======================================================================================================================
def plot_p(x, y, L, z3):

    import plotly.graph_objects as go

    def trace(x, y, z3):

        colors = ['r', 'g', 'b']
        xx, yy = np.meshgrid(np.log10(x), np.log10(y))
        z = [np.array([[oo[i] for oo in o] for o in z3]) for i in range(2)]

        return  go.Surface(x=xx, y=yy, z=z[0], showscale=False, surfacecolor=colors[0], opacity=0.9, visible=False), \
                go.Surface(x=xx, y=yy, z=z[1], showscale=False, surfacecolor=colors[1], opacity=0.9, visible=False), \
                go.Surface(x=xx, y=yy, z=z[2], showscale=False, surfacecolor=colors[2], opacity=0.9, visible=False) \

    # Create figure
    fig = go.Figure()

    # Add traces, one for each slider step
    for i in range(len(L)):
        fig.add_trace(trace(x, y, z3[i]))

    fig.show()

# ======================================================================================================================
def plot_p_simple(x, y, L, z3):

    import plotly.graph_objects as go

    def trace(x, y, z3):

        colors = {'r' : [[0, 'rgb(255,0,0)'], [1, 'rgb(255,0,0)']],
                  'g' : [[0, 'rgb(0,255,0)'], [1, 'rgb(0,255,0)']],
                  'b' : [[0, 'rgb(0,0,255)'], [1, 'rgb(0,0,255)']]}
        xx, yy = np.meshgrid(np.log10(x), np.log10(y))
        z = [np.array([[oo[i] for oo in o] for o in z3]) for i in range(3)]

        return [go.Surface(x=xx, y=yy, z=z[0], showscale=False, colorscale=colors['r'], opacity=0.9, visible=True),
                go.Surface(x=xx, y=yy, z=z[1], showscale=False, colorscale=colors['g'], opacity=0.9, visible=True),
                go.Surface(x=xx, y=yy, z=z[2], showscale=False, colorscale=colors['b'], opacity=0.9, visible=True)]

    # Create figure
    fig = go.Figure()

    for ss in trace(x, y, z3[0]):
        fig.add_trace(ss)
    fig.update_layout(title_text="Node numbers")
    fig.show()

# ======================================================================================================================
def plot_phase_equil(x123, c2, L, figsize):

    """ Plot the phase space of node numbers.
    """

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib import colors as crs

    fig = plt.figure(figsize=figsize)
    fig.suptitle(f'Node numbers for L={L} at equilibrium.')
    ax = fig.add_subplot(1, 1, 1, projection='3d')

    line_colors = [cm.jet(c) for c in np.linspace(0.0, 1.0, len(x123[0]))]
    x = [[np.array([[oo[i] for oo in o] for o in xx]) for i in range(3)] for xx in x123]

    for xx in x:
        for c, x0, x1, x2 in zip(line_colors, *xx):
            ax.plot(x0, x1, x2, lw=0.2, color=c, marker='o', ms=3, fillstyle=None)


    ax.set_xlabel('nodes deg. 1')
    ax.set_ylabel('nodes deg. 2')
    ax.set_zlabel('nodes deg. 3')

    plt.colorbar(cm.ScalarMappable(norm=crs.Normalize(c2[0], c2[-1]), cmap=cm.jet),
                 ax=ax, label='c2: end-to-side fusion reduced rate')
    fig.tight_layout()
    plt.show()

# ======================================================================================================================
def is_stable(x, b, a1, a2, L):

    """ True iff the ODE system solution 'x' is asymptotically stable.
        'b', 'a1', 'a2' and 'L' are the parameter values.
    """

    import scipy.linalg as la

    n1, n2 = a1.shape[0], a2.shape[0]

    x = [np.array([[oo[i] for oo in o] for o in x]) for i in range(3)]
    jac = np.array([[_jacobian(x[0][i,j], x[2][i,j], a1[i], a2[j], b, L)
                    for j in range(n2)]
                   for i in range(n1)])
    lamda = np.array([[la.eigvals(oo) for oo in o] for o in jac])

    return [[np.all(lamda[i,j,:].real < 0.)
             for j in range(n2)]
            for i in range(n1)]


# ======================================================================================================================
def plot_stability(b, c1, c2, L, st):

    """ Visualize the stability 'st':
        blue marker -> is stable
        red marker -> is unstable
        'b', 'c1', 'c2' and 'L' are the parameter values.
    """

    import matplotlib.pyplot as plt

    fig = plt.figure()
    fig.suptitle(f'Solution stability for L={L}, b={b:.3g}')


    ax = fig.add_subplot(1, 1, 1)
    x, y = np.meshgrid(np.log10(c1), np.log10(c2))
    c = np.where(st, 'b', 'r')
    ax.scatter(x.flatten(), y.flatten(), s=10, c=c.flatten())
    ax.set_xlabel('log(c1)')
    ax.set_ylabel('log(c2)')

    plt.show()


# ======================================================================================================================
def plot_time_evol(b, a1, a2, L, sol, t, figsize=None):

    """ Plot time-dependent solution 'sol' versus time 't' for each of the node types.
        'b', 'a1', 'a2' and 'L' are the parameter values.
    """

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=figsize)
    fig.suptitle(f'Time evolution of node numbers for b={b:.3g}, a1={a1:.3g}, a2={a2:.3g}, L={L}')
    colors = ('r', 'g', 'b')

    for i in range(3):
        ax = fig.add_subplot(1, 3, i+1)
        for s in sol:
            ax.plot(t, s[i], c=colors[i], lw=0.5)
        ax.set_xlabel('time (sec)')
        ax.set_ylabel('number of nodes')
        ax.set_title(f'degree {i+1}')
        ax.set_facecolor('0.9')
        ax.grid(True)

    fig.tight_layout()
    plt.show()


#######################################################################################################################
if __name__ == '__main__':

    # Analysis of the deterministc approximation to the mitochondrial graph.

    # Initialize the parameters:

    def a(m):
        """Initialize fusion rate constants."""
        init = -np.floor(2 * m / 3)
        return np.logspace(init, init + m - 1, num=m, base=2)

    m1, m2, mL = 57, 57, 1      # grid dimensions
    b = 1                       # fission rate constant
    a1, a2 = a(m1), a(m2)       # fusion rate constant
    c1, c2 = a1/b, a2/b         # reduced rates

    L = [10000, 30000]          # total number of edges in the graph

    # Find the equilibrium solution as t -> Inf and plot it.
    x = [[[node_numbers_equil(cc1, cc2, l) for cc1 in c1] for cc2 in c2] for l in L]
#    [plot_node_numbers(c1, c2, l, xx, figsize=(15, 5)) for xx, l in zip(x, L)]

    # Plot the solution in phase coordinates:
    plot_phase_equil(x, c2, L, figsize=(10,10))

    # Examine the steady state solution.
    # The equilibrium is asymptotically stable if real parts of all eigenvalues of the Jacobian are strictly negative.
    st = [is_stable(xx, b, a1, a2, l) for xx, l in zip(x, L)]
    # Plot the stability map indicating stable and unstable solutions with blue and red markers respectively:
#    [plot_stability(b, c1, c2, l, s) for s, l in zip(st, L)]

    # Examine transient bechavior.
    # Slove the ODEs directly for specific parameters and plot the results:

    Lt = L[0]                                   # graph total mass (in edges)
    bt = b                                      # fission rate constant
    a1t = a1[20]                                # end-to-end fusion rate constant
    a2t = a2[30]                                # end-to-side fusion rate constant
    tspan = [0., 20.]                          # time interval
    tsol = np.linspace(tspan[0], tspan[1], 100) # time points for plotting
    # initial values:
    x1_0 = np.linspace(0, Lt, 10)
    x3_0 = (Lt - x1_0) / 2

    x123 = []
    for x10, x30 in zip(x1_0, x3_0):
        # new scipy ivp solver: requires scipy >= 1.4:
        sol = spi.solve_ivp(eqs, t_span=tspan, y0=[x10, x30], args=(b, a1t, a2t, Lt), dense_output=True)
        x13 = sol.sol(tsol)
        x123.append([x13[0,:], x2(x13, Lt), x13[1,:]])

    plot_time_evol(b, a1t, a2t, Lt, x123, tsol, figsize=(16, 5))




