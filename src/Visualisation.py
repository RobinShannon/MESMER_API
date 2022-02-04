import numpy as np
from scipy.interpolate import BSpline
from scipy.interpolate import splrep
import matplotlib.pyplot as plt

def plotPES(me):

    # Setup empy lists holding scatter plot and line plot objects
    scatterArrayX= []
    scatterArrayY= []
    lines = []

    # set up an empty matplotlib object
    fig = plt.figure()
    axs = fig.gca()
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.spines['bottom'].set_visible(False)
    ax = plt.axes()

    # Setup counter to space species along x-axis, every time a new species is added increse the counter and add the
    # x-axis value of the new species to encountered_dict
    x_axis_counter = 0
    encountered_dict = {}

    # Loop over reactions in me object
    for i in me.reactions_dict.keys():
        # For each reaction we fit a 3 point spline to the reactant, ts and product energies.
        # Need two vectors; x-axis = [x1,x2,x3] and y-axis = [y1,y2,y3]
        # Get x-axis value for reactants and then get y-axis value from the combined zpe of the reactants
        r = me.reactions_dict[i].reacs[0].name
        if r in encountered_dict:
            x1 = encountered_dict[r]
        else:
            x1 = x_axis_counter
            encountered_dict[r] = x1
            x_axis_counter += 10
        y1 = 0
        for reac in me.reactions_dict[i].reacs:
            y1 += reac.zpe

        # Do the same for the products
        p = me.reactions_dict[i].prods[0].name
        if p in encountered_dict:
            x3 = encountered_dict[p]
        else:
            x3 = x_axis_counter
            encountered_dict[p] = x3
            x_axis_counter += 10
        y3 = 0
        for prod in me.reactions_dict[i].prods:
            y3 += prod.zpe


        # x-axis for TS is mid point between x1 and x2
        if x1 > x3:
            x2 = x3 + ((x1 - x3) / 2)
        else:
            x2 = x1 + ((x3 - x1) / 2)

        # get ts energy if present, if not set it to between the reactant and product energy
        if me.reactions_dict[i].ts == 'none':
            if y1 > y3:
                y2 = y3 +  ((y1 - y3) / 4)
            else:
                y2 = y1 + ((y3 - y1) / 4)
        else:
            y2 = me.reactions_dict[i].ts.zpe
            encountered_dict[me.reactions_dict[i].ts.name] = x2

        # Now we have the x and y vector for the spline
        x = [x1,x2,x3]
        y = [y1,y2,y3]

        # Make sure x-axis is numerically increasing
        if x1 > x3:
            x.reverse()
            y.reverse()

        if me.reactions_dict[i].ts != 'none':
            x1_5 = x[0] + ((x[1] - x[0]) / 2)
            x2_5 = x[1] + ((x[2] - x[1]) / 2)
            y1_5 = y[0] + ((y[1] - y[0])*0.75)
            y2_5 = y[1] + ((y[2] - y[1])*0.75)
            x1 = [x[0], x1_5, x[1]]
            x2 = [x[1], x2_5, x[2]]
            y1 = [y[0], y1_5, y[1]]
            y2 = [y[1], y2_5, y[2]]

            t, c, k = splrep(x1, y1, s=0, k=2)
            x1_smooth = np.linspace(min(x1), max(x1), 150)
            y1_spline = BSpline(t, c, k, extrapolate=False)
            y1_smooth = y1_spline(x1_smooth)
            t, c, k = splrep(x2, y2, s=0, k=2)
            x2_smooth = np.linspace(min(x2), max(x2), 150)
            y2_spline = BSpline(t, c, k, extrapolate=False)
            y2_smooth = y2_spline(x2_smooth)
            x_smooth = np.concatenate([x1_smooth,x2_smooth])
            y_smooth = np.concatenate([y1_smooth,y2_smooth])
        else:
            t, c, k = splrep(x, y, s=0, k=2)
            x_smooth = np.linspace(min(x), max(x), 150)
            y_spline = BSpline(t, c, k, extrapolate=False)
            y_smooth = y_spline(x_smooth)

        # Make a line plot of x_smooth and y_smooth and add the line to the list of lines in the plot
        line = ax.plot(x_smooth, y_smooth, c='teal', linewidth=3.0, alpha=1, zorder=1)
        lines.append(line[0])

    # Now for the scatter plot
    # Loop over molecule in me_mols_dict and only create scatter point if that species is in encountered_dict
    for mol in me.mols_dict.values():
        if mol.name in encountered_dict:
            scatterArrayX.append(encountered_dict[mol.name])
            scatterArrayY.append(mol.zpe)

    s = ax.scatter(scatterArrayX, scatterArrayY, s=100, c='darkorange', alpha=1.0, zorder=1000)
    s.set_edgecolors('darkred')
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    plt.savefig('myFigure.png')