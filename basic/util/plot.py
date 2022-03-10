import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


def get_color_mapping_df(ordered_labels, cmap_name='gist_ncar'):
    '''Get a dataframe that maps a list of ordered labels
    to a list of colors.
    '''
    cmap = matplotlib.cm.get_cmap(cmap_name)
    n_labels = len(ordered_labels)

    colors = [matplotlib.colors.to_hex(cmap(i / n_labels)) for i in range(n_labels)]
    return pd.DataFrame({'labels': ordered_labels, 'colors':colors})


def make_frameless_scatter_plot(X, Y, colors, size=0.1):
    '''Make a scatter plot without any axis frames.
    This function is useful for making plots like UMAPs.
    '''
    fig, ax = plt.subplots()
    ax.scatter(X, Y, c=colors, 
           s=size, marker='.', edgecolor='none', rasterized=True)

    ax.set_aspect('equal')
    ax.set_axis_off()
    return ax

def make_graph_plot(X, Y, edges, colors, size=0.1, linewidth=0.005):
    '''Plot a graph.'''
    fig, ax = plt.subplots()
    
    # Plot the edges
    for e in edges:
        ex = [X[e[0]], X[e[1]]] 
        ey = [Y[e[0]], Y[e[1]]]
        ax.plot(ex, ey, color='black', alpha=0.3, linewidth=linewidth, rasterized=True, zorder=1)
    
    # Plot the nodes
    ax.scatter(X, Y, c=colors, 
           s=size, marker='.', edgecolor='none', rasterized=True, zorder=2)
    
    ax.set_aspect('equal')
    ax.set_axis_off()
