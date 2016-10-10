import math
import numpy as np
import matplotlib.pyplot as plt
tan = np.tan


def geom_create():
    #Wing parameters where r - root, t - tip, b - span, le - leading edge, te - trailing edge, S - planform area, AR - aspect ratio 
    r_chord = 1
    taper_ratio = 2
    t_chord = r_chord/taper_ratio
    b = 10
    sweep_angle = 20 
    #Wing Co-ordinates
    r_le = (0.0, 0.0)
    r_te = (r_chord, 0.0)
    t_le = (0.5*b*tan(math.radians(sweep_angle)), 0.5*b) 
    t_te = (0.5*b*tan(math.radians(sweep_angle)) + t_chord, 0.5*b)
    S = 0.5*b*t_le[0] + b*(t_te[0]-t_le[0]) + 0.5*b*(r_chord-t_te[0])
    AR = b**2/S
    return r_le, r_te, t_le, t_te, sweep_angle

def mesh_create(r_le, r_te, t_le, t_te, sweep_angle):
    le_ygrid = np.linspace(r_le[1], t_le[1], 20)
    le_xgrid = list(map(lambda y : y*tan(math.radians(sweep_angle)), le_ygrid))
    le_xy_grid = [(x,y) for x,y in zip(le_xgrid, le_ygrid)]
    te_ygrid = le_ygrid[:]
    te_xgrid = list(map(lambda y : r_te[0]-(y*((r_te[0]-t_te[0])/t_te[1])), te_ygrid))
    te_xy_grid = [(x,y) for x,y in zip(te_xgrid, te_ygrid)]
    #Stores start point((x,y) in le) and end point((x,y) in te) of chord variation along span
    chord_endpoints = [[pos_le, pos_te] for pos_le,pos_te in zip(le_xy_grid, te_xy_grid)]
    #Stores grid points from le to te along span direction
    planform_grid = map(lambda (pos_le, pos_te) : [[x,y] for x,y in zip(np.linspace(pos_le[0], pos_te[0], 5), np.linspace(pos_le[1], pos_te[1], 5))], chord_endpoints)   
    return np.asarray(planform_grid)

def plot_mesh(planform_grid):
    plt.figure(1)
    for chord_line in planform_grid:
        chord_line_x = chord_line[:,0]
        chord_line_y = chord_line[:,1]
        plt.plot(chord_line_x, chord_line_y, 'r-')
    for i in range(len(planform_grid[0])):
        span_line = planform_grid[:,i]
        span_line_x = span_line[:,0]
        span_line_y = span_line[:,1]
        plt.plot(span_line_x, span_line_y, 'r-')
        



