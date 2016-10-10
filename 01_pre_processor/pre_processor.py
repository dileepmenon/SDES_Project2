import math
import numpy as np
import matplotlib.pyplot as plt
tan = np.tan


def geom_create():
    # Wing parameters where r - root, t - tip, b - span, le - leading edge, te - trailing edge, S - planform area, AR - aspect ratio 
    r_chord = 1
    taper_ratio = 2
    t_chord = r_chord/taper_ratio
    b = 10
    sweep_angle = 20 
    # Wing edge co-ordinates
    r_le = (0.0, 0.0)
    r_te = (r_chord, 0.0)
    t_le = (0.5*b*tan(math.radians(sweep_angle)), 0.5*b) 
    t_te = (0.5*b*tan(math.radians(sweep_angle)) + t_chord, 0.5*b)
    S = 0.5*b*t_le[0] + b*(t_te[0]-t_le[0]) + 0.5*b*(r_chord-t_te[0])
    AR = b**2/S
    return r_le, r_te, t_le, t_te, sweep_angle


def mesh_create(r_le, r_te, t_le, t_te, sweep_angle):
    # Stores (x,y) locations of grid points along le and te
    le_ygrid = np.linspace(r_le[1], t_le[1], 20)
    le_xgrid = list(map(lambda y : y*tan(math.radians(sweep_angle)), le_ygrid))
    le_xy_grid = [(x,y) for x,y in zip(le_xgrid, le_ygrid)]
    te_ygrid = le_ygrid[:]
    te_xgrid = list(map(lambda y : r_te[0]-(y*((r_te[0]-t_te[0])/t_te[1])), te_ygrid))
    te_xy_grid = [(x,y) for x,y in zip(te_xgrid, te_ygrid)]
    # Stores start point and end point of chord variation along span
    chord_endpoints = [[pos_le, pos_te] for pos_le,pos_te in zip(le_xy_grid, te_xy_grid)]
    # Stores grid points from le to te along span direction
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


def midpoint(p1, p2):
    return ((p1[0]+p2[0])/2.0, (p1[1]+p2[1])/2.0) 


class panel:
    def __init__(self, (x1,y1), (x2,y2), (x3,y3), (x4,y4)):
        # Attributes of panel
        self.pos1 = (x1,y1)
        self.pos2 = (x2,y2)
        self.pos3 = (x3,y3)
        self.pos4 = (x4,y4)
        self.mid_pos1_pos4 = midpoint(self.pos1, self.pos4) 
        self.mid_pos2_pos3 = midpoint(self.pos2, self.pos3)
        self.mid_panel = midpoint(self.mid_pos1_pos4, self.mid_pos2_pos3) 
        # Control point of panel
        self.pos_cp = midpoint(self.mid_panel, self.mid_pos2_pos3)
        # Co-ordinates of bound and trailing vortices
        self.pos_c_by_4 = midpoint(self.mid_pos1_pos4, self.mid_panel)        
        self.pos_trail_vor1 = (self.pos_c_by_4[0], self.pos1[1])
        self.pos_trail_vor2 = (self.pos_c_by_4[0], self.pos4[1])
    
    def downwash_bound_vortex(self, panel_2):
        pass
    
    def downwash_trailing_vortex_1(self, panel_2):
        pass
    
    def downwash_trailing_vortex_2(self, panel_2):
        pass
    
    def panel_cp_distance(self, panel_2):
        pass
        

def create_panels(planform_grid):
    pass
    
