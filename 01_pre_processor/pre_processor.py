import math
import numpy as np
import matplotlib.pyplot as plt
tan = np.tan
pi = np.pi
sqrt = np.sqrt


def geom_create():
    # Wing parameters where r - root, t - tip, b - span, le - leading edge, te - trailing edge, S - planform area, AR - aspect ratio 
    r_chord = 1
    taper_ratio = 1
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
    planform_grid = list(map(lambda (pos_le, pos_te) : [[x,y] for x,y in zip(np.linspace(pos_le[0], pos_te[0], 5), np.linspace(pos_le[1], pos_te[1], 5))], chord_endpoints))   
    return np.asarray(planform_grid)


def plot_mesh(planform_grid):
    # Plots the mesh created on the surface of the given wing geometry
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
        plt.title('Discretisation of planform into panels')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid()

def midpoint(p1, p2):
    return ((p1[0]+p2[0])/2.0, (p1[1]+p2[1])/2.0) 


class Panel:
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
    
    def downwash_bound_vortex(self, panel_j):
        # Calculates the velocity induced on panel control point by bound vortex of panel j
        trail_vor1 = panel_j.pos_trail_vor1
        trail_vor2 = panel_j.pos_trail_vor2
        a = (1.0/(4.0*pi))*(1.0/coeff1(trail_vor1, trail_vor2))
        b = self.coeff2(trail_vor1, trail_vor2, 0)/self.panel_cp_distance(trail_vor1)
        c = self.coeff2(trail_vor1, trail_vor2, 1)/self.panel_cp_distance(trail_vor2)
        vel = a*(b-c)
        return vel

    def downwash_trailing_vortex1(self, panel_j):
        # Calculates the velocity induced on panel control point by first trailing vortex of panel j
        trail_vor1 = panel_j.pos_trail_vor1
        trail_vor2 = panel_j.pos_trail_vor2
        a = (1.0/(4.0*pi))*(1/(trail_vor1[1]-self.pos_cp[1]))
        b = (self.pos_cp[0]-trail_vor1[0])/self.panel_cp_distance(trail_vor1)
        vel = a*(1.0+b)
        return vel
    
    def downwash_trailing_vortex2(self, panel_j):
        # Calculates the velocity induced on panel control point by second trailing vortex of panel j
        trail_vor1 = panel_j.pos_trail_vor1
        trail_vor2 = panel_j.pos_trail_vor2
        a = (-1.0/(4.0*pi))*(1.0/(trail_vor2[1]-self.pos_cp[1]))
        b = (self.pos_cp[0]-trail_vor2[0])/self.panel_cp_distance(trail_vor2)
        vel = a*(1.0+b)
        return vel
    
    def panel_cp_distance(self, pos_trail_vor):
        # Calculates distance between control point of panel with the trailing vortex co-ordinates of panel j 
        return sqrt((self.pos_cp[0]-pos_trail_vor[0])**2 + (self.pos_cp[1]-pos_trail_vor[1])**2) 
    
    def coeff1(self, pos_vor1, pos_vor2):
        return (self.pos_cp[0]-pos_vor1[0])*(self.pos_cp[1]-pos_vor2[1]) - (self.pos_cp[0]-pos_vor2[0])*(self.pos_cp[1]-pos_vor1[1])
    
    def coeff2(self, pos_vor1, pos_vor2, num):
        pos_vor_num = [pos_vor1, pos_vor2]
        return (pos_vor2[0]-pos_vor1[0])*(self.pos_cp[0]-pos_vor_num[num][0]) - (pos_vor2[1]-pos_vor1[1])*(self.pos_cp[1]-pos_vor_num[num][1])
        
        
def create_panels(planform_grid):
    # Initializes the panels class and stores the list of panel classes 
    list_of_starboard_panels = []
    list_of_port_panels = []
    for i,j in zip(planform_grid[:-1], planform_grid[1:]):
        for k in range(len(i)-1):
            list_of_starboard_panels.append(Panel((i[k][0],i[k][1]),(i[k+1][0],i[k+1][1]),(j[k+1][0],j[k+1][1]),(j[k][0],j[k][1])))
            list_of_port_panels.append(Panel((i[k][0],-1*i[k][1]),(i[k+1][0],-1*i[k+1][1]),(j[k+1][0],-1*j[k+1][1]),(j[k][0],-1*j[k][1])))
    list_of_all_panels = list_of_starboard_panels+list_of_port_panels
    return list_of_all_panels
    
