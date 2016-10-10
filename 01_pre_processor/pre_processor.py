import math
import numpy as np
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
    return r_le, r_te, t_le, t_te
