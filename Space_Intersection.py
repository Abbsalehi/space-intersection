import sympy as sym
import numpy as np
from sympy import Matrix, diff

X, Y, Z = sym.symbols('X Y Z')
XL, YL, ZL = sym.symbols('XL YL ZL')

# IOP
x0 = 3.244
y0 = 2.344
f = 5.039
pixel_size = 1.61797e-3

# matching points
x1 = 3.55815646
y1 = 1.714258792
x2 = 3.676853814
y2 = 4.634558571
# initial X, Y, Z values
G_init = Matrix([[2463520], [7430870], [15]])
G = G_init
# EOPs for two images
EOPs = np.array([[0.34359, 0.241779, -2.190935, 2463539.499, 7430884.88, 73.023236],
                 [0.474031, 0.087711, -2.472626, 2463539.493, 7430916.1, 73.109874]])


def Form_M_matrix(Omega, Phi, Kappa):
# Form rotation matrix M
    M = Matrix([[sym.cos(Phi) * sym.cos(Kappa),
                 -sym.cos(Phi) * sym.sin(Kappa), sym.sin(Phi)],
                [sym.cos(Omega) * sym.sin(Kappa) + sym.sin(Omega) * sym.sin(
                    Phi) * sym.cos(Kappa),
                 sym.cos(Omega) * sym.cos(Kappa) - sym.sin(Omega) * sym.sin(
                     Phi) * sym.sin(Kappa),
                 -sym.sin(Omega) * sym.cos(Phi)],
                [sym.sin(Omega) * sym.sin(Kappa) - sym.cos(Omega) * sym.sin(
                    Phi) * sym.cos(Kappa),
                 sym.sin(Omega) * sym.cos(Kappa) + sym.cos(Omega) * sym.sin(
                     Phi) * sym.sin(Kappa),
                 sym.cos(Omega) * sym.cos(Phi)]])
    M = M.transpose()
    return M


def Form_Partial_Deriviative(M):
#Partial deriative with respect with X Y Z
    x_1 = x0 - (f * ((X - XL) * M[0, 0] + (Y - YL) * M[0, 1] + (Z - ZL) * M[0, 2]) / (
        ((X - XL) * M[2, 0] + (Y - YL) * M[2, 1] + (Z - ZL) * M[2, 2])))
    y_1 = y0 - (-f * ((X - XL) * M[1, 0] + (Y - YL) * M[1, 1] + (Z - ZL) * M[1, 2]) / (
        ((X - XL) * M[2, 0] + (Y - YL) * M[2, 1] + (Z - ZL) * M[2, 2])))

    dxy = Matrix([[diff(x_1, 'X'), diff(x_1, 'Y'), diff(x_1, 'Z')], [diff(y_1, 'X'), diff(y_1, 'Y'), diff(y_1, 'Z')]])
    return dxy


condition = True
i = 0

while condition:

    Omega = np.radians(EOPs[0, 0])
    Phi = np.radians(EOPs[0, 1])
    Kappa = np.radians(EOPs[0, 2])
    XL = EOPs[0, 3]
    YL = EOPs[0, 4]
    ZL = EOPs[0, 5]

    M = Form_M_matrix(Omega, Phi, Kappa)

    x_1 = x0 - (f * ((X - XL) * M[0, 0] + (Y - YL) * M[0, 1] + (Z - ZL) * M[0, 2]) / (
        ((X - XL) * M[2, 0] + (Y - YL) * M[2, 1] + (Z - ZL) * M[2, 2])))
    y_1 = y0 - (-f * ((X - XL) * M[1, 0] + (Y - YL) * M[1, 1] + (Z - ZL) * M[1, 2]) / (
        ((X - XL) * M[2, 0] + (Y - YL) * M[2, 1] + (Z - ZL) * M[2, 2])))

    dxy1 = Form_Partial_Deriviative(M)

    Omega = np.radians(EOPs[1, 0])
    Phi = np.radians(EOPs[1, 1])
    Kappa = np.radians(EOPs[1, 2])
    XL = EOPs[1, 3]
    YL = EOPs[1, 4]
    ZL = EOPs[1, 5]

    M = Form_M_matrix(Omega, Phi, Kappa)
    d_x_2 = x0 - (f * ((X - XL) * M[0, 0] + (Y - YL) * M[0, 1] + (Z - ZL) * M[0, 2]) / (
        ((X - XL) * M[2, 0] + (Y - YL) * M[2, 1] + (Z - ZL) * M[2, 2])))
    d_y_2 = y0 - (-f * ((X - XL) * M[1, 0] + (Y - YL) * M[1, 1] + (Z - ZL) * M[1, 2]) / (
        ((X - XL) * M[2, 0] + (Y - YL) * M[2, 1] + (Z - ZL) * M[2, 2])))
    dxy2 = Form_Partial_Deriviative(M)
    x_2_cap = d_x_2.subs({'X': G[0], 'Y': G[1], 'Z': G[2]})
    y_2_cap = d_y_2.subs({'X': G[0], 'Y': G[1], 'Z': G[2]})

    dxy1 = dxy1.subs({X: G[0], Y: G[1], Z: G[2]})
    dxy2 = dxy2.subs({X: G[0], Y: G[1], Z: G[2]})

    A = Matrix([dxy1, dxy2])

    delta_x_1 = x1 - x_1.subs({X: G[0], Y: G[1], Z: G[2]})
    delta_y_1 = y1 - y_1.subs({X: G[0], Y: G[1], Z: G[2]})
    delta_x_2 = x2 - d_x_2.subs({X: G[0], Y: G[1], Z: G[2]})
    delta_y_2 = y2 - d_y_2.subs({X: G[0], Y: G[1], Z: G[2]})

    delta_l = Matrix([[delta_x_1], [delta_y_1], [delta_x_2], [delta_y_2]])
    delta_X = ((np.transpose(A) * A).inv()) * (np.transpose(A) * delta_l)

    G = G + delta_X
    i += 1
    print(i)
    if np.max(np.abs(delta_X)) < 0.0001:
        # print optimaized X Y Z for the GCP
        print("Optimized XYZ =", G)
        condition = False
        break
