import numpy as np
import cmath

# Constants
dbl = np.float64
pi = np.pi
mi0 = 4.0e-7 * pi
ii = complex(0.0, 1.0)

# Globals
coordenadas = None
nos_elemento = None
nos_hetero = None
interv_x = interv_z = dx = dz = area_elem = None
a = np.zeros(3, dtype=dbl)
b = np.zeros(3, dtype=dbl)
c = np.zeros(3, dtype=dbl)
nos = nos_x = nos_z = h_elements = n_elements = None
rho = mi = espess = rho_ap = None
omega = hz = hx = prof_hetero = rho_elemento = mi_elemento = mi_hetero = rho_hetero = None
nc = None

#===================================================================================
def get_data():
    global interv_x, interv_z, dx, dz, area_elem, nos, nos_x, nos_z, n_elements, sB
    with open('data.txt', 'r') as file:
        interv_x = dbl(file.readline())
        interv_z = dbl(file.readline())
        nos_x = int(file.readline())
        nos_z = int(file.readline())

    dx = interv_x / (nos_x - 1)
    dz = interv_z / (nos_z - 1)
    nos = nos_x * nos_z
    n_elements = (nos_x - 1) * (nos_z - 1) * 2
    area_elem = (dx * dz) / 2
    sB = nos_z + 1

#===================================================================================
def get_model():
    global nc, omega, rho, mi, espess, hetero_z, hetero_x, prof_hetero, mi_hetero, rho_hetero
    with open('modelo.txt', 'r') as file:
        frequencia = dbl(file.readline())
        nc = int(file.readline())
        rho = np.zeros(nc, dtype=dbl)
        mi = np.zeros(nc, dtype=dbl)
        espess = np.zeros(nc - 1, dtype=dbl)
        
        for i in range(nc):
            rho[i] = dbl(file.readline())
        for i in range(nc):
            mi[i] = dbl(file.readline())
        for i in range(nc - 1):
            espess[i] = dbl(file.readline())

        hetero_z = dbl(file.readline())
        hetero_x = dbl(file.readline())
        prof_hetero = dbl(file.readline())
        mi_hetero = dbl(file.readline())
        rho_hetero = dbl(file.readline())

    omega = 2.0 * pi * frequencia

#===================================================================================
def grid_2d():
    global coordenadas, nos_elemento
    coordenadas = np.zeros((nos, 2), dtype=dbl)
    nos_elemento = np.zeros((n_elements, 3), dtype=int)
    
    j = 0
    for k in range(nos_x):
        for i in range(nos_z):
            coordenadas[j, 0] = k * dx
            coordenadas[j, 1] = i * dz
            j += 1

    k = 0
    temp = nos_z
    for j in range(0, n_elements, 2):
        nos_elemento[j, 0] = k
        nos_elemento[j, 1] = k + 1
        nos_elemento[j, 2] = k + nos_z
        k += 1
        if k == temp:
            k += 1
            temp += nos_z

    k = 1
    temp = nos_z + 1
    for j in range(1, n_elements, 2):
        nos_elemento[j, 0] = k
        nos_elemento[j, 1] = k + nos_z
        nos_elemento[j, 2] = k + nos_z - 1
        k += 1
        if k == temp:
            k += 1
            temp += nos_z

#===================================================================================
def hetero_ident():
    global nos_hetero, h_elements
    hetero_zz = int(hetero_z / dz)
    hetero_xx = int(hetero_x / dx)
    h_elements = (hetero_zz * hetero_xx) * 2

    nos_hetero = np.zeros((h_elements, 3), dtype=int)

    elz = 2 * (nos_z - 1)
    el_prof_hetero = 2 * (int(prof_hetero / dz))

    elx = int((nos_x - 1) / 2)
    elehx = int(int(hetero_x / dx) / 2)
    elx = elx - elehx
    elemento1 = elx * elz + el_prof_hetero + 1

    a = elemento1
    b = a + (2 * hetero_zz) - 1
    k = 0
    for j in range(a, b + 1):
        jj = j
        for i in range(k, h_elements, 2 * hetero_zz):
            nos_hetero[i, 0] = nos_elemento[jj, 0]
            nos_hetero[i, 1] = nos_elemento[jj, 1]
            nos_hetero[i, 2] = nos_elemento[jj, 2]
            jj += elz
        k += 1

#===================================================================================
def ondaplana(l):
    Hy = np.zeros(3, dtype=complex)
    Ex = np.zeros(3, dtype=complex)
    
    u_j = np.zeros(nc, dtype=complex)
    imped_intr = np.zeros(nc, dtype=complex)
    imped_ap = np.zeros(nc, dtype=complex)
    Hj = np.ones(nc, dtype=complex)

    for k in range(nc):
        permeabilidade = mi[k] * mi0
        u_j[k] = cmath.sqrt(ii * omega * permeabilidade / rho[k])
        imped_intr[k] = u_j[k] * rho[k]

    imped_ap[nc - 1] = imped_intr[nc - 1]

    for xx in range(3):
        Hy[xx] = Hj[nc - 1] * cmath.exp(-u_j[nc - 1] * coordenadas[nos_elemento[l, xx], 1])
        Ex[xx] = imped_intr[nc - 1] * Hy[xx]

    return Hy, Ex

#===================================================================================
def parameters(i, j):
    global rho_elemento, mi_elemento
    fonte_local1 = np.zeros(3, dtype=complex)
    fonte_local2 = np.zeros(3, dtype=complex)
    
    a[0] = coordenadas[nos_elemento[i, 1], 0] * coordenadas[nos_elemento[i, 2], 1] - coordenadas[nos_elemento[i, 2], 0] * coordenadas[nos_elemento[i, 1], 1]
    a[1] = coordenadas[nos_elemento[i, 2], 0] * coordenadas[nos_elemento[i, 0], 1] - coordenadas[nos_elemento[i, 0], 0] * coordenadas[nos_elemento[i, 2], 1]
    a[2] = coordenadas[nos_elemento[i, 0], 0] * coordenadas[nos_elemento[i, 1], 1] - coordenadas[nos_elemento[i, 1], 0] * coordenadas[nos_elemento[i, 0], 1]

    b[0] = coordenadas[nos_elemento[i, 1], 1] - coordenadas[nos_elemento[i, 2], 1]
    b[1] = coordenadas[nos_elemento[i, 2], 1] - coordenadas[nos_elemento[i, 0], 1]
    b[2] = coordenadas[nos_elemento[i, 0], 1] - coordenadas[nos_elemento[i, 1], 1]

    c[0] = coordenadas[nos_elemento[i, 2], 0] - coordenadas[nos_elemento[i, 1], 0]
    c[1] = coordenadas[nos_elemento[i, 0], 0] - coordenadas[nos_elemento[i, 2], 0]
    c[2] = coordenadas[nos_elemento[i, 1], 0] - coordenadas[nos_elemento[i, 0], 0]

    if np.array_equal(nos_elemento[i], nos_hetero[j]):
        Hy, Ex = ondaplana(i)
        delta_mi = (mi_hetero - mi[nc - 1]) * mi0
        delta_sigma = 1.0 / rho_hetero - 1.0 / rho[nc - 1]
        rho_elemento = rho_hetero
        mi_elemento = mi_hetero * mi0
        j += 1
    else:
        delta_mi = 0.0
        delta_sigma = 0.0
        rho_elemento = rho[nc - 0]

