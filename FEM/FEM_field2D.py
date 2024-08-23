import numpy as np
import logging
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
    with open('/data.txt', 'r') as file:
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
        rho_elemento = rho[nc - 1]
        mi_elemento = mi[nc - 1] * mi0

    fonte_local1[0] = -ii * omega * delta_mi * area_elem / 12.0 * (2 * Hy[0] + Hy[1] + Hy[2]) + rho_elemento * delta_sigma / 6.0 * (Ex[0] + Ex[1] + Ex[2]) * c[0]
    fonte_local1[1] = -ii * omega * delta_mi * area_elem / 12.0 * (Hy[0] + 2 * Hy[1] + Hy[2]) + rho_elemento * delta_sigma / 6.0 * (Ex[0] + Ex[1] + Ex[2]) * c[1]
    fonte_local1[2] = -ii * omega * delta_mi * area_elem / 12.0 * (Hy[0] + Hy[1] + 2 * Hy[2]) + rho_elemento * delta_sigma / 6.0 * (Ex[0] + Ex[1] + Ex[2]) * c[2]

    fonte_local2[0] = rho_elemento * delta_sigma / 6.0 * c[0]
    fonte_local2[1] = rho_elemento * delta_sigma / 6.0 * c[1]
    fonte_local2[2] = rho_elemento * delta_sigma / 6.0 * c[2]

    return fonte_local1, fonte_local2

#===================================================================================
def cc_origin(Fonte, Global):
    borda = nos - (nos_z - 1)
    Fonte[:nos_z] = 0.0 + 0.0j
    Fonte[borda:] = 0.0 + 0.0j

    for i in range(nos_z, nos - (nos_z - 1), nos_z):
        Fonte[i] = 0.0 + 0.0j
        Global[sB - 1, i] = 10.0 + 10.0j

    for j in range(2 * nos_z, nos, nos_z):
        Fonte[j] = 0.0 + 0.0j
        Global[sB - 1, j] = 10.0 + 10.0j

    Global[sB - 1, :nos_z] = 10.0 + 10.0j
    Global[sB - 1, borda:] = 10.0 + 10.0j

#===================================================================================
def cc(Fonte, Global):
    borda = nos - (nos_z - 1)
    Fonte[:nos_z] = 0.0 + 0.0j
    Fonte[borda:] = 0.0 + 0.0j

    for i in range(nos_z, nos - (nos_z - 1), nos_z):
        Fonte[i] = 0.0 + 0.0j
        Global[i, :] = 0.0 + 0.0j
        Global[i, i] = 1.0 + 1.0j

    for j in range(2 * nos_z, nos, nos_z):
        Fonte[j] = 0.0 + 0.0j
        Global[j, :] = 0.0 + 0.0j
        Global[j, j] = 1.0 + 1.0j

    for i in range(nos_z):
        Global[i, :] = 0.0 + 0.0j
        Global[i, i] = 1.0 + 1.0j

    for j in range(borda, nos):
        Global[j, :] = 0.0 + 0.0j
        Global[j, j] = 1.0 + 1.0j

#===================================================================================
def GlobalMat(e, no, nE, sB, gl, nosE, matE):
    matG = np.zeros((sB, no), dtype=complex)
    for i in range(gl):
        k = nosE[e, i]
        for j in range(gl):
            m = sB - nosE[e, j] + k
            if m < sB:
                n = sB + k - m
                matG[m, n] += matE[j, i]
    return Gmat

#===================================================================================
def matGfactor(matG, nG, sB):
    for j in range(nG):
        q = max(0, j - nG + sB)
        for i in range(sB - 1, q - 1, -1):
            k = j - i + sB
            if abs(matG[i, k]) != 0.0:
                c = matG[i, k] / matG[sB - 1, j]
                m = sB
                n = k - 1
                for p in range(i, q - 1, -1):
                    m -= 1
                    n += 1
                    matG[m, n] -= c * matG[p, n]
                matG[i, k] = c

#===================================================================================
def solve4MagX(matG, nG, sB, vecG):
    for j in range(nG):
        k = max(0, j - nG + sB)
        for m in range(sB - 1, k - 1, -1):
            n = j - m + sB
            vecG[n] -= matG[m, n] * vecG[j]
        vecG[j] /= matG[sB - 1, j]

    for j in range(1, nG):
        i = nG - j
        k = max(0, i - nG + sB)
        for m in range(sB - 1, k - 1, -1):
            n = i - m + sB
            vecG[i] -= matG[m, n] * vecG[n]

#===================================================================================
def soma_hyp(Fonte):
    Hy_z = np.zeros(nos_z, dtype=complex)
    Ex_z = np.zeros(nos_z, dtype=complex)
    n_elements_z = (nos_z - 1) * 2
    j = 0
    for i in range(0, n_elements_z, 2):
        Hy, Ex = ondaplana(i)
        Hy_z[j] = Hy[0]
        Ex_z[j] = Ex[0]
        j += 1
    i = n_elements_z
    Hy, Ex = ondaplana(i)
    Hy_z[j] = Hy[0]
    Ex_z[j] = Ex[0]
    j = 0
    while True:
        for i in range(nos_z):
            Fonte[j] += Hy_z[i]
            j += 1
        if j >= nos:
            break

#===================================================================================
def resistividade_aparente(Fonte):
    Hy = np.zeros(nos_x, dtype=complex)
    Ex = np.zeros(nos_x, dtype=complex)
    j = 0
    for i in range(0, nos, nos_z):
        Ex[j] = -(Fonte[i + 2] - Fonte[i]) / (2.0 * dz) * rho[nc - 1]
        Hy[j] = Fonte[i + 1]
        j += 1
    rho_ap = np.zeros(nos_x, dtype=dbl)
    for i in range(nos_x):
        rho_ap[i] = 1.0 / omega / mi0 * abs(Ex[i] / Hy[i]) ** 2

    with open('hy.txt', 'w') as f:
        for i in range(nos_x):
            f.write(f"{Hy[i].real} {Hy[i].imag}\n")

#===================================================================================
def print_result(Fonte):
    mesh = np.zeros((nos_x, nos_z), dtype=complex)
    n = 0
    for i in range(nos_x):
        for j in range(nos_z):
            mesh[i, j] = Fonte[n + j]
        n += nos_z

    with open('solucao.txt', 'w') as f:
        for i in range(nos_x):
            f.write(" ".join(f"{mesh[i, j].real:.4f}" for j in range(nos_z)) + "\n")

    with open('x.txt', 'w') as f:
        for i in range(nos_x):
            f.write(f"{i * dx}\n")

    with open('z.txt', 'w') as f:
        for i in range(nos_z):
            f.write(f"{i * dz}\n")

    with open('rho_ap.txt', 'w') as f:
        for i in range(nos_x):
            f.write(f"{i * dx - interv_x / 2.0} {rho_ap[i]}\n")

    with open('teste.txt', 'w') as f:
        for i in range(int(nos_x / 2) * nos_z, int(nos_x / 2) * nos_z + nos_z):
            f.write(f"{abs(Fonte[i])}\n")

#===================================================================================
def plota_malha():
    with open('n_elements.txt', 'w') as f:
        f.write(f"{n_elements}\n")
    with open('h_elements.txt', 'w') as f:
        f.write(f"{h_elements}\n")
    with open('malha.txt', 'w') as f:
        for i in range(n_elements):
            f.write(" ".join(f"{nos_elemento[i, j]:7d}" for j in range(3)) + "\n")
    with open('hetero.txt', 'w') as f:
        for i in range(h_elements):
            f.write(" ".join(f"{nos_hetero[i, j]:7d}" for j in range(3)) + "\n")
    with open('coordenadas.txt', 'w') as f:
        for i in range(nos):
            f.write(" ".join(f"{coordenadas[i, j]:10.3f}" for j in range(2)) + "\n")


def main():
    # Step 1: Get data
    get_data()
    print("Data loaded successfully.")

    # Step 2: Get model parameters
    get_model()
    print("Model parameters loaded successfully.")

    # Step 3: Generate the grid
    grid_2d()
    print("Grid generated successfully.")

    # Step 4: Identify heterogeneity
    hetero_ident()
    print("Heterogeneity identified successfully.")

    # Step 5: Solve for Hy and Ex
    Fonte = np.zeros(nos, dtype=complex)
    Global = np.zeros((sB, nos), dtype=complex)
    
    cc_origin(Fonte, Global)
    print("Source and global matrix initialized.")

    cc(Fonte, Global)
    print("Boundary conditions applied.")

    print_result(Fonte)
    print("Results printed successfully.")

    logging.basicConfig(level=logging.INFO)
    logging.info("Starting the FEM simulation.")
    fem_model = FiniteElementModel()
    
    fem_model.load_data()
    logging.info("Data loaded successfully.")
    
    fem_model.generate_grid()
    logging.info("Grid generated successfully.")

    logging.info("Simulation completed.")



if __name__ == "__main__":
    main()

