# Optimization Functions

import numpy as np
from scipy.linalg import block_diag
from sys import float_info


# Define machine epsilon
eps = float_info.epsilon

def process_cone_info(coneinfo):
    cone_orders = coneinfo.cone_orders
    cone_exprs = coneinfo.cone_exprs
    cone_vars = coneinfo.cone_vars

    num_cones = len(cone_orders) # == len(cone_vars)
    cone_start_idxs = np.hstack([0, np.cumsum(cone_vars[:-1])])
    # Detect start of SOCs by finding where order is first >1
    soc_idx = np.argmax(cone_orders>1)
    l_idx = np.sum(cone_vars[:soc_idx]) # index corresponding to start of SOCs

    coneinfo.num_cones = num_cones
    coneinfo.cone_start_idxs = cone_start_idxs
    coneinfo.l_idx = l_idx

    # Set up \textbb{e} vector
    e_socs = [0]*num_cones
    for i in range(0, num_cones):
        e_socs[i] = np.ones((cone_vars[i],1))
        if cone_orders[i] > 1:
            # Assume SOC
            e_socs[i][1:,0] = 0
    e = np.vstack(e_socs)
    m = int(e.T@e) # total cone order

    return coneinfo, e, m


def pdinit(A, b, c, G, h, e, coneinfo):
    cone_orders =     coneinfo.cone_orders
    cone_exprs =      coneinfo.cone_exprs
    cone_vars =       coneinfo.cone_vars
    num_cones =       coneinfo.num_cones
    cone_start_idxs = coneinfo.cone_start_idxs
    l_idx =           coneinfo.l_idx

    nx = np.shape(A)[1]
    ny = np.shape(A)[0]
    nz = np.shape(G)[0]
    ns = np.shape(G)[0]
    nt = 1
    nk = 1

    # PRIMAL INITIALIZATION
    init_A_p = np.block([[np.zeros((nx,nx)), A.T, G.T],
                         [A, np.zeros((ny,ny)), np.zeros((ny,ns))],
                         [G, np.zeros((ns,ny)), -np.eye(ns)]])
    init_b_p = np.block([[np.zeros((nx,1))],
                         [b],
                         [h]])
    init_x_p = np.linalg.solve(init_A_p, init_b_p)
    init_x_hat = init_x_p[:nx]
    init_y = init_x_p[nx:nx+ny]
    init_s_tilde = -init_x_p[nx+ny:]

    alpha_p = -99999
    for i in range(0, l_idx):
        start = cone_start_idxs[i] # current cone starting index
        stop = start + cone_vars[i] # current cone stopping index (uninclusive)
        alpha_p = max(alpha_p, -init_s_tilde[start:stop][0])

    for i in cone_start_idxs[l_idx:]:
        tmp_gap = (init_s_tilde[i]-np.linalg.norm(init_s_tilde[i+1:i+cone_vars[i]]))[0]
        alpha_p = max(alpha_p, -tmp_gap)

    if alpha_p < 0:
        init_s_hat = init_s_tilde
    else:
        init_s_hat = init_s_tilde + (1+alpha_p)*e

    # DUAL INITIALIZATION
    init_A_d = np.block([[np.zeros((nx,nx)), A.T, G.T],
                         [A, np.zeros((ny,ny)), np.zeros((ny,nz))],
                         [G, np.zeros((nz,ny)), -np.eye(nz)]])
    init_b_d = np.block([[-c],
                         [np.zeros((ny,1))],
                         [np.zeros((nz,1))]])
    init_x_d = np.linalg.solve(init_A_d, init_b_d)
    init_x = init_x_d[:nx]
    init_y_hat = init_x_d[nx:nx+ny]
    init_z_tilde = init_x_d[nx+ny:]
    alpha_d = -99999
    for i in range(0, l_idx):
        start = cone_start_idxs[i] # current cone starting index
        stop = start + cone_vars[i] # current cone stopping index (uninclusive)
        alpha_d = max(alpha_d, -init_z_tilde[start:stop][0])
        
    for i in cone_start_idxs[l_idx:]:
        tmp_gap = (init_z_tilde[i]-np.linalg.norm(init_z_tilde[i+1:i+cone_vars[i]]))[0]
        alpha_d = max(alpha_d, -tmp_gap)

    if alpha_d < 0:
        init_z_hat = init_z_tilde
    else:
        init_z_hat = init_z_tilde + (1+alpha_d)*e

    return init_x_hat, init_y_hat, init_s_hat, init_z_hat

def find_alpha_lam(lam, del_s_tilde, del_z_tilde, coneinfo, printall=False):
    # Calculate alpha for del_sa (or W^(-T)@del_sa) and del_za (or W@del_za)
    # Find largest step size which has feasible lambda

    cone_vars =       coneinfo.cone_vars
    num_cones =       coneinfo.num_cones
    cone_start_idxs = coneinfo.cone_start_idxs
    l_idx =           coneinfo.l_idx

    # Calculate alpha for del_sc (or W^(-T)@del_sc) and del_zc (or W@del_zc)
    alpha_ks = [0]*num_cones
    for i in range(0, l_idx):
        start = cone_start_idxs[i] # current cone starting index
        stop = start + cone_vars[i] # current cone stopping index (uninclusive)
        del_sck_tilde = del_s_tilde[start:stop]
        del_zck_tilde = del_z_tilde[start:stop]
        lamk = lam[start:stop]

        lamkinv = 1/lamk
        rho_k = lamkinv * del_sck_tilde
        sigma_k = lamkinv * del_zck_tilde
        max_k = max(eps, -min(rho_k), -min(sigma_k))
        if max_k == 0:
            alpha_ks[i] = 1
        else:
            alpha_ks[i] = 1/max_k
    for i in range(l_idx, num_cones):
        start = cone_start_idxs[i] # current cone starting index
        stop = start + cone_vars[i] # current cone stopping index (uninclusive)
        del_sck_tilde = del_s_tilde[start:stop]
        del_zck_tilde = del_z_tilde[start:stop]
        lamk = lam[start:stop]

        J = block_diag(*(1, -np.eye(cone_vars[i]-1)))
        lamk_bar = lamk/np.sqrt(lamk.T@J@lamk)
        tmp1 = np.block([[lamk_bar.T@J@del_sck_tilde],
                         [del_sck_tilde[1:]-(lamk_bar.T@J@del_sck_tilde+del_sck_tilde[0,0])/(lamk_bar[0,0]+1)*lamk_bar[1:]]])
        tmp2 = np.block([[lamk_bar.T@J@del_zck_tilde],
                         [del_zck_tilde[1:]-(lamk_bar.T@J@del_zck_tilde+del_zck_tilde[0,0])/(lamk_bar[0,0]+1)*lamk_bar[1:]]])
        rho_k = 1/np.sqrt(lamk.T@J@lamk) * tmp1
        sigma_k = 1/np.sqrt(lamk.T@J@lamk) * tmp2
        max_k = max(eps, -(rho_k[0,0]-np.linalg.norm(rho_k[1:])),
                    -(sigma_k[0,0]-np.linalg.norm(sigma_k[1:])))
        if max_k == 0:
            alpha_ks[i] = 1
        else:
            alpha_ks[i] = 1/max_k

        if printall:
            print(alpha_ks)

    return min(alpha_ks)

def check_termination(x_hat, y_hat, z_hat, t_hat, k_hat, h, b, c):
    if t_hat==0 and k_hat>0:
        # Then h.T@z + b.T@y + c.T@x < 0, so we must have h.T@z + b.T@y < 0
        # or c.T@x < 0 or both
        if h.T@z_hat + b.T@y_hat + c.T@x_hat < 0:
            print(f"Primal problem is infeasible!")
            # Since G.T@z + A.T@y = 0,  z >= 0,  h.T@z + b.T@y < 0
        elif h.T@z_hat + b.T@y_hat < 0:
            print(f"Dual problem is infeasible!")
            # Since G@x + s = 0,  Ax = 0,  s >= 0,  c.T@x < 0
    elif t_hat==0 and k_hat==0:
        print(f"t_hat and k_hat are both zero!")
        print(f"No conclusion can be made about the results!")
