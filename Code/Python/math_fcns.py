# Math Functions

import numpy as np


def coneprod(u, v, order):
	# u and v are vectors
	# order is the order of the cone that u and v are in
	if order == 1:
		uv = u*v # Standard multiplication for order 1 cone
	elif order > 1:
		uv = np.block([[u.T@v],
					   [u[0,0]*v[1:] + v[0,0]*u[1:]]])
	else:
		print(f"Order of {order} is not plausible! Setting u o v = 0.")
		uv = u*0 + v*0
	return uv

def coneinv(u, cone_start_idxs, cone_orders, cone_vars):
    num_cones = len(cone_start_idxs)
    uinv = [0]*num_cones
    for i in range(0, num_cones):
        start = cone_start_idxs[i] # current cone starting index
        stop = start + cone_vars[i] # current cone stopping index (uninclusive)
        uk = u[start:stop]
        if cone_orders[i]==1:
            uinv[i] = 1/uk
        elif cone_orders[i]>1:
            lk0 = uk[0,0]
            lk1 = uk[1:]
            tmp1 = np.block([[lk0, -lk1.T],
                             [-lk1, 1/lk0*((lk0**2-lk1.T@lk1)*np.eye(cone_vars[i]-1)+lk1@lk1.T)]])
            uinv[i] = 1/(lk0**2-lk1.T@lk1)*tmp1
        else:
            print(f"Order of {order} is not plausible! Using uinv=0.")
    return uinv