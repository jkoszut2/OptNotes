import numpy as np
from scipy.linalg import block_diag
from math_fcns import coneprod, coneinv
from prob_data import A, b, c, G, h, cone_orders, cone_exprs, cone_vars, maxit
from opt_fcns import process_cone_info, pdinit, find_alpha_lam, check_termination

# Extract size information
nx = np.shape(A)[1]
ny = np.shape(A)[0]
nz = np.shape(G)[0]
ns = np.shape(G)[0]
nt = 1
nk = 1

# Cone information class
class ConeInfo(object):
	pass

# Initialize cone info object
coneinfo = ConeInfo()
coneinfo.cone_orders = cone_orders
coneinfo.cone_exprs = cone_exprs
coneinfo.cone_vars = cone_vars

# Extract additional cone info from given data
coneinfo, e, m = process_cone_info(coneinfo)
num_cones = coneinfo.num_cones
cone_start_idxs = coneinfo.cone_start_idxs
l_idx = coneinfo.l_idx

# Obtain initial values (assuming now that they are not given)
init_x_hat, init_y_hat, init_s_hat, init_z_hat = pdinit(A, b, c, G, h, e, coneinfo)
init_k_hat, init_t_hat = [1, 1]
mu_hat = (init_s_hat.T@init_z_hat+init_k_hat*init_t_hat)/(m+1) # s_hat.T@z_hat = lam.T@lam

# For record keeping
x_hat_hist = np.zeros((nx, maxit))
s_hat_hist = np.zeros((ns, maxit))
y_hat_hist = np.zeros((ny, maxit))
z_hat_hist = np.zeros((nz, maxit))
k_hat_hist = np.zeros((nk, maxit))
t_hat_hist = np.zeros((nt, maxit))
mu_hat_hist = np.zeros((1, maxit))
w_bar_hist = np.zeros((ns, maxit))

# STEP 1 - EVALUATE RESIDUALS, GAP, AND STOPPING CRITERIA
x_hat = init_x_hat
s_hat = init_s_hat
y_hat = init_y_hat
z_hat = init_z_hat
k_hat = init_k_hat
t_hat = init_t_hat
x_hat_hist[:,[0]] = x_hat # [] around k to extract 2D col as opposed to 1D array
s_hat_hist[:,[0]] = s_hat
y_hat_hist[:,[0]] = y_hat
z_hat_hist[:,[0]] = z_hat
k_hat_hist[:,[0]] = k_hat
t_hat_hist[:,[0]] = t_hat
mu_hat_hist[:,[0]] = mu_hat
# Residuals
# resx0 = max(1.0, np.sqrt(c.T@c))
# resy0 = max(1.0, np.sqrt(b.T@b))
# resz0 = max(1.0, misc.snrm2(h, dims))

# Print info header
print(f"  i|  alphaa|  alphac|       k|     t|   pobj|   dobj")

k=0
while k <(maxit-1) and k_hat>1e-7:
	x_hat = x_hat_hist[:,[k]]
	s_hat = s_hat_hist[:,[k]]
	y_hat = y_hat_hist[:,[k]]
	z_hat = z_hat_hist[:,[k]]
	k_hat = k_hat_hist[:,[k]]
	t_hat = t_hat_hist[:,[k]]
	mu_hat = mu_hat_hist[:,[k]]
	v1 = np.block([[np.zeros((nx, 1))],
				   [np.zeros((ny, 1))],
				   [s_hat],
				   [k_hat]])
	m1 = np.block([[np.zeros((nx,nx)),        A.T       ,         G.T      , c],
				   [       -A        , np.zeros((ny,ny)), np.zeros((ny,nz)), b],
				   [       -G        , np.zeros((nz,ny)), np.zeros((nz,nz)), h],
				   [     -c.T        ,       -b.T       ,        -h.T      , 0]])
	v2 = np.block([[x_hat_hist[:,[k]]],
				   [y_hat_hist[:,[k]]],
				   [z_hat_hist[:,[k]]],
				   [t_hat_hist[:,[k]]]])
	res = v1 - m1@v2
	rx = res[:nx]
	ry = res[nx:nx+ny]
	rz = res[nx+ny:nx+ny+nz]
	rt = res[nx+ny+nz:nx+ny+nz+nt]

	if k == 0:
		# Compute Nesterov-Todd scaling W at s_hat and z_hat
		# For LP cone Q1, W = sqrt(s/z)
		W_lp = [0]*l_idx
		for i in range(0, l_idx):
			start = cone_start_idxs[i] # current cone starting index
			stop = start + cone_vars[i] # current cone stopping index (uninclusive)
			W_lp[i] = np.sqrt(s_hat[start:stop]/z_hat[start:stop])
			w_bar_hist[start:stop,[k]] = np.nan # w_bar is not defined for LP case
		# For SOC cone Qm, W = multiple steps
		W_so = [0]*(num_cones-l_idx)
		for i in range(l_idx, num_cones):
			start = cone_start_idxs[i] # current cone starting index
			stop = start + cone_vars[i] # current cone stopping index (uninclusive)

			# Extract current z and s values
			zk = z_hat[start:stop]
			sk = s_hat[start:stop]

			# Calculate z_bar and s_bar
			tmp_zfrac = (zk[0]**2 - zk[1:].T@zk[1:])
			tmp_sfrac = (sk[0]**2 - sk[1:].T@sk[1:])
			zk_bar = zk/np.sqrt(tmp_zfrac)
			sk_bar = sk/np.sqrt(tmp_sfrac)

			# Calculate gamma and w
			gammak = np.sqrt((1+zk_bar.T@sk_bar)/2)
			J = block_diag(*(1, -np.eye(cone_vars[i]-1)))
			wk_bar = 1/(2*gammak)*(sk_bar+J@zk_bar)
			w_bar_hist[start:stop,[k]] = wk_bar

			# Calculate eta and W_so
			curr_eta = (tmp_sfrac/tmp_zfrac)**0.25
			tmp3 = np.eye(cone_vars[i]-1) + 1/(1+wk_bar[0]) * wk_bar[1:]@wk_bar[1:].T
			W_so[i-l_idx] = curr_eta*np.block([[wk_bar[0] , wk_bar[1:].T],
											   [wk_bar[1:], tmp3]])

		W = block_diag(*(W_lp+W_so)) # W = blkdiag(W1, W2, ..., WN)
		lam = W@z_hat # Should equal W^(-T)@s_hat or just W^(-1)@s_hat if W=W.T
		laminv = coneinv(lam, cone_start_idxs, cone_orders, cone_vars)
		Winv = np.linalg.inv(W)
		lam_check = Winv@s_hat


	# STEP 1.5 CHECK TERMINATION CRITERIA
	check_termination(x_hat, y_hat, z_hat, t_hat, k_hat, h, b, c)

	# STEP 2 - SOLVE FOR THE AFFINE DIRECTION
	KKT = np.block([[np.zeros((nx,nx)),         A.T      ,         G.T      ],
					[       -A        , np.zeros((ny,ny)), np.zeros((ny,nz))],
					[       -G        , np.zeros((nz,ny)),        W.T@W     ]])
	v1a = np.block([[-c],
					[-b],
					[-h]])
	v2a = np.block([[rx],
					[ry],
					[rz-s_hat]])
	tmp_sol1 = np.linalg.solve(KKT, v1a)
	x1a = tmp_sol1[:nx]
	y1a = tmp_sol1[nx:nx+ny]
	z1a = tmp_sol1[nx+ny:nx+ny+nz]
	tmp_sol2 = np.linalg.solve(KKT, v2a)
	x2a = tmp_sol2[:nx]
	y2a = tmp_sol2[nx:nx+ny]
	z2a = tmp_sol2[nx+ny:nx+ny+nz]
	dta = rt
	dka = k_hat*t_hat
	del_ta = (dta - dka/t_hat + c.T@x2a + b.T@y2a + h.T@z2a) \
			/ (k_hat/t_hat - c.T@x1a - b.T@y1a - h.T@z1a)
	del_xa = x2a + del_ta*x1a
	del_ya = y2a + del_ta*y1a
	del_za = z2a + del_ta*z1a
	del_sa = -W.T@lam - W.T@W@del_za
	del_ka = -(dka + k_hat*del_ta)/t_hat
	del_sa_tilde = Winv.T@del_sa
	del_za_tilde = W@del_za

	# STEP 3 - SOLVE FOR THE STEP SIZE AND CENTERING PARAMETER
	alpha = 1
	# Calculate alpha for del_ta
	if del_ta < 0:
		alpha = min(alpha, -t_hat/del_ta)
	# Calculate alpha for del_ka
	if del_ka < 0:
		alpha = min(alpha, -k_hat/del_ka)
	alphaa_lam = find_alpha_lam(lam, del_sa_tilde, del_za_tilde, coneinfo)
	alphaa = float(min(alpha, alphaa_lam))
	sigma = (1-alphaa)**3

	# STEP 4 - SOLVE FOR THE COMBINED DIRECTION
	# Very similar to STEP 3. Reuse some of the material from there.
	# Have to manually calculate ds term by splitting up vectors into cones
	# ds = lam O lam + (W^(-T)@del_sa) O (W@del_za) + sigma*mu_hat*e
	ds = np.zeros((sum(cone_vars),1))

	# supp_term is -W.T@(lam inv ds) term in (45b) and (47) CVXOPT
	supp_term = np.zeros((sum(cone_vars),1))
	for i in range(0, num_cones):
		start = cone_start_idxs[i] # current cone starting index
		stop = start + cone_vars[i] # current cone stopping index (uninclusive)

		del_sak = del_sa[start:stop]
		del_zak = del_za[start:stop]
		lamk = lam[start:stop]
		ek = e[start:stop]
		Wk = W[start:stop,start:stop]
		Winvk = Winv[start:stop,start:stop]

		ds[start:stop] = coneprod(lamk,lamk,cone_orders[i]) \
						  + coneprod(Winvk.T@del_sak,Wk@del_zak,cone_orders[i]) \
						  - sigma*mu_hat*ek
		supp_term[start:stop] = -Wk.T@(laminv[i]@ds[start:stop])

	# Solution of first KKT system stays the same
	x1c = x1a
	y1c = y1a
	z1c = z1a
	# Solution of second KKT system changes
	v2c = np.block([[(1-sigma)*rx],
					[(1-sigma)*ry],
					[(1-sigma)*rz+supp_term]])
	tmp_sol2c = np.linalg.solve(KKT, v2c)
	x2c = tmp_sol2c[:nx]
	y2c = tmp_sol2c[nx:nx+ny]
	z2c = tmp_sol2c[nx+ny:nx+ny+nz]

	dtc = (1-sigma)*rt
	dkc = k_hat*t_hat + del_ka*del_ta - sigma*mu_hat
	del_tc = (dtc - dkc/t_hat + c.T@x2c + b.T@y2c + h.T@z2c) \
			 / (k_hat/t_hat - c.T@x1c - b.T@y1c - h.T@z1c)
	del_xc = x2c + del_tc*x1c
	del_yc = y2c + del_tc*y1c
	del_zc = z2c + del_tc*z1c
	del_sc = np.zeros((sum(cone_vars),1))
	for i in range(0, num_cones):
		a = cone_start_idxs[i] # current cone starting index
		o = a + cone_vars[i] # current cone stopping index (uninclusive)
		del_sc[a:o] = -W[a:o,a:o].T@(laminv[i]@ds[a:o]+W[a:o,a:o]@del_zc[a:o])

	del_kc = -(dkc + k_hat*del_tc)/t_hat
	del_sc_tilde = Winv.T@del_sc
	del_zc_tilde = W@del_zc

	# STEP 5 - SOLVE FOR THE STEP SIZE AND UPDATE ITERATES AND SCALING MATRICES
	# Solve for step size
	alpha = 1
	# Logic below only works if t_hat and k_hat are >= 0 which is how the
	# algorithm is defined
	# Calculate alpha for del_ta
	if del_tc < 0:
		alpha = min(alpha, -t_hat/del_tc)
	# Calculate alpha for del_ka
	if del_kc < 0:
		alpha = min(alpha, -k_hat/del_kc)
	# Calculate alpha for del_sc (or W^(-T)@del_sc) and del_zc (or W@del_zc)
	alphac_lam = find_alpha_lam(lam, del_sc_tilde, del_zc_tilde, coneinfo)
	alphac = float(min(alpha, alphac_lam))*0.98

	# Update iterates
	s_hat = s_hat + alphac*del_sc
	k_hat = k_hat + alphac*del_kc
	x_hat = x_hat + alphac*del_xc
	y_hat = y_hat + alphac*del_yc
	z_hat = z_hat + alphac*del_zc
	t_hat = t_hat + alphac*del_tc
	mu_hat = (s_hat.T@z_hat+k_hat*t_hat)/(m+1)
	# For record keeping
	x_hat_hist[:,[k+1]] = x_hat # [] around k to extract 2D col as opposed to 1D array
	s_hat_hist[:,[k+1]] = s_hat
	y_hat_hist[:,[k+1]] = y_hat
	z_hat_hist[:,[k+1]] = z_hat
	k_hat_hist[:,[k+1]] = k_hat
	t_hat_hist[:,[k+1]] = t_hat
	mu_hat_hist[:,[k+1]] = mu_hat

	# Print info
	pobj = c.T@x_hat
	dobj = -b.T@y_hat - h.T@z_hat
	print(f"{k:>3}|{alphaa:>8.3}|{alphac:>8.3}|{k_hat[0,0]:>8.3}|" + \
		  f"{t_hat[0,0]:>6.3}|{pobj[0,0]:>7.3}|{dobj[0,0]:>7.3}|")

	# Update scaling matrices
	# Compute Nesterov-Todd scaling W at s_hat and z_hat
	# For LP cone Q1, W = sqrt(s/z)
	W_lp_next = [0]*l_idx
	lam_lp_next = [0]*l_idx
	s_tilde_next = lam + alphac*Winv.T@del_sc
	z_tilde_next = lam + alphac*W@del_zc
	for i in range(0, l_idx):
		start = cone_start_idxs[i] # current cone starting index
		stop = start + cone_vars[i] # current cone stopping index (uninclusive)
		wk = np.diagonal(W_lp[i]).reshape(-1,1) # This is only true for the LP case!
		W_lp_next[i] = np.sqrt(s_tilde_next[start:stop]/z_tilde_next[start:stop]) * wk
		lam_lp_next[i] = np.sqrt(s_tilde_next[start:stop]*z_tilde_next[start:stop])
		w_bar_hist[start:stop,[k+1]] = np.nan # w_bar is not defined for LP case
	# For SOC cone Qm, W = multiple steps
	W_so_next = [0]*(num_cones-l_idx)
	lam_so_next = [0]*(num_cones-l_idx)
	for i in range(l_idx, num_cones):
		start = cone_start_idxs[i] # current cone starting index
		stop = start + cone_vars[i] # current cone stopping index (uninclusive)

		# Extract current z_hat and s_hat values
		zk_hat = z_hat_hist[start:stop,[k]]
		sk_hat = s_hat_hist[start:stop,[k]]
		zk_tilde_next = z_tilde_next[start:stop]
		sk_tilde_next = s_tilde_next[start:stop]

		# For testing
		zk_hat_next = z_hat_hist[start:stop,[k+1]]
		sk_hat_next = s_hat_hist[start:stop,[k+1]]
		tmp_z_hatfrac = (zk_hat_next[0]**2 - zk_hat_next[1:].T@zk_hat_next[1:])
		tmp_s_hatfrac = (sk_hat_next[0]**2 - sk_hat_next[1:].T@sk_hat_next[1:])
		zk_hat_bar_next = zk_hat_next/np.sqrt(tmp_z_hatfrac)
		sk_hat_bar_next = sk_hat_next/np.sqrt(tmp_s_hatfrac)

		# Calculate z_tilde_bar_next and s_tilde_bar_next
		tmp_zfrac = (zk_tilde_next[0]**2 - zk_tilde_next[1:].T@zk_tilde_next[1:])
		tmp_sfrac = (sk_tilde_next[0]**2 - sk_tilde_next[1:].T@sk_tilde_next[1:])
		zk_tilde_bar_next = zk_tilde_next/np.sqrt(tmp_zfrac)
		sk_tilde_bar_next = sk_tilde_next/np.sqrt(tmp_sfrac)

		# Calculate gamma and w
		gammak_next = np.sqrt((1+zk_tilde_bar_next.T@sk_tilde_bar_next)/2)
		# above is same as np.sqrt((1+zk_hat_bar_next.T@zk_hat_bar_next)/2)
		J = block_diag(*(1, -np.eye(cone_vars[i]-1)))
		wk_bar_next = 1/(2*gammak_next)*(sk_hat_bar_next+J@zk_hat_bar_next)
		w_bar_hist[start:stop,[k+1]] = wk_bar_next # w_bar is not defined for LP case

		# Calculate eta and W_so
		curr_eta = (tmp_s_hatfrac/tmp_z_hatfrac)**0.25
		tmp3 = np.eye(cone_vars[i]-1) + 1/(1+wk_bar_next[0]) * wk_bar_next[1:]@wk_bar_next[1:].T
		W_so_next[i-l_idx] = curr_eta*np.block([[wk_bar_next[0] , wk_bar_next[1:].T],
										   [wk_bar_next[1:], tmp3]])

		# NEW
		ek = e[start:stop]
		qk_bar = 1/(2*gammak_next)*(sk_tilde_bar_next+J@zk_tilde_bar_next)
		wk_bar = w_bar_hist[start:stop,[k]]
		vk = 1/np.sqrt(2*(wk_bar[0,0]+1))*(wk_bar+ek)

		# Calculate lam
		lamk0_bar_next = gammak_next
		uk = sk_tilde_bar_next - J@zk_tilde_bar_next
		d = (vk[0,0]*(vk.T@uk) - uk[0,0]/2)/(2*vk[0,0]*vk.T@qk_bar - qk_bar[0,0]+1)
		Wk_bar = 2*vk@vk.T - J
		# Might need sk_tilde_bar instead of sk_tilde_bar_next (same for z)
		lamk1_bar_next = Wk_bar@((1-d/gammak_next)/2*sk_tilde_bar_next-(1+d/gammak_next)/2*J@zk_tilde_bar_next)
		lamk_bar_next = np.block([[lamk0_bar_next],
								  [lamk1_bar_next[1:]]])
		lam_so_next[i-l_idx] = np.sqrt(np.sqrt(sk_tilde_next.T@J@sk_tilde_next) * \
									   np.sqrt(zk_tilde_next.T@J@zk_tilde_next)) * lamk_bar_next

	# W = blkdiag(W1, W2, ..., WN)
	# embed()
	W = block_diag(*(W_lp_next+W_so_next))
	Winv = np.linalg.inv(W)
	lam = np.vstack(lam_lp_next+lam_so_next)
	laminv = coneinv(lam, cone_start_idxs, cone_orders, cone_vars)

	# Repeat (go back to step 1.5 where termination criteria is checked)
	k += 1


print("\nSolution History")
print(np.round(x_hat_hist[:,:k+1].T/t_hat_hist[:,:k+1].T, 5))
print("\n")