// Problem Data
// Original version before removal of dynamic memory allocation

#include "cpg_workspace.h" 

// PROBLEM SETUP
// const double A[1][3] = {{1, 1, 0}};
// const double b[1][1] = {{5}};
// const double c[3][1] = {{0},
// 			  	        {0},
// 			  	        {1}};
// const double G[4][4] = {{ 0, -1,  0},
// 			            { 0,  0, -1},
// 			            {-1,  0,  0},
// 			            { 0, -1,  0}};
// const double h[4][1] = {{-2},
// 			            { 0},
// 			            { 0},
// 			            { 0}};

const int num_cones = 2;
const int cone_orders[num_cones] = {1, 3}; // order of each cone
const int cone_exprs[num_cones] = {1, 1}; // number of expressions in each cone
const int num_cone_vars = 4; // sum of each cone's order
// cone_vars calculated in main
// Examples:
// 2*x1 + 2*x2 >= 0 is of order 1 and has 1 expression
// x1>=0, x2>=0 is of order 1 and has 2 expressions
// ||x||_2 >=||{y,z}^T||_2 is of order 3 and has 1 expression

const int maxit = 20; // maximum number of iterations

// #define ROWSIZE(mat) (sizeof(mat) / sizeof(mat[0]));
// #define COLSIZE(mat) (sizeof(mat[0]) / sizeof(mat[0][0]));

// const int nx = COLSIZE(A);
// const int ny = ROWSIZE(A);
// const int nz = COLSIZE(G);
// const int ns = COLSIZE(G);
// const int nt = 1;
// const int nk = 1;

// NEW FORMAT

// Vector containing flattened user-defined parameters
c_float cpg_params_vec[2] = {
(c_float)0.00000000000000000000,
(c_float)1.00000000000000000000,
};

// Sparse mappings from user-defined to canonical parameters
c_int canon_h_map_i[2] = {
0,
1,
};
c_int canon_h_map_p[5] = {
0,
2,
2,
2,
2,
};
c_float canon_h_map_x[2] = {
(c_float)-1.00000000000000000000,
(c_float)-2.00000000000000000000,
};
csc canon_h_map = {2, 4, 2, canon_h_map_p, canon_h_map_i, canon_h_map_x, -1};

// Canonical parameters
c_float canon_c[3] = {
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)1.00000000000000000000,
};

c_float canon_c_ECOS[3] = {
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)1.00000000000000000000,
};

c_int canon_A_i[2] = {
0,
0,
};
c_int canon_A_p[4] = {
0,
1,
2,
2,
};
c_float canon_A_x[2] = {
(c_float)1.00000000000000000000,
(c_float)1.00000000000000000000,
};
csc canon_A = {2, 1, 3, canon_A_p, canon_A_i, canon_A_x, -1};

c_int canon_A_ECOS_i[2] = {
0,
0,
};
c_int canon_A_ECOS_p[4] = {
0,
1,
2,
2,
};
c_float canon_A_ECOS_x[2] = {
(c_float)1.00000000000000000000,
(c_float)1.00000000000000000000,
};
csc canon_A_ECOS = {2, 1, 3, canon_A_ECOS_p, canon_A_ECOS_i, canon_A_ECOS_x, -1};

c_float canon_b[1] = {
(c_float)5.00000000000000000000,
};

c_float canon_b_ECOS[1] = {
(c_float)5.00000000000000000000,
};

c_int canon_G_i[4] = {
2,
0,
3,
1,
};
c_int canon_G_p[4] = {
0,
1,
3,
4,
};
c_float canon_G_x[4] = {
(c_float)-1.19999999999999995559,
(c_float)-1.00000000000000000000,
(c_float)-1.30000000000000004441,
(c_float)-1.10000000000000008882,
};
csc canon_G = {4, 4, 3, canon_G_p, canon_G_i, canon_G_x, -1};

c_int canon_G_ECOS_i[4] = {
2,
0,
3,
1,
};
c_int canon_G_ECOS_p[4] = {
0,
1,
3,
4,
};
c_float canon_G_ECOS_x[4] = {
(c_float)-1.19999999999999995559,
(c_float)-1.00000000000000000000,
(c_float)-1.30000000000000004441,
(c_float)-1.10000000000000008882,
};
csc canon_G_ECOS = {4, 4, 3, canon_G_ECOS_p, canon_G_ECOS_i, canon_G_ECOS_x, -1};

c_float canon_h[4] = {
(c_float)-2.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
};

c_float canon_h_ECOS[4] = {
(c_float)-2.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
};

// Struct containing canonical parameters
Canon_Params_t Canon_Params = {
.c = (c_float *) &canon_c,
.d = 0.00000000000000000000,
.A = &canon_A,
.b = (c_float *) &canon_b,
.G = &canon_G,
.h = (c_float *) &canon_h,
};

Canon_Params_t Canon_Params_ECOS = {
.c = (c_float *) &canon_c_ECOS,
.d = 0.00000000000000000000,
.A = &canon_A_ECOS,
.b = (c_float *) &canon_b_ECOS,
.G = &canon_G_ECOS,
.h = (c_float *) &canon_h_ECOS,
};

// Struct containing flags for outdated canonical parameters
Canon_Outdated_t Canon_Outdated = {
.c = 0,
.d = 0,
.A = 0,
.b = 0,
.G = 0,
.h = 0,
};

// User-defined variables
c_float var1[3] = {
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
(c_float)0.00000000000000000000,
};

// Struct containing primal solution
CPG_Prim_t CPG_Prim = {
.var1 = (c_float *) &var1,
};

// Dual variables associated with user-defined constraints
// Struct containing dual solution
CPG_Dual_t CPG_Dual = {
.d0 = 0,
.d1 = 0,
.d2 = 0,
};

// Struct containing solver info
CPG_Info_t CPG_Info = {
.obj_val = 0,
.iter = 0,
.status = 0,
.pri_res = 0,
.dua_res = 0,
};

// Struct containing solution and info
CPG_Result_t CPG_Result = {
.prim = &CPG_Prim,
.dual = &CPG_Dual,
.info = &CPG_Info,
};

// Struct containing solver settings
Canon_Settings_t Canon_Settings = {
.feastol = 1e-8,
.abstol = 1e-8,
.reltol = 1e-8,
.feastol_inacc = 1e-4,
.abstol_inacc = 5e-5,
.reltol_inacc = 5e-5,
.maxit = 100,
};

// ECOS array of SOC dimensions
c_int ecos_q[1] = {
3,
};

// workspace
pwork* workspace = 0;

// ECOS exit flag
c_int ecos_flag = -99;