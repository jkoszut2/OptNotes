#include <stdio.h>
#include <stdlib.h>
#include "prob_workspace.h"
#include "ecos.h"
#include "math_fcns.h"
#include "opt_fcns.h"
#include "amd.h"
#include "amd_internal.h"
#include "splamm.h"
#include "ldl.h"


static c_int initialized = 0;

int main() {

	int cone_vars[num_cones]; // number of variables in each cone
	for (int i=0; i<(sizeof(cone_vars)/sizeof(int)); i++) {
	    cone_vars[i] = cone_orders[i] * cone_exprs[i];
	}
	printf("Value for cone_vars[%d]: is %d\n", 1, cone_vars[1]);

    const int n = 3; // = sizeof(u)/sizeof(double)
	double u[n] = {3, 2, 1};
	double v[n] = {4, 3, -1};
	double uv[n] = {0};
    coneprod2(u, v, sizeof(u)/sizeof(double), uv);

	for (int i=0; i<n; i++) {
	    printf("Value for uv[%d]: is %f\n", i, uv[i]);
	}

	printf("Value for coneinfo.cone_vars[1]: is %d\n", cone_vars[1]);

	printf("Machine epsilon is %f\n", DBL_EPS);

	int soc_idx = 0;                 // cone index corresponding to start of SOCs
	int l_idx = 0;                   // variable index corresponding to start of SOCs
	int m = 0;                       // total cone order
	process_cone_info(cone_vars, cone_start_idxs, &l_idx, e_socs, &m);

	for (int i=0; i<num_cone_vars; i++) {
	    printf("Value for e_socs[%d]: is %d\n", i, e_socs[i]);
	}
	printf("Value for m is is %d\n", m);


 	if (!initialized) {
 		// ECOS_setup(idxint n, idxint m, idxint p, idxint l, idxint ncones,
 				   // idxint* q, idxint nexc,
                   // pfloat* Gpr, idxint* Gjc, idxint* Gir,
                   // pfloat* Apr, idxint* Ajc, idxint* Air,
                   // pfloat* c, pfloat* h, pfloat* b)
 		// Canon_Params_ECOS.G->x[1] is already dereferenced
 		// *Canon_Params_ECOS.G->x gives same as Canon_Params_ECOS.G->x[0]
 		printf("Value for Canon_Params_ECOS.G->x[1] is %f\n", Canon_Params_ECOS.G->x[1]);
 		// printf("Memory address for setupspace is %p\n", setupspace);
		// setupspace = ECOS_setup(3, 4, 1, 1, 1, (int *) &ecos_q, 0, Canon_Params_ECOS.G->x,
								  // Canon_Params_ECOS.G->p, Canon_Params_ECOS.G->i,
								  // Canon_Params_ECOS.A->x, Canon_Params_ECOS.A->p,
								  // Canon_Params_ECOS.A->i, Canon_Params_ECOS.c,
								  // Canon_Params_ECOS.h, Canon_Params_ECOS.b);

 		// ++Start++ of ECOS_setup (should be replaced by codegen)
 		// Parts have been moved to "PROCESSED" PROBLEM DATA STARTS BELOW section in prob_data.h
 		const idxint n = 3; // number of variables
		const idxint m = 4; // number of inequality constraints
		const idxint p = 1; // number of equality constraints
		const idxint l = 1; // dimension of lp cone (positive orthant)
		const idxint ncones = 1; // number of so cones
		const idxint* q = (int *) &ecos_q; // dim of each so cone
		const int soc0size = 3; // from codegen - dim of 0-th so cone
		idxint nexc = 0; // number of exponential cones
 		pfloat* Gpr = Canon_Params_ECOS.G->x; // each data point
 		idxint* Gjc = Canon_Params_ECOS.G->p; // cumsum elements per col
		idxint* Gir = Canon_Params_ECOS.G->i; // row of each data point
		pfloat* Apr = Canon_Params_ECOS.A->x;
		idxint* Ajc = Canon_Params_ECOS.A->p;
		idxint* Air = Canon_Params_ECOS.A->i;
		pfloat* c = Canon_Params_ECOS.c;
		pfloat* h = Canon_Params_ECOS.h;
		pfloat* b = Canon_Params_ECOS.b;

		PRINTTEXT("\n");
	    PRINTTEXT("    Primal variables (n): %d\n", (int)n);
		PRINTTEXT("Equality constraints (p): %d\n", (int)p);
		PRINTTEXT("     Conic variables (m): %d\n", (int)m);
	    PRINTTEXT("         Size of LP cone: %d\n", (int)l);
	    PRINTTEXT("          Number of SOCs: %d\n", (int)ncones);
	    for(int i=0; i<ncones; i++ ){
	        PRINTTEXT("    Size of SOC #%02d: %d\n", (int)(i+1), (int)q[i]);
	    }
	    PRINTTEXT("\n");

		// // Allocate memory for work data structure
		static pwork setupspace;

		// // Set dimensions
		setupspace.n = n; // number of primal variables (x)
		setupspace.m = m; // number of inequality constraints (s,z)
		setupspace.p = p; // number of equality constraints (y)
	    setupspace.D = l + ncones; // number of cone constraints

		// Allocate memory for problem variables
		// static to make sure variables exist outside scope of init block
	    static pfloat setupspace_x[n];
	    static pfloat setupspace_y[m];
	    static pfloat setupspace_z[m];
	    static pfloat setupspace_s[m];
	    static pfloat setupspace_lambda[m];
	    static pfloat setupspace_dsaff_by_W[m];
	    static pfloat setupspace_dsaff[m];
	    static pfloat setupspace_dzaff[m];
	    static pfloat setupspace_saff[m];
	    static pfloat setupspace_zaff[m];
		static pfloat setupspace_W_times_dzaff[m];
	    setupspace.x =             setupspace_x;
	    setupspace.y =             setupspace_y;
	    setupspace.z =             setupspace_z;
	    setupspace.s =             setupspace_s;
	    setupspace.lambda =        setupspace_lambda;
	    setupspace.dsaff_by_W =    setupspace_dsaff_by_W;
	    setupspace.dsaff =         setupspace_dsaff;
	    setupspace.dzaff =         setupspace_dzaff;
	    setupspace.saff =          setupspace_saff;
	    setupspace.zaff =          setupspace_zaff;
		setupspace.W_times_dzaff = setupspace_W_times_dzaff;

	    // Allocate memory for best iterates so far
  	    // If don't use static for arrays, then they lose scope after init block
  	    // and can run into memory acccess issue (seg fault) during main IPM loop
	    static pfloat setupspace_best_x[n];
	    static pfloat setupspace_best_y[p];
	    static pfloat setupspace_best_z[m];
	    static pfloat setupspace_best_s[m];
	    static stats* setupspace_best_info;
	    setupspace.best_x =    setupspace_best_x;
	    setupspace.best_y =    setupspace_best_y;
	    setupspace.best_z =    setupspace_best_z;
	    setupspace.best_s =    setupspace_best_s;
	    setupspace.best_info = setupspace_best_info;

		// Allocate memory for cone struct
		static cone setupspace_C;

		// Allocate memory for LP cone
	    static lpcone C_lpc;
	    static pfloat lpc_w[l];
	    static pfloat lpc_v[l];
	    static idxint lpc_kkt_idx[l];
	    C_lpc.p = l;
		C_lpc.w = lpc_w;
		C_lpc.v = lpc_v;
		C_lpc.kkt_idx = lpc_kkt_idx;
		setupspace_C.lpc = &C_lpc;

	    static socone C_soc[ncones];
	    idxint cidx = 0, conesize = 0;
	    // Codegen will need to produce below for all SOCs
	    // Can't use for loop unless perform dynamic memory allocation
	    // Right now only have one SOC so it is simple
	    // for(int i=0; i<ncones; i++ ){
        conesize = (idxint)q[0];
        C_soc[0].p = conesize;
        C_soc[0].a = 0;
		C_soc[0].eta = 0;
        static pfloat soc0_q[soc0size-1];
        static pfloat soc0_skbar[soc0size];
        static pfloat soc0_zkbar[soc0size];
        static idxint soc0_colstart[soc0size];
        C_soc[0].q = soc0_q;
		C_soc[0].skbar = soc0_skbar;
		C_soc[0].zkbar = soc0_zkbar;
        C_soc[0].colstart = soc0_colstart;
        cidx += conesize;
	    // }
	    setupspace_C.nsoc = ncones;
	    setupspace_C.soc = C_soc;
		setupspace.C = &setupspace_C;

		printf("cidx is %d\n", cidx);
		printf("l is %d\n", l);
		printf("m is %d\n", m);
	    /* The number of conic variables has to equal l+sum(q) else terminate*/
	    pwork* setupspace_ptr = &setupspace;
	    if(cidx+l!=m)
	    {
	        PRINTTEXT("Number of conic variables does not match l+sum(q)\n");
	        setupspace_ptr = NULL;
	    }

		// Allocate memory for info struct
		static stats setupspace_stats = {0};
	    setupspace.info = &setupspace_stats;

		/* settings */
		static settings setupspace_stgs = {0};
		setupspace.stgs = &setupspace_stgs;
		setupspace.stgs->maxit = MAXIT;
		setupspace.stgs->gamma = GAMMA;
		setupspace.stgs->delta = DELTA;
	    setupspace.stgs->eps = EPS;
		setupspace.stgs->nitref = NITREF;
		setupspace.stgs->abstol = ABSTOL;
		setupspace.stgs->feastol = FEASTOL;
		setupspace.stgs->reltol = RELTOL;
	    setupspace.stgs->abstol_inacc = ATOL_INACC;
		setupspace.stgs->feastol_inacc = FTOL_INACC;
		setupspace.stgs->reltol_inacc = RTOL_INACC;
	    setupspace.stgs->verbose = VERBOSE;

	    // Written settings
	    setupspace.c = c;
	    setupspace.h = h;
	    setupspace.b = b;

	    /* Store problem data */
	    // n is equal to num primal vars which is = to num of cols in A and G
	    // indexing last col index gives total num of vars in "x" vector
	    // Ajc is "col num nnz lookup reference"
	    // e.g. Ajc[3]-Ajc[2] = num nnz in col index 2
	    // Ajc[n] is total num of elements
	    // Air is "row index lookup
	    // Apr is the data
  
	    // Set up A matrix (sparse matrix)
	    // setupspace->A = ecoscreateSparseMatrix(p, n, Ajc[n], Ajc, Air, Apr);
  	    spmat tmpMatA;
	    tmpMatA.m = p;
	    tmpMatA.n = n;
	    tmpMatA.nnz = Ajc[n];	
	    tmpMatA.jc = Ajc;
        tmpMatA.ir = Air;
        tmpMatA.pr = Apr;
	    if (tmpMatA.jc)
	    	tmpMatA.jc[n] = Ajc[n];
        setupspace.A = &tmpMatA;

	    // Set up G matrix (sparse matrix)
	    // setupspace->G = ecoscreateSparseMatrix(m, n, Gjc[n], Gjc, Gir, Gpr);
  	    spmat tmpMatG;
	    tmpMatG.m = m;
	    tmpMatG.n = n;    
	    tmpMatG.nnz = Gjc[n];	
	    tmpMatG.jc = Gjc;
        tmpMatG.ir = Gir;
        tmpMatG.pr = Gpr;
	    if (tmpMatG.jc)
	    	tmpMatG.jc[n] = Gjc[n];
        setupspace.G = &tmpMatG;

	    // Transpose A
	    idxint AtoAt[setupspace.A->nnz];
		spmat* A = setupspace.A;
		// spmat* At = newSparseMatrix(A->n, A->m, A->nnz);
  	    spmat tmpAt;
  	    static idxint At_jc[p+1]; // can't use A->m+1 as index with static
  	    static idxint At_ir[n];
  	    static pfloat At_pr[n];
	    tmpAt.m = A->n;
	    tmpAt.n = A->m;
	    tmpAt.nnz = A->nnz;	
	    tmpAt.jc = At_jc;
        tmpAt.ir = At_ir;
        tmpAt.pr = At_pr;
	    if (tmpAt.jc)
	    	tmpAt.jc[A->m] = A->nnz;
        spmat* At = &tmpAt;
		if (A->nnz != 0) // Assumed true but kept bc keeps w and q only local scope
		{
			// Using m instead of A->m enables initialization
			idxint w[p] = {0}; // array of size rows(A)

			/* row count: how often does row k occur in A? */
			// w[row_idx] gives num of nnz's in row_idx
			for(int k=0; k < A->nnz; k++ ) { w[A->ir[k]]++; }

			// Row becomes column for transposed matrix so can use row to calc jc
			/* row pointers: cumulative sum of w gives At->jc */
			// spla_cumsum(At->jc, w, A->m);
			idxint cumsum = 0;
			for(int k=0; k<A->m; k++) {
				At->jc[k] = cumsum;
				cumsum += w[k];
				w[k] = At->jc[k];
			}
			for(int i=0; i<A->n+1; i++){printf("A->jc[%d] is %d\n", i, A->jc[i]);}
			for(int i=0; i<A->m+1; i++){printf("At->jc[%d] is %d\n", i, At->jc[i]);}
			for(int i=0; i<A->m; i++){printf("w[%d] is %d\n", i, w[i]);}

			/* now walk through A and copy data to right places and set row counter */
			idxint q;
			// for each col j in A
			for(int j=0; j < A->n; j++ ){
				// for each absolute elem k in the j-th col of A
				for(int k = A->jc[j]; k < A->jc[j+1]; k++ ){
					// A->ir[k] is row of k-th data point
					// w[idx] is current elem number in idx-th col
					// w[idx] gets incremented after q is assigned to it
					q = w[A->ir[k]]++;
					At->ir[q] = j; // row of q-th data point is j
					At->pr[q] = A->pr[k]; // q-th pt of Gt is k-th pt of G
					AtoAt[k] = q; // same note as above
				}
			}
		}

	    // Transpose G
		idxint GtoGt[setupspace.G->nnz];
		spmat* G = setupspace.G;
		// spmat* Gt = newSparseMatrix(G->n, G->m, G->nnz);
  	    spmat tmpGt;
  	    static idxint Gt_jc[m+1];
  	    static idxint Gt_ir[n];
  	    static pfloat Gt_pr[n];
	    tmpGt.m = G->n;
	    tmpGt.n = G->m;
	    tmpGt.nnz = G->nnz;	
	    tmpGt.jc = Gt_jc;
        tmpGt.ir = Gt_ir;
        tmpGt.pr = Gt_pr;
	    if (tmpGt.jc)
	    	tmpGt.jc[G->m] = G->nnz;
        spmat* Gt = &tmpGt;

		if (G->nnz != 0) // Assumed true but kept bc keeps w and q only local scope
		{
			// Using m instead of G->m enables initialization
			idxint w[m] = {0}; // array of size rows(G)

			/* row count: how often does row k occur in G? */
			for(int k=0; k < G->nnz; k++ ) { w[G->ir[k]]++; }

			// Row becomes column for transposed matrix so can use row to calc jc
			/* row pointers: cumulative sum of w gives Gt->jc */
			// spla_cumsum(Gt->jc, w, G->m);
			idxint cumsum = 0;
			for(int k=0; k<G->m; k++) {
				Gt->jc[k] = cumsum;
				cumsum += w[k];
				w[k] = Gt->jc[k];
			}
			for(int i=0; i<G->n+1; i++){printf("G->jc[%d] is %d\n", i, G->jc[i]);}
			for(int i=0; i<G->m+1; i++){printf("Gt->jc[%d] is %d\n", i, Gt->jc[i]);}
			for(int i=0; i<G->m; i++){printf("w[%d] is %d\n", i, w[i]);}

			/* now walk through G and copy data to right places and set row counter */
			idxint q;
			// for each col j in G
			for(int j=0; j < G->n; j++ ){
				// for each absolute elem k in the j-th col of G
				for(int k = G->jc[j]; k < G->jc[j+1]; k++ ){
					// G->ir[k] is row of k-th data point
					// w[idx] is current elem number in idx-th col
					// w[idx] gets incremented after q is assigned to it
					q = w[G->ir[k]]++;
					Gt->ir[q] = j; // row of q-th data point is j
					Gt->pr[q] = G->pr[k]; // q-th pt of Gt is k-th pt of G
					GtoGt[k] = q; // same note as above
				}
			}
		}

	    // Set up (upper part of) KKT system
		/* nK = Dimension of KKT matrix
	     *    =  n (number of variables)
	     *     + p (number of equality constraints)
	     *     + m (number of inequality constraints)
	     */
	    static const idxint nK = n+m+p;
		idxint AttoK[setupspace.A->nnz];
		idxint GttoK[setupspace.G->nnz];
	    idxint Sign[nK];
		spmat *KU; // &KU is a pointer to a pointer
		// createKKT_U(Gt, At, setupspace.C, &Sign, &KU, AttoK, GttoK);

	    PRINTTEXT("Non-zeros in KKT matrix: %d\n", (int) NNZK);

        // ++Start++ of KKTUpper setup
		/* Allocate memory for KKT matrix */
		static pfloat Kpr[NNZK];
		static idxint Kir[NNZK];
		static idxint Kjc[nK+1];

	    // Populate Sign vector
		/* Set signs for regularization of (1,1) block */
	    for (idxint ks=0; ks < n; ks++ ){Sign[ks] = +1;} // (1,1) block
	    for (idxint ks=n; ks < n+p; ks++){Sign[ks] = -1;} // (2,2) block
	    for (idxint ks=n+p; ks < n+p+m; ks++) {Sign[ks] = -1;} // (3,3) block

	    /* count the number of non-zero entries in K */
	    idxint k = 0;

	    /* (1,1) block: the first n columns are empty */
	#if STATICREG == 0
	    for (idxint j=0; j<n; j++) {Kjc[j] = 0;}
	#else
	    for (idxint j=0; j<n; j++) {
	        Kjc[j] = j;
	        Kir[j] = j;
	        Kpr[k++] = DELTASTAT;
	    }
	#endif

	    /* Fill upper triangular part of K with values */
	    /* (1,2) block: A' */
		idxint i = 0; /* counter for non-zero entries in A or G, respectively */
		idxint row_stop = 0; // gets reused later
		idxint row = 0; // gets reused later
		for(idxint j=0; j<p; j++ ){
	        /* A' */
	        // diff in cumsum col counter gives # of rows of data in curr col
			row = At->jc[j];
			row_stop = At->jc[j+1];
			if( row <= row_stop ){
				Kjc[n+j] = k;
				while( row++ < row_stop ){
					Kir[k] = At->ir[i];
					Kpr[k] = At->pr[i];
					AttoK[i++] = k++;
				}
			}
	    // Same reg as earlier but now for (2,2) block
	#if STATICREG == 1
	        Kir[k] = n+j;
	        Kpr[k++] = -DELTASTAT;
	#endif
	    }
		/* (1,3) and (3,3) block: [G'; 0; -Vinit]
	     * where
	     *
	     *   Vinit = blkdiag(I, blkdiag(I,1,-1), ...,  blkdiag(I,1,-1));
	     *                        ^ #number of second-order cones ^
	     *
	     * Note that we have to prepare the (3,3) block accordingly
	     * (put zeros for init but store indices that are used in KKT_update
	     * of cone module)
	     */

		/* LP cone */
		i = 0; /* counter for non-zero entries in A or G, respectively */
		for(idxint j=0; j < setupspace.C->lpc->p; j++ ){
	        /* copy in G' */
	        // diff in cumsum col counter gives # of rows of data in curr col
			row = Gt->jc[j];
			row_stop = Gt->jc[j+1];
			if( row <= row_stop ){
				Kjc[n+p+j] = k;
				while( row++ < row_stop ){
					Kir[k] = Gt->ir[i];
					Kpr[k] = Gt->pr[i];
					GttoK[i++] = k++;
				}
			}
	        /* -I for LP-cone */
			setupspace.C->lpc->kkt_idx[j] = k;
			Kir[k] = n+p+j;
			Kpr[k] = -1.0;
	        k++;
		}

	    /* Second-order cones - copy in G' and set up the scaling matrix
	     * which has a dense structure (only upper half part is shown):
	     *
	     *                     [ a        b*q'     ]
	     *  - W^2 = -V = eta^2 [ b*q  I + c*(q*q') ]
	     *
	     * where    I: identity of size conesize
	     *          q: vector of size consize - 1
	     *      a,b,c: scalars
	     *
	     * NOTE: only the upper triangular part (with the diagonal elements)
	     *       is copied in here.
	     */
		idxint cone_strt = setupspace.C->lpc->p;
		conesize = 0; // WARNING: also used above in different context
	    for(idxint l=0; l < setupspace.C->nsoc; l++ ){

	        /* size of the cone */
	        conesize = setupspace.C->soc[l].p;

	        /* go column-wise about it */
			for(idxint j=0; j < conesize; j++ ){
	            // diff in cumsum col counter gives # of rows of data in curr col
	            row = Gt->jc[cone_strt+j];
	            row_stop = Gt->jc[cone_strt+j+1];
	            if( row <= row_stop ){
	                Kjc[n+p+cone_strt+j] = k;
	                while( row++ < row_stop ){
	                    Kir[k] = Gt->ir[i];
	                    Kpr[k] = Gt->pr[i];
	                    GttoK[i++] = k++;
	                }
	            }

	            /* first elements - record where this column starts */
	            Kir[k] = n+p+cone_strt;
	            Kpr[k] = -1.0;
	            setupspace.C->soc[l].colstart[j] = k;
	            k++; // important to keep counting for next Kjc

	            /* the rest of the column */
	            for (idxint r=1; r<=j; r++) {
	                Kir[k] = n+p+cone_strt+r;
	                Kpr[k] = -1.0;
	                k++; // important to keep counting for next Kjc
	            }
	        }

	        /* prepare index for next cone */
			cone_strt += setupspace.C->soc[l].p;
		}

	    PRINTTEXT("CREATEKKT: Written %d KKT entries\n", (int)k);
	    PRINTTEXT("CREATEKKT: Size of KKT matrix: %d\n", (int)nK);

		// KU = ecoscreateSparseMatrix(nK, nK, NNZK, Kjc, Kir, Kpr);
  	    spmat tmpKU;
	    tmpKU.m = nK;
	    tmpKU.n = nK;
	    tmpKU.nnz = NNZK;
	    tmpKU.jc = Kjc;
        tmpKU.ir = Kir;
        tmpKU.pr = Kpr;
	    if (tmpKU.jc)
	    	tmpKU.jc[nK] = NNZK;
        KU = &tmpKU;
        // ++End++ of KKTUpper setup

	    /* Save a mapping from data in A and G to the KKT matrix */
	    idxint AtoK[setupspace.A->nnz];
	    setupspace.AtoK = AtoK;
	    for(int i=0; i<setupspace.A->nnz; i++){ setupspace.AtoK[i] = AttoK[AtoAt[i]]; }
	    idxint GtoK[setupspace.G->nnz];
	    setupspace.GtoK = GtoK;
	    for(int i=0; i<setupspace.G->nnz; i++){ setupspace.GtoK[i] = GttoK[GtoGt[i]]; }

		// Set up KKT system related data
	    // (L comes later after symbolic factorization)
	    // STATICREG is used in preproc and kkt
	    // Currently set to 1

	    // PRINTTEXT("KKT matrix is: %d\n", (int)KKT);
	    PRINTTEXT("Dimension of KKT matrix: %d\n", (int)nK);
	    PRINTTEXT("Non-zeros in KKT matrix: %d\n", (int)KU->nnz);

	    // Allocate memory for KKT system
	    kkt tmpKKT;
	    static pfloat kkt_D[nK];
    	static idxint kkt_Parent[nK];
		static idxint kkt_Pinv[nK];
		static pfloat kkt_work1[nK];
		static pfloat kkt_work2[nK];
	    static pfloat kkt_work3[nK];
	    static pfloat kkt_work4[nK];
	    static pfloat kkt_work5[nK];
	    static pfloat kkt_work6[nK];
		static idxint kkt_Flag[nK];
		static idxint kkt_Pattern[nK];
		static idxint kkt_Lnz[nK];
		static pfloat kkt_RHS1[nK];
		static pfloat kkt_RHS2[nK];
		static pfloat kkt_dx1[n];
		static pfloat kkt_dx2[n];
		static pfloat kkt_dy1[p];
		static pfloat kkt_dy2[p];
		static pfloat kkt_dz1[m];
		static pfloat kkt_dz2[m];
	    static idxint kkt_Sign[nK];

  	    spmat tmpkkt_PKPt;
  	    static idxint PKPT_jc[nK+1];
  	    static idxint PKPT_ir[NNZK];
  	    static pfloat PKPT_pr[NNZK];
	    tmpkkt_PKPt.m = nK;
	    tmpkkt_PKPt.n = nK;
	    tmpkkt_PKPt.nnz = NNZK;	
	    tmpkkt_PKPt.jc = PKPT_jc;
        tmpkkt_PKPt.ir = PKPT_ir;
        tmpkkt_PKPt.pr = PKPT_pr;
	    if (tmpkkt_PKPt.jc)
	    	tmpkkt_PKPt.jc[nK] = NNZK;
        spmat* kkt_PKPt = &tmpkkt_PKPt;

	    // spmat* kkt_PKPt = newSparseMatrix(nK, nK, KU->nnz);
		static idxint kkt_PK[NNZK];

	    tmpKKT.D        = kkt_D;
    	tmpKKT.Parent   = kkt_Parent;
		tmpKKT.Pinv     = kkt_Pinv;
		tmpKKT.work1    = kkt_work1;
		tmpKKT.work2    = kkt_work2;
	    tmpKKT.work3    = kkt_work3;
	    tmpKKT.work4    = kkt_work4;
	    tmpKKT.work5    = kkt_work5;
	    tmpKKT.work6    = kkt_work6;
		tmpKKT.Flag     = kkt_Flag;
		tmpKKT.Pattern  = kkt_Pattern;
		tmpKKT.Lnz      = kkt_Lnz;
		tmpKKT.RHS1     = kkt_RHS1;
		tmpKKT.RHS2     = kkt_RHS2;
		tmpKKT.dx1      = kkt_dx1;
		tmpKKT.dx2      = kkt_dx2;
		tmpKKT.dy1      = kkt_dy1;
		tmpKKT.dy2      = kkt_dy2;
		tmpKKT.dz1      = kkt_dz1;
		tmpKKT.dz2      = kkt_dz2;
	    tmpKKT.Sign     = kkt_Sign;
	    tmpKKT.PKPt     = kkt_PKPt;
		tmpKKT.PK       = kkt_PK;
		setupspace.KKT = &tmpKKT;
		// printf("Value for kkt_Flag[%d]: is %d\n", 1, kkt_Flag[1]);
		// printf("Value for setupspace->KKT->kkt_Flag[%d]: is %p\n", 1, setupspace->KKT->Flag);
	    // printf("Address of Sign1 is %p\n", setupspace->KKT->Sign);
	    // printf("Address of Sign1[0] is %p\n", &setupspace->KKT->Sign[0]);

	    /* calculate ordering of KKT matrix using AMD */
		idxint P[nK];

		double Control [AMD_CONTROL], Info [AMD_INFO];
		AMD_defaults(Control);
		idxint amd_result = AMD_order(nK, KU->jc, KU->ir, P, Control, Info);

		if( amd_result == AMD_OK ){
			// do nothing
		} else {
			PRINTTEXT("Problem in AMD ordering, exiting.\n");
	        AMD_info(Info);
	        return NULL;
		}

		/* calculate inverse permutation and permutation mapping of KKT matrix */
		// pinv(nK, P, setupspace.KKT->Pinv);
		{
			idxint i;
			for( i=0; i<nK; i++ ){ setupspace.KKT->Pinv[P[i]] = i; }
		}
		idxint *Pinv = setupspace.KKT->Pinv;

		permuteSparseSymmetricMatrix(KU, setupspace.KKT->Pinv, setupspace.KKT->PKPt, setupspace.KKT->PK);

		/* permute sign vector */
	    for(int i=0; i<nK; i++ ){ setupspace.KKT->Sign[Pinv[i]] = Sign[i]; }

		/* symbolic factorization */
		idxint Ljc[nK+1];
	    // Allocate memory for cholesky factor L
 		// pr (x) is each data point
 		// jc (p) is cumsum elements per col (size n+1 b/c starts at 0)
		// ir (i) is row of each data point
		LDL_symbolic2(
			setupspace.KKT->PKPt->n,    /* A and L are n-by-n, where n >= 0 */
			setupspace.KKT->PKPt->jc,   /* input of size n+1, not modified */
			setupspace.KKT->PKPt->ir,	 /* input of size nz=Ap[n], not modified */
			Ljc,					 /* output of size n+1, not defined on input */
			setupspace.KKT->Parent,	 /* output of size n, not defined on input */
			setupspace.KKT->Lnz,		 /* output of size n, not defined on input */
			setupspace.KKT->Flag		 /* setupspace of size n, not defn. on input or output */
		);


	    // Created Cholesky factor (L) of K in KKT struct
		idxint lnz = Ljc[nK]; // number of elements up to and inc. nK-th col
		idxint *Lir = (idxint *)MALLOC(lnz*sizeof(idxint));
		pfloat *Lpr = (pfloat *)MALLOC(lnz*sizeof(pfloat));
		setupspace.KKT->L = ecoscreateSparseMatrix(nK, nK, lnz, Ljc, Lir, Lpr);


		/* permute KKT matrix - we work on this one from now on */
		permuteSparseSymmetricMatrix(KU, setupspace.KKT->Pinv, setupspace.KKT->PKPt, NULL);

		/* get memory for residuals */
		pfloat setupspace_rx[n];
		pfloat setupspace_ry[n];
		pfloat setupspace_rz[n];
		setupspace.rx = setupspace_rx;
		setupspace.ry = setupspace_ry;
		setupspace.rz = setupspace_rz;

	    /* clean up */
	    setupspace.KKT->P = P;
 	 // ++End++ of ECOS_setup (should be replaced by codegen)

		initialized = 1;

		// workspace = &setupspace; // doesn't allow setting ptr to NULL
		workspace = setupspace_ptr;
	}
	printf("workspace->KKT->PKPt->n is %d\n", workspace->KKT->PKPt->n);
	// Defn's from ECOS_solve
	idxint exitcode = ECOS_FATAL, interrupted = 0;
	pfloat pres_prev = (pfloat)ECOS_NAN;

	// ++Start++ of init
	// Workspace Initialization
	idxint* Pinv = workspace->KKT->Pinv;

    /* Initialize KKT matrix */
    kkt_init(workspace->KKT->PKPt, workspace->KKT->PK, workspace->C);

    /* initialize RHS1 */
	idxint k = 0, j = 0;
	for(int i=0; i<workspace->n; i++ ){ workspace->KKT->RHS1[workspace->KKT->Pinv[k++]] = 0; }
	for(int i=0; i<workspace->p; i++ ){ workspace->KKT->RHS1[workspace->KKT->Pinv[k++]] = workspace->b[i]; }
	for(int i=0; i<workspace->C->lpc->p; i++ ){ workspace->KKT->RHS1[workspace->KKT->Pinv[k++]] = workspace->h[i]; j++; }
	for(int l=0; l<workspace->C->nsoc; l++ ){
		for(int i=0; i < workspace->C->soc[l].p; i++ ){ workspace->KKT->RHS1[workspace->KKT->Pinv[k++]] = workspace->h[j++]; }
	}

	/* initialize RHS2 */
	for(int i=0; i<workspace->n; i++ ){ workspace->KKT->RHS2[workspace->KKT->Pinv[i]] = -workspace->c[i]; }
	for(int i=workspace->n; i<workspace->KKT->PKPt->n; i++ ){ workspace->KKT->RHS2[workspace->KKT->Pinv[i]] = 0; }

	/* get scalings of problem data */
	pfloat rx = norm2(workspace->c, workspace->n); workspace->resx0 = MAX(1, rx);
	pfloat ry = norm2(workspace->b, workspace->p); workspace->resy0 = MAX(1, ry);
	pfloat rz = norm2(workspace->h, workspace->m); workspace->resz0 = MAX(1, rz);

	/* Factor KKT matrix - this is needed in all 3 linear system solves */
    idxint KKT_FACTOR_RETURN_CODE = kkt_factor(workspace->KKT, workspace->stgs->eps, workspace->stgs->delta);

    /* check if factorization was successful, exit otherwise */
	if(  KKT_FACTOR_RETURN_CODE != KKT_OK ){
        PRINTTEXT("\nProblem in factoring KKT system, aborting.");
        return ECOS_FATAL;
    }

    // PRIMAL INITIALIZATION
	workspace->info->nitref1 = kkt_solve(workspace->KKT, workspace->A, workspace->G, workspace->KKT->RHS1, workspace->KKT->dx1, workspace->KKT->dy1, workspace->KKT->dz1, workspace->n, workspace->p, workspace->m, workspace->C, 1, workspace->stgs->nitref);
	/* Copy out initial value of x */
	for(int i=0; i<workspace->n; i++ ){ workspace->x[i] = workspace->KKT->dx1[i]; }
	/* Copy out -r into temporary variable */
	for(int i=0; i<workspace->m; i++ ){ workspace->KKT->work1[i] = -workspace->KKT->dz1[i]; }
	/* Bring variable to cone */
	bring2cone(workspace->C, workspace->KKT->work1, workspace->s );

	// DUAL INITIALIZATION
	workspace->info->nitref2 = kkt_solve(workspace->KKT, workspace->A, workspace->G, workspace->KKT->RHS2, workspace->KKT->dx2, workspace->KKT->dy2, workspace->KKT->dz2, workspace->n, workspace->p, workspace->m, workspace->C, 1, workspace->stgs->nitref);
    /* Copy out initial value of y */
	for(int i=0; i<workspace->p; i++ ){ workspace->y[i] = workspace->KKT->dy2[i]; }
	/* Bring variable to cone */
	bring2cone(workspace->C, workspace->KKT->dz2, workspace->z );

	/* Prepare RHS1 - before this line RHS1 = [0; b; h], after it holds [-c; b; h] */
	for(int i=0; i<workspace->n; i++){ workspace->KKT->RHS1[Pinv[i]] = -workspace->c[i]; }

	/*
	 * other variables
	 */
	workspace->kap = 1.0;
	workspace->tau = 1.0;

	workspace->info->step = 0;
	workspace->info->step_aff = 0;
	workspace->info->dinf = 0;
	workspace->info->pinf = 0;
	// ++End++ of init

	for (int i=0; i<workspace->n; i++) {printf("Init value for x[%d] is %e\n", i, (pfloat)workspace->x[i]);}
	for (int i=0; i<workspace->p; i++) {printf("Init value for y[%d] is %e\n", i, (pfloat)workspace->y[i]);}
	for (int i=0; i<workspace->m; i++) {printf("Init value for s[%d] is %e\n", i, (pfloat)workspace->s[i]);}
	for (int i=0; i<workspace->m; i++) {printf("Init value for z[%d] is %e\n", i, (pfloat)workspace->z[i]);}

	// =====================================
	// ===== Core IPM Resides Here =========
	// =====================================
	// for (int k=0; k<maxit; k++) {
	for (int k=0; k<10; k++) {
		computeResiduals(workspace);
		printf("Residual for nx is %e\n", (pfloat)workspace->nx);
		updateStatistics(workspace);

		// STEP 1.5 CHECK TERMINATION CRITERIA HERE
		// check_termination(x_hat, y_hat, z_hat, t_hat, k_hat, h, b, c)

		// ++Start++ of updateScalings
		cone* C = workspace->C;
		pfloat* s = workspace->s;
		pfloat* z = workspace->z;
		pfloat* lambda = workspace->lambda;

		idxint i, l, k, p; /*, pm1; */
		pfloat sres, zres, snorm, znorm, gamma, one_over_2gamma;
		pfloat* sk;
		pfloat* zk;
	    pfloat a, c, d, w, temp, divisor; /*, b; */
	    pfloat u0, u0_square, u1, v1, d1, c2byu02_d, c2byu02;

		// Update w for LP cones
		for (int i=0; i < C->lpc->p; i++ ){
			C->lpc->v[i] = SAFEDIV_POS(s[i], z[i]);
			C->lpc->w[i] = sqrt(C->lpc->v[i]); // w = sqrt(s/z)
		}

		// Update w for SO cones
		int soc_start_idx = C->lpc->p;
		for(int l=0; l < C->nsoc; l++ ){

			/* indices and variables */
			sk = s+soc_start_idx; // just s is a pointer to s[0]
			zk = z+soc_start_idx; // just z is a pointer to z[0]
			p = C->soc[l].p; /* pm1 = p-1; */

			// Check residuals
			sres = socres(sk, p);
			zres = socres(zk, p);

			/* normalize variables */
			snorm = sqrt(sres);
			znorm = sqrt(zres);
			for(int i=0; i<p; i++ ){ C->soc[l].skbar[i] = SAFEDIV_POS(sk[i],snorm); }
			for(int i=0; i<p; i++ ){ C->soc[l].zkbar[i] = SAFEDIV_POS(zk[i],znorm); }
			C->soc[l].eta_square = SAFEDIV_POS(snorm,znorm);
			C->soc[l].eta = sqrt(C->soc[l].eta_square);

			/* Normalized Nesterov-Todd scaling point */
			// Calculate gamma and 1/2gamma
			gamma = 1.0;
			for(int i=0; i<p; i++){ gamma += C->soc[l].skbar[i]*C->soc[l].zkbar[i]; }
			gamma = sqrt(0.5*gamma);
			one_over_2gamma = SAFEDIV_POS(0.5,gamma);
			// Define a=wkbar[0]
			a = one_over_2gamma*(C->soc[l].skbar[0] + C->soc[l].zkbar[0]);
			// Define q=wkbar[1:] and w=wkbar[1:].T@wkbar[1:]
			w = 0;
			for(int i=1; i<p; i++ ){
				C->soc[l].q[i-1] = one_over_2gamma*(C->soc[l].skbar[i] - C->soc[l].zkbar[i]);
				w += C->soc[l].q[i-1]*C->soc[l].q[i-1];
			}
	        C->soc[l].w = w;
	        C->soc[l].a = a;

			/* pre-compute variables needed for KKT matrix (kkt_update uses those) */
	        temp = 1.0 + a; // = wkbar[0] + 1.0
	        /* b = SAFEDIV_POS(1.0,temp); */
	        c = 1.0 + a + SAFEDIV_POS(w,temp); // wkbar[0] + 1.0 + wkbar[1:].T@wkbar[1:]/(wkbar[0] + 1.0)
	        divisor = temp*temp; // (wkbar[0] + 1.0)^2
	        // Below is 1 + 2/(wkbar[0] + 1.0) + wkbar[1:].T@wkbar[1:]/(wkbar[0] + 1.0)^2
	        d = 1 + SAFEDIV_POS(2,temp) + SAFEDIV_POS(w,divisor);

	        C->soc[l].c = c;
	        C->soc[l].d = d;

			/* increase offset for next cone */
			soc_start_idx += C->soc[l].p;
		}

		// Compute lambda = W*z
		idxint j, cone_start;
		pfloat zeta, factor;

		/* LP cone */
		// C->lpc->p is dimension of LP cone
		for(int i=0; i < C->lpc->p; i++ )
			lambda[i] = C->lpc->w[i] * z[i];

		/* Second-order cone */
		cone_start = C->lpc->p;
		// C->nsoc is number of SO cones
		for(int l=0; l < C->nsoc; l++ ){

			// zeta = q'*z1 (remember q=wkbar[1:])
			zeta = 0; // = wkbar[1:].T@zk[1:]
			// soc[l].p is dimension of l-th SO cone
			for(int i=1; i < C->soc[l].p; i++ )
				zeta += C->soc[l].q[i-1] * z[cone_start + i];

			// factor = z0 + zeta / (1+a)
			//		  = z0 + wkbar[1:].T@zk/(1+wkbar[0])
			factor = z[cone_start] + SAFEDIV_POS(zeta,(1+C->soc[l].a));

			// second pass (on k): write out result
			// lamk0 = eta*(wkbar[0]*zk0 + wkbar[1:].T@zk)
			//       = eta*(wkbar.T@zk)
			lambda[cone_start] = C->soc[l].eta*(C->soc[l].a*z[cone_start] + zeta); /* lambda[0] */
			for(int i=1; i < C->soc[l].p; i++ ){
				j = cone_start+i;
				// lambdak1 = eta*(zkbar[1:].T@wkbar[1:])
				lambda[j] = C->soc[l].eta*(z[j] + factor*C->soc[l].q[i-1]);
			}

			cone_start += C->soc[l].p;
		}
		// ++End++ of updateScalings

		// STEP 2 - SOLVE FOR THE AFFINE DIRECTION
		/* Update KKT matrix with scalings */
		kkt_update(workspace->KKT->PKPt, workspace->KKT->PK, workspace->C);

        /* factor KKT matrix */
        KKT_FACTOR_RETURN_CODE = kkt_factor(workspace->KKT, workspace->stgs->eps, workspace->stgs->delta);

	    /* check if factorization was successful, exit otherwise */
		if(  KKT_FACTOR_RETURN_CODE != KKT_OK ){
	       if( workspace->stgs->verbose ) PRINTTEXT("\nProblem in factoring KKT system, aborting.");
	        return ECOS_FATAL;
	    }

		/* Solve for RHS1, which is used later also in combined direction */
		workspace->info->nitref1 = kkt_solve(workspace->KKT, workspace->A, workspace->G, workspace->KKT->RHS1, workspace->KKT->dx1, workspace->KKT->dy1, workspace->KKT->dz1, workspace->n, workspace->p, workspace->m, workspace->C, 0, workspace->stgs->nitref);

		/* AFFINE SEARCH DIRECTION (predictor, need dsaff and dzaff only) */
		RHS_affine(workspace);
		workspace->info->nitref2 = kkt_solve(workspace->KKT, workspace->A, workspace->G, workspace->KKT->RHS2, workspace->KKT->dx2, workspace->KKT->dy2, workspace->KKT->dz2, workspace->n, workspace->p, workspace->m, workspace->C, 0, workspace->stgs->nitref);

		/* dtau_denom = kap/tau - (c'*x1 + by1 + h'*z1); */
		pfloat dtau_denom = workspace->kap/workspace->tau - eddot(workspace->n, workspace->c, workspace->KKT->dx1) - eddot(workspace->p, workspace->b, workspace->KKT->dy1) - eddot(workspace->m, workspace->h, workspace->KKT->dz1);

        /* dtauaff = (dt + c'*x2 + by2 + h'*z2) / dtau_denom; */
		pfloat dtauaff = (workspace->rt - workspace->kap + eddot(workspace->n, workspace->c, workspace->KKT->dx2) + eddot(workspace->p, workspace->b, workspace->KKT->dy2) + eddot(workspace->m, workspace->h, workspace->KKT->dz2)) / dtau_denom;

		/* dzaff = dz2 + dtau_aff*dz1 */
        /* let dz2   = dzaff  we use this in the linesearch for unsymmetric cones*/
        /* and w_times_dzaff = Wdz_aff*/
        /* and dz2 = dz2+dtau_aff*dz1 will store the unscaled dz*/
		for( i=0; i<workspace->m; i++ ){ workspace->KKT->dz2[i] = workspace->KKT->dz2[i] + dtauaff*workspace->KKT->dz1[i]; }
		scale(workspace->KKT->dz2, workspace->C, workspace->W_times_dzaff);

		/* W\dsaff = -W*dzaff -lambda; */
		for( i=0; i<workspace->m; i++ ){ workspace->dsaff_by_W[i] = -workspace->W_times_dzaff[i] - workspace->lambda[i]; }

		/* dkapaff = -(bkap + kap*dtauaff)/tau; bkap = kap*tau*/
		pfloat dkapaff = -workspace->kap - workspace->kap/workspace->tau*dtauaff;

		/* Line search on W\dsaff and W*dzaff */
		workspace->info->step_aff = lineSearch(workspace->lambda, workspace->dsaff_by_W, workspace->W_times_dzaff, workspace->tau, dtauaff, workspace->kap, dkapaff, workspace->C, workspace->KKT);

		/* Centering parameter */
        pfloat sigma = 1.0 - workspace->info->step_aff;
        sigma = sigma*sigma*sigma;
        if( sigma > SIGMAMAX ) sigma = SIGMAMAX;
        if( sigma < SIGMAMIN ) sigma = SIGMAMIN;
        workspace->info->sigma = sigma;

		// STEP 4 - SOLVE FOR THE COMBINED DIRECTION
		// Very similar to STEP 3. Reuse some of the material from there.
		// Have to manually calculate ds term by splitting up vectors into cones
		// ds = lam O lam + (W^(-T)@del_sa) O (W@del_za) + sigma*mu_hat*e
		RHS_combined(workspace);
		workspace->info->nitref3 = kkt_solve(workspace->KKT, workspace->A, workspace->G, workspace->KKT->RHS2, workspace->KKT->dx2, workspace->KKT->dy2, workspace->KKT->dz2, workspace->n, workspace->p, workspace->m, workspace->C, 0, workspace->stgs->nitref);

  		/* bkap = kap*tau + dkapaff*dtauaff - sigma*info.mu; */
		pfloat bkap = workspace->kap*workspace->tau + dkapaff*dtauaff - sigma*workspace->info->mu;

		/* dtau = ((1-sigma)*rt - bkap/tau + c'*x2 + by2 + h'*z2) / dtau_denom; */
		pfloat dtau = ((1-sigma)*workspace->rt - bkap/workspace->tau + eddot(workspace->n, workspace->c, workspace->KKT->dx2) + eddot(workspace->p, workspace->b, workspace->KKT->dy2) + eddot(workspace->m, workspace->h, workspace->KKT->dz2)) / dtau_denom;

		/* dx = x2 + dtau*x1;     dy = y2 + dtau*y1;       dz = z2 + dtau*z1; */
		for(int i=0; i < workspace->n; i++ ){ workspace->KKT->dx2[i] += dtau*workspace->KKT->dx1[i]; }
		for(int i=0; i < workspace->p; i++ ){ workspace->KKT->dy2[i] += dtau*workspace->KKT->dy1[i]; }
		for(int i=0; i < workspace->m; i++ ){ workspace->KKT->dz2[i] += dtau*workspace->KKT->dz1[i]; }

		/*  ds_by_W = -(lambda \ bs + conelp_timesW(scaling,dz,dims)); */
		/* note that at this point workspace->dsaff_by_W holds already (lambda \ ds) */
		scale(workspace->KKT->dz2, workspace->C, workspace->W_times_dzaff);
		for(int i=0; i < workspace->m; i++ ){ workspace->dsaff_by_W[i] = -(workspace->dsaff_by_W[i] + workspace->W_times_dzaff[i]); }

		/* dkap = -(bkap + kap*dtau)/tau; */
		pfloat dkap = -(bkap + workspace->kap*dtau)/workspace->tau;

		// STEP 5 - SOLVE FOR THE STEP SIZE AND UPDATE ITERATES AND SCALING MATRICES
		/* Line search on combined direction */
		workspace->info->step = lineSearch(workspace->lambda, workspace->dsaff_by_W, workspace->W_times_dzaff, workspace->tau, dtau, workspace->kap, dkap, workspace->C, workspace->KKT) * workspace->stgs->gamma;

        /* Bring ds to the final unscaled form */
	    /* ds = W*ds_by_W */
		scale(workspace->dsaff_by_W, workspace->C, workspace->dsaff);


		/* Update variables */
		for(int i=0; i < workspace->n; i++ ){ workspace->x[i] += workspace->info->step * workspace->KKT->dx2[i]; }
		for(int i=0; i < workspace->p; i++ ){ workspace->y[i] += workspace->info->step * workspace->KKT->dy2[i]; }
		for(int i=0; i < workspace->m; i++ ){ workspace->z[i] += workspace->info->step * workspace->KKT->dz2[i]; }
		for(int i=0; i < workspace->m; i++ ){ workspace->s[i] += workspace->info->step * workspace->dsaff[i]; }
		workspace->kap += workspace->info->step * dkap;
		workspace->tau += workspace->info->step * dtau;

	}
	

	// Unscale variables by tau
    for(int i=0; i < workspace->n; i++ ){ workspace->x[i] /= workspace->tau; }
    for(int i=0; i < workspace->p; i++ ){ workspace->y[i] /= workspace->tau; }
	for(int i=0; i < workspace->m; i++ ){ workspace->z[i] /= workspace->tau; }
	for(int i=0; i < workspace->m; i++ ){ workspace->s[i] /= workspace->tau; }
	
	for (int i=0; i<workspace->n; i++) {printf("Final value for x[%d] is %e\n", i, (pfloat)workspace->x[i]);}
	for (int i=0; i<workspace->p; i++) {printf("Final value for y[%d] is %e\n", i, (pfloat)workspace->y[i]);}
	for (int i=0; i<workspace->m; i++) {printf("Final value for s[%d] is %e\n", i, (pfloat)workspace->s[i]);}
	for (int i=0; i<workspace->m; i++) {printf("Final value for z[%d] is %e\n", i, (pfloat)workspace->z[i]);}

    return 0;
}
