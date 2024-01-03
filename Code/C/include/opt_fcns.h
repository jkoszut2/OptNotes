// Optimization Functions


#define DBL_EPS 2.2204460492503131e-16 // double

int process_cone_info(int *cone_vars, int *cone_start_idxs, int *l_idx,
                      int *e_socs, int *m) {

    
    // Find indices corresponding to starts of cones
    for (int i=1; i<num_cones; i++) {
        cone_start_idxs[i] = cone_start_idxs[i-1] + cone_vars[i-1];
    }

    // Detect start of SOCs by finding where order is first >1
    int soc_idx = 0;
    for (int i=0; i<num_cones; i++) {
        if (soc_idx==0 & cone_orders[i]>1) {
            soc_idx = i;
        }
    }

    // Calculate variable index corresponding to start of SOCs
    for (int i=0; i<soc_idx; i++) {
        *l_idx += cone_vars[i];
    }
    
    // Set up \textbb{e} vector
    for (int i=0; i<num_cones; i++) {
        int cone_start_idx = cone_start_idxs[i];
        e_socs[cone_start_idx] = 1;
        if((cone_orders[i] == 1) & (cone_vars[i]>1)) {
            for (int j=1; j<cone_vars[i]; j++) {
                // Order-1 cone has all 1's for ek
                e_socs[cone_start_idx+j] = 1;
            }
        }
    }

    // Calculate total cone order
    for (int i=0; i<num_cone_vars; i++) {
        *m += e_socs[i];
    }

    return 0;
}