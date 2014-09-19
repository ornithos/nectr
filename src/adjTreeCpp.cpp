#include <R.h>
#include <Rinternals.h>

extern "C" {

class stack
{
// Stack explicitly built to house **non-negative** integers.
  private:
    int *arr;
    int top;
	int max_stack;

  public:
	 stack(const int stacksize)         //Constructor
	 {
	    top=-1;      //Sets the Top Location to -1 indicating an empty stack
		arr = new int[stacksize];
		max_stack = stacksize;
	 }

	 void push(int a) {
		if(top + 1 < max_stack) {
			arr[++top] = a;
		 }
	 }

	 int pop() {
		if(top >= 0) {
			int data=arr[top];
			arr[top] = 0;
			top--;
			return data;
		}
		else {
			return -1;
		}
	 }

	 ~stack() //Destructor
	 {
		 delete[] arr;
		 arr = 0;
	 }
};

SEXP doTurnTraversalCpp(SEXP m_nb, SEXP v_type, SEXP v_ord, SEXP v_clsset)
{
	//VARIABLE DECLARATIONS AND POINTERS ##################################################

    int start, clsfound_i;
	int *pm_nb, *pv_type, *pv_ord, *pv_clsset, *pout, *pv_found, *pv_setlkp, *pv_temp;
	int n = nrows(m_nb);
	int p = ncols(m_nb);
	pm_nb = INTEGER(m_nb);
	pv_type = INTEGER(v_type);
	pv_ord = INTEGER(v_ord);
	pv_clsset = INTEGER(v_clsset);
	
	//Create output and temporary vector
	SEXP out = PROTECT(allocVector(INTSXP, n));
	SEXP v_found = PROTECT(allocVector(INTSXP, n));
	SEXP v_setlkp = PROTECT(allocVector(INTSXP, n));
	SEXP v_temp = PROTECT(allocVector(INTSXP, n));
	
	//pointer to array in temporary vectors
	pout = INTEGER(out);
	pv_found = INTEGER(v_found);
	pv_setlkp = INTEGER(v_setlkp);
	pv_temp = INTEGER(v_temp);
	
	//initialise created vectors
	memset(pout, 0, n * sizeof(int));
	memset(pv_found, 0, n * sizeof(int));
	memset(pv_setlkp, 0, n * sizeof(int));
	
	//Create node and col stacks
	stack node_stack(n + p);
	stack col_stack(n + p);
	
	//  Miscellaneous variables for traversal
	int j, cls, x, lbl;
	int r = 0;
	int new_r = 0;	
	bool bInternal, bVisited;
	
	int unique_bleed = 0;
	bool bUnique = true;
	// ####################################################################################
	
	
    //outer loop - will fire from each node in order as the starting node.
    for(int i = 0; i < n; i++)
    {	
        //Set # clusters 'found' during traversal to 0 (as beginning new run)
		// and start to the beginning node.
		clsfound_i = 0;
        start = pv_ord[i];
		
		cls = i+1;
		r = start;
		j = 0;
		bVisited = false;		//Visited datapoint within same traversal - for travelling back down the stack
		
		//Run the traversal from this node ('start') - providing starting node not already labelled.
		// The 'within cluster' loop
        if (pout[start] == 0) {

			while(true) {
				// Label of node that we have traversed to.
				x = pout[r];

				if(!bVisited) {
					lbl = pv_clsset[r];
					bInternal = pv_type[r] == 1;
			
					//If node already labelled by a previous traversal, we want this branch of the
					// recursion to halt. (If part of a different cluster (and internal), we add this
					// to the 'found' array, and that cluster will be subsumed by the current)
					if (x != 0) {
						if (x != cls && bInternal) {
							if (pv_setlkp[x-1] == pv_setlkp[cls-1]) {
								pv_found[clsfound_i++] = x;
							}
						}
					
						// Since found cluster label, continue up parent stack.
						do {						//BACK UP PARENT STACK
							r = node_stack.pop();
							if(r == -1) goto end_traversal_loop;
							j = col_stack.pop();
						} while (j == p);
						bVisited = true;
						continue;
					}

					//ClusterSet label of row (lbl_new). If both the current Clusterset and the row ClusterSet exist (non-zero),
					// and they are not the same (ie demarcated as different clusters), stop the traversal and do not assign to current cls.
					int lbl_new = pv_clsset[r];
					if ((lbl > 0) && (lbl_new > 0) && (lbl != lbl_new)) {
						do {						//BACK UP PARENT STACK
							r = node_stack.pop();
							if(r == -1) goto end_traversal_loop;
							j = col_stack.pop();
						} while (j == p);
						bVisited = true;
						continue;
					}

					//Label node as part of cluster
					//if(cls == 1) Rprintf("%d labelled in cls 1\n", r+1);
					pout[r] = cls;

					// If ClusterSet label for row exists, and there is not a clash (tested above), assign current cluster to this ClusterSet.
					if(lbl_new > 0) pv_setlkp[cls-1] =  lbl_new;
				}
				else { bInternal = true; }

				//If more neighbours left in current row, loop until find next non-null, add current row and col to stack
				// and traverse with this new row.
				if (bInternal) {
						
					//Find next non-null neighbour
					while(j < p) {
						new_r = pm_nb[n*j++ + r] - 1;
						if(new_r != -1) break;
					}
						
					//If found a non-null neighbour, add current pos to stack and continue outer loop with new row.
					if(new_r != -1) {
						node_stack.push(r);
						col_stack.push(j);
						r = new_r;
						j = 0;
						bVisited = false;
						continue;
					}
				}
				
				//At this point, either the last neighbour in row has been traversed, no non-null neighbours are left,
				// or point is Internal. For all of these cases, we need to traverse up the stack.

				do {						//BACK UP PARENT STACK
					r = node_stack.pop();
					if(r == -1) goto end_traversal_loop;
					j = col_stack.pop();
				} while (j == p);
				bVisited = true;

			}
		}
end_traversal_loop:
		
		//If we have bled into existing clusters as stored in v_found, relabel those as current.
        if (clsfound_i > 0) {
				
			//Extract unique clusters found from pv_found
			unique_bleed = 0;
			for(int k = 0; k < clsfound_i; k++) {
				bUnique = true;
				for(int m = 0; m < unique_bleed; m++) {
					if(pv_found[k] == pv_temp[m]) bUnique = false;
				}
				if(bUnique) pv_temp[unique_bleed++] = pv_found[k];
			}
			
            for(int k = 0; k < n; k++) {
                for(int m = 0; m < unique_bleed; m++) {
                    if(pout[k] == pv_temp[m]) {
                        pout[k] = i + 1;
                    }
                }
            }
        }
    }
	
	UNPROTECT(4);
	return out;
}

}