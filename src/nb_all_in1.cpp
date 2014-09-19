#include "R.h"

extern "C" {

////////////////////////////////////////////////////////
// MERGE SORT
////////////////////////////////////////////////////////

void merge(int* a, int* indx, int* aux, int lo, int mid, int hi)
{
	//copy original array to auxiliary
    for(int i = lo; i <= hi; i++) {
        aux[i] = indx[i];
    }

    //do merge
    int j = lo;
    int k = mid+1;
    for(int i = lo; i <= hi; i++) {
        if (j > mid) {
            indx[i] = aux[k++];
        } else if (k > hi) {
            indx[i] = aux[j++];
        } else if (a[aux[j]] > a[aux[k]]) {
            indx[i] = aux[k++];
        } else {
            indx[i] = aux[j++];
        }
    }

	
}

 
void orderix(int* a, int* indx, int* aux, int lo, int hi)
{
    if (hi <= lo) {return;}
    int mid = lo + (hi - lo)/2;
    orderix(&a[0], &indx[0], &aux[0], lo, mid);
    orderix(&a[0], &indx[0], &aux[0], mid+1, hi);
    merge(&a[0], &indx[0], &aux[0], lo, mid, hi);
}



void morder(int* a, int size, int* out)
{
	int* aux = new int[size];
	
	for(int i = 0; i < size; i++) {
		out[i] = i;
	}

	orderix(&a[0], &out[0], &aux[0], 0, size - 1);

	delete [] aux;
}

////////////////////////////////////////////////////////
// CALC NEIGHBOURS	
////////////////////////////////////////////////////////

void reorderMatrix(int* m, int* cpy, int* order, int n, int p)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < p; j++)
		{
			cpy[j*n + i] = m[j*n + order[i]];
		}
	}

	for(int i = 0; i < p*n; i++) {
		m[i] = cpy[i];
	}

	return;
}

void reorderVec(int* v, int* cpy, int* order, int n)
{
	for(int i = 0; i < n; i++)
	{
		cpy[i] = v[order[i]];
	}

	for(int i = 0; i <n; i++) {
		v[i] = cpy[i];
	}
	return;
}

void calculateNeighbours(int* iv_id, int* im_pos, int* im_nbs, double* iv_ss, int* iv_type, int* nP, int* pP, double* phiP)
{
	//Initialise input variables
	int* v_id = iv_id;
	int* m_pos = im_pos;
	int* m_nbs = im_nbs;
	int* v_type = iv_type;
	double* v_ss = iv_ss;
	int n = *nP;
	int p = *pP;
	double phi = *phiP;

	//Initialise temporary variables
	bool clsenb = false;
	int d = 0;
	int itmp = 0;
	double dtmp = 0;
	int* col = new int[n];
	int* new_order = new int[n];
	int* id_order = new int[n];
	int* v_id_copy = new int[n];
	int* m_pos_copy = new int[p*n];	
	int* m_nbs_tmp = new int[2*n];
	double* v_tmp = new double[n];

	//initialise new arrrays as required.
	for(int i = 0; i < n; i++) {
		m_nbs_tmp[2*i] = 0;
		m_nbs_tmp[2*i + 1] = 0;
		v_tmp[i] = 0;
	}

	//	***********************************************************
	//	--Outer Loop: for dimension d, get left and right neighbours:
	//  ***********************************************************
	//	  start with column 1 (zero indexed in C)
		for(d = 0; d < p; d++)
		{

	//		--get the order of this column (1. copying contents of column to feed into mergesort).
	//		-- (2. Create order, returned in 'new_order')
			for(int i = 0; i < n; i++) {
				col[i] = m_pos[d*n + i];
			}
			morder(col, n, new_order);


	//		--Perform re-ordering of 'm_nbs' and 'v_id' based on this order.
	//		--Need inverse of v_id for returning results - returned in 'id_order'
			reorderMatrix(m_pos, m_pos_copy, new_order, n, p);
			reorderVec(v_id, v_id_copy, new_order, n);
			morder(v_id, n, id_order);

	//		*********************************************************
	//		--Inner Loop: This is the actual data generation: the 'Get Neighbours' part of the script.
	//		*********************************************************
			for(int i = 0; i < n-1; i++) {

		//		--Insert above and below neighbours from the current order into m_nbs output.
				m_nbs_tmp[i] = v_id[i+1];
				m_nbs_tmp[n + i+1] = v_id[i];

		//		--Calculate metrics by assessing the entire row (1:p)
				for(int c = 0; c < p; c++) {

					//Get position of neighbour in dimension c. By looking at Right neighbour @ n+1, and left
					// neighbour @ n, we can re-use a single distance value, as the only difference will be the sign.
					// we're not concerned with the sign. itmp = distance. Define nb as close:= abs(dist) < 1.
					
					itmp = m_pos[i+1 + c*n] - m_pos[i + c*n];
					((itmp > 1) || (itmp < -1)) ? clsenb = false : clsenb = true;

			//		.....For 'Right' Neighbour................
					if (!clsenb) {
						m_nbs_tmp[i] = 0;
					}
			//		Add square distance to sum of sqs.
					v_tmp[i] += itmp*itmp;

			//		.....Repeat for 'Left' Neighbour............
					if (!clsenb) {
						m_nbs_tmp[n + i+1] = 0;
					}
					v_tmp[i+1] += itmp*itmp;	
				}

	//		*********************************************************		
	//		--End of inner loop
	//		*********************************************************
			}


	//		--Calculate root sum of squared distance and add to root sum of square distance from other
	//		  dimensional orderings. First and last entry only have one neighbour, so to avoid changing the
	//		  denominator later on, we double these. Minimum value to add is 1 to avoid DIV/0 Errors.
			for(int i = 1; i < n-1; i++)
			{
				dtmp = v_tmp[id_order[i]];
				if (dtmp < 1) {
					v_ss[i] +=1;
				}
				else {
					v_ss[i] += sqrt(dtmp);
				}
			}
	//		(first and last - separately to avoid p*n 'if' comparisons in the for loop.)
			dtmp = v_tmp[id_order[0]];
			(dtmp < 1) ? v_ss[0] += 2 : v_ss[0] += 2*sqrt(dtmp);
			dtmp = v_tmp[id_order[n-1]];
			(dtmp < 1) ? v_ss[n-1] += 2 : v_ss[n-1] += 2*sqrt(dtmp);

	//		Clear sum of squares array to start from scratch for next dim.
			for(int i = 1; i <n; i++) {
				v_tmp[i] = 0;
			}

	//		Input neighbours into output array. Need to reorder based on original id now.
			for(int i = 0; i < n; i++)
			{
				m_nbs[2*d*n + i] = m_nbs_tmp[id_order[i]];
				m_nbs[(2*d+1)*n + i] = m_nbs_tmp[id_order[i] + n];
			}
	
	//	*********************************************************		
	//	--End of outer loop
	//	*********************************************************
		}
		
	
	//	*****Final operation: Density Calculation.***************
	//	--Compute density = p/n. If density > phi, then point is interior; type = 1.
		
		for(int i = 0; i < n; i++)
		{
			if (p/v_ss[i] > phi) {
				v_type[i] = 1;
			}
			else {
				v_type[i] = 0;
			}
		}

	delete[] col;
	delete[] new_order;
	delete[] id_order;
	delete[] v_id_copy;
	delete[] m_pos_copy;
	delete[] m_nbs_tmp;
	delete[] v_tmp;
	return;
}
}