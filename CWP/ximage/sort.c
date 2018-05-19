/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
SORT - Functions to sort arrays of data or arrays of indices

hpsort		sort an array of floats by the heap sort method
qkisort		sort an array of indices i[] so that 
		a[i[0]] <= a[i[1]] <= ... <= a[i[n-1]]
		uses the "quick sort" method
qkifind		partially sort an array of indices i[] so that the 
		index i[m] has the value it would have if the entire
		array of indices were sorted such that 
		a[i[0]] <= a[i[1]] <= ... <= a[i[n-1]]
		uses the "quick sort" method
qksort		sort an array of floats such that a[0] <= a[1] <= ... <= a[n-1]
		uses the "quick sort" method
qkfind		partially sort an array of floats  so that the element a[m] has
		the value it would have if the entire array were sorted
		such that a[0] <= a[1] <= ... <= a[n-1]
		uses the "quick sort" method

******************************************************************************
Function Prototypes:
void hpsort (int n, float a[]);
void qkisort (int n, float a[], int i[]);
void qkifind (int m, int n, float a[], int i[]);
void qksort (int n, float a[]);
void qkfind (int m, int n, float a[]);

******************************************************************************
hpsort:
Input:
n		number of elements in a
a		array[n] to be sorted

Output:
a		array[n] sorted

******************************************************************************
qkisort:
Input:
n		number of elements in array a
a		array[n] elements
i		array[n] indices to be sorted

Output:
i		array[n] indices sorted

******************************************************************************
qkifind:
Input:
m		index of element to be found
n		number of elements in array a
a		array[n] elements
i		array[n] indices to be partially sorted

Output:
i		array[n] indices partially sorted sorted

******************************************************************************
qksort:
Input:
n		number of elements in array a
a		array[n] containing elements to be sorted

Output:
a		array[n] containing sorted elements

******************************************************************************
qkfind:
Input:
m		index of element to be found
n		number of elements in array a
a		array[n] to be partially sorted

Output:
a		array[n] partially sorted


******************************************************************************
Notes:
hpsort:
The heap sort algorithm is, at worst, N log_2 N, and in most cases
is 20% faster.  Adapted from Standish.

qkisort, qkifind, qksort, qkfind:
n must be less than 2^NSTACK, where NSTACK is defined above.

qkisort:
This function is adapted from procedure quicksort by
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321;
the main difference is that recursion is accomplished
explicitly via a stack array for efficiency; also, a simple
insertion sort is used to sort subarrays too small to be
partitioned efficiently.

qkifind:
This function is adapted from procedure find by
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321.

qksort:
This function is adapted from procedure quicksort by
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321;
the main difference is that recursion is accomplished
explicitly via a stack array for efficiency; also, a simple
insertion sort is used to sort subarrays too small to be
partitioned efficiently.

qkfind:
This function is adapted from procedure find by Hoare 1961.

******************************************************************************
Reference:
hpsort:
Standish, T. A., Data Structure Techniques, p. 91.
See also Press, W. A., et al., Numerical Recipes in C.

quick sort:
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 12/26/89
*****************************************************************************/
/**************** end self doc ********************************/

#include "cwp.h"

void
hpsort (int n, float a[])
/*****************************************************************************
sort an array so that a[0] <= a[1] <= ... <= a[n-1]
******************************************************************************
Input:
n		number of elements in a
a		array[n] to be sorted

Output:
a		array[n] sorted
******************************************************************************
Notes:
Adapted from Standish, T. A., Data Structure Techniques, p. 91.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 12/26/89
*****************************************************************************/
{
	int kroot,klast,kparent,kchild;
	float aparent;
	
	/* initialize indices of root node and last node in heap */
	kroot = n/2-1;
	klast = n-1;
	
	/* loop until array is sorted */
	while (klast>0) {
	
		/* if in the heap building phase */
		if (kroot>=0) {
		
			/* set index of parent to be added to the heap */
			kparent = kroot;

			/* set value of parent to be added to the heap */
			aparent = a[kroot--];
		
		/* else, the tree is a heap and in the sorting phase */
		} else {
			
			/* set index of parent at which to begin sifting */
			kparent = 0;
	
			/* set parent value to last value in heap */
			aparent = a[klast];
			
			/* copy top of heap to sorted elements at end */
			a[klast--] = a[0];
		}

		/* sift parent down until greater than both children */
		for (kchild=kparent*2+1; kchild<=klast; kchild=kparent*2+1) {

			/* get index of greater of two children */
			if (kchild<klast && a[kchild+1]>a[kchild]) kchild++;

			/* if greater child is greater than parent */
			if (a[kchild]>aparent) {
			
				/* promote the greater child */
				a[kparent] = a[kchild];

				/* demote the parent */
				kparent = kchild;
			
			/* else, if parent is greater than children, break */
			} else 
				break;
		}
	
		/* set value of parent (which may have been demoted) */
		a[kparent] = aparent;
	}
}



#define NSTACK 50	/* maximum sort length is 2^NSTACK */
#define NSMALL 7	/* size of array for which insertion sort is fast */
#define FM 7875		/* constants used to generate random pivots */
#define FA 211
#define FC 1663

static void
qkipart (float a[], int i[], int p, int q, int *j, int *k)
/*****************************************************************************
quicksort partition (FOR INTERNAL USE ONLY):
take the value x of a random element from the subarray a[p:q] of
a[0:n-1] and rearrange indices in the subarray i[p:q] in such a way
that there exist integers j and k with the following properties:
  p <= j < k <= q, provided that p < q
  a[i[l]] <= x,  for p <= l <= j
  a[i[l]] == x,  for j < l < k
  a[i[l]] >= x,  for k <= l <= q
note that this effectively partitions the subarray with bounds
[p:q] into lower and upper subarrays with bounds [p:j] and [k:q]
******************************************************************************
Input:
a		array[p:q]
i		array[p:q] of indices to be rearranged
p		lower bound of subarray; must be less than q
q		upper bound of subarray; must be greater then p

Output:
i		array[p:q] of indices rearranged
j		upper bound of lower output subarray
k		lower bound of upper output subarray
******************************************************************************
Notes:
This function is adapted from procedure partition by
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int pivot,left,right,temp;
	float apivot;
	static long int seed=0L;
 
	/* choose random pivot element between p and q, inclusive */
	seed = (seed*FA+FC)%FM;
	pivot = p+(q-p)*(float)seed/(float)FM;
	if (pivot<p) pivot = p;
	if (pivot>q) pivot = q;
	apivot = a[i[pivot]];

	/* initialize left and right pointers and loop until break */
	for (left=p,right=q;;) {
		/*
		 * increment left pointer until either
		 * (1) an element greater than the pivot element is found, or
		 * (2) the upper bound of the input subarray is reached
		 */
		while (a[i[left]]<=apivot && left<q) left++;
 
		/*
		 * decrement right pointer until either
		 * (1) an element less than the pivot element is found, or
		 * (2) the lower bound of the input subarray is reached
		 */
		while (a[i[right]]>=apivot && right>p) right--;
 
		/* if left pointer is still to the left of right pointer */
		if (left<right) {

			/* exchange left and right indices */
			temp = i[left];
			i[left++] = i[right];
			i[right--] = temp;
		} 
		/* else, if pointers are equal or have crossed, break */
		else break;
	}
	/* if left pointer has not crossed pivot */
	if (left<pivot) {

		/* exchange indices at left and pivot */
		temp = i[left];
		i[left++] = i[pivot];
		i[pivot] = temp;
	}
	/* else, if right pointer has not crossed pivot */
	else if (pivot<right) {

		/* exchange indices at pivot and right */
		temp = i[right];
		i[right--] = i[pivot];
		i[pivot] = temp;
	}
	/* left and right pointers have now crossed; set output bounds */
	*j = right;
	*k = left;
}

static void
qkiinss (float a[], int i[], int p, int q)
/*****************************************************************************
quicksort insertion sort (FOR INTERNAL USE ONLY):
Sort a subarray of indices bounded by p and q so that
a[i[p]] <= a[i[p+1]] <= ... <= a[i[q]]
******************************************************************************
Input:
a		subarray[p:q] containing elements
i		subarray[p:q] containing indices to be sorted
p		lower bound of subarray; must be less than q
q		upper bound of subarray; must be greater then p

Output:
i		subarray[p:q] of indices sorted
******************************************************************************
Notes:
Adapted from Sedgewick, R., 1983, Algorithms, Addison Wesley, p. 96.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int j,k,ij;
	float aij;

	for (j=p+1; j<=q; j++) {
		for (ij=i[j],aij=a[ij],k=j; k>p && a[i[k-1]]>aij; k--)
			i[k] = i[k-1];
		i[k] = ij;
	}
}

void
qkisort (int n, float a[], int i[])
/*****************************************************************************
Sort an array of indices i[] so that 
a[i[0]] <= a[i[1]] <= ... <= a[i[n-1]]
******************************************************************************
Input:
n		number of elements in array a
a		array[n] elements
i		array[n] indices to be sorted

Output:
i		array[n] indices sorted
******************************************************************************
Notes:
n must be less than 2^NSTACK, where NSTACK is defined above.

This function is adapted from procedure quicksort by
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321;
the main difference is that recursion is accomplished
explicitly via a stack array for efficiency; also, a simple
insertion sort is used to sort subarrays too small to be
partitioned efficiently.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int pstack[NSTACK],qstack[NSTACK],j,k,p,q,top=0;

	/* initialize subarray lower and upper bounds to entire array */
	pstack[top] = 0;
	qstack[top++] = n-1;

	/* while subarrays remain to be sorted */
	while(top!=0) {

		/* get a subarray off the stack */
		p = pstack[--top];
		q = qstack[top];

		/* while subarray can be partitioned efficiently */
		while(q-p>NSMALL) {

			/* partition subarray into two subarrays */
			qkipart(a,i,p,q,&j,&k);

			/* save larger of the two subarrays on stack */
			if (j-p<q-k) {
				pstack[top] = k;
				qstack[top++] = q;
				q = j;
			} else {
				pstack[top] = p;
				qstack[top++] = j;
				p = k;
			}
		}
		/* use insertion sort to finish sorting small subarray */
		qkiinss(a,i,p,q);
	}
}

void
qkifind (int m, int n, float a[], int i[])
/*****************************************************************************
Partially sort an array of indices i[] so that the index i[m] has the
value it would have if the entire array of indices were sorted such that 
a[i[0]] <= a[i[1]] <= ... <= a[i[n-1]]
******************************************************************************
Input:
m		index of element to be found
n		number of elements in array a
a		array[n] elements
i		array[n] indices to be partially sorted

Output:
i		array[n] indices partially sorted sorted
******************************************************************************
Notes:
This function is adapted from procedure find by
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int j,k,p,q;

	/* initialize subarray lower and upper bounds to entire array */
	p = 0;  q = n-1;

	/* while subarray can be partitioned efficiently */
	while(q-p>NSMALL) {

		/* partition subarray into two subarrays */
		qkipart(a,i,p,q,&j,&k);

		/* if desired value is in lower subarray, then */
		if (m<=j)
			q = j;

		/* else, if desired value is in upper subarray, then */
		else if (m>=k)
			p = k;
		
		/* else, desired value is between j and k */
		else
			return;
	}
			
	/* completely sort the small subarray with insertion sort */
	qkiinss(a,i,p,q);
}



/*#define NSTACK 50	 maximum sort length is 2^NSTACK */
/*#define NSMALL 7	 size of array for which insertion sort is fast */
/*#define FM 7875	 constants used to generate random pivots */
/*#define FA 211	*/
/*#define FC 1663	*/


static void
qkpart (float a[], int p, int q, int *j, int *k)
/*****************************************************************************
quicksort partition (FOR INTERNAL USE ONLY):
Take the value x of a random element from the subarray a[p:q] of
a[0:n-1] and rearrange the elements in this subarray in such a way
that there exist integers j and k with the following properties:
  p <= j < k <= q, provided that p < q
  a[l] <= x,  for p <= l <= j
  a[l] == x,  for j < l < k
  a[l] >= x,  for k <= l <= q
Note that this effectively partitions the subarray with bounds
[p:q] into lower and upper subarrays with bounds [p:j] and [k:q].
******************************************************************************
Input:
a		array[p:q] to be rearranged
p		lower bound of subarray; must be less than q
q		upper bound of subarray; must be greater then p

Output:
a		array[p:q] rearranged
j		upper bound of lower output subarray
k		lower bound of upper output subarray
******************************************************************************
Notes:
This function is adapted from procedure partition by
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int pivot,left,right;
	float apivot,temp;
	static long int seed=0L;
 
	/* choose random pivot element between p and q, inclusive */
	seed = (seed*FA+FC)%FM;
	pivot = p+(q-p)*(float)seed/(float)FM;
	if (pivot<p) pivot = p;
	if (pivot>q) pivot = q;
	apivot = a[pivot];

	/* initialize left and right pointers and loop until break */
	for (left=p,right=q;;) {
		/*
		 * increment left pointer until either
		 * (1) an element greater than the pivot element is found, or
		 * (2) the upper bound of the input subarray is reached
		 */
		while (a[left]<=apivot && left<q) left++;
 
		/*
		 * decrement right pointer until either
		 * (1) an element less than the pivot element is found, or
		 * (2) the lower bound of the input subarray is reached
		 */
		while (a[right]>=apivot && right>p) right--;
 
		/* if left pointer is still to the left of right pointer */
		if (left<right) {
			/* exchange left and right elements */
			temp = a[left];
			a[left++] = a[right];
			a[right--] = temp;
		} 
		/* else, if pointers are equal or have crossed, break */
		else break;
	}
	/* if left pointer has not crossed pivot */
	if (left<pivot) {

		/* exchange elements at left and pivot */
		temp = a[left];
		a[left++] = a[pivot];
		a[pivot] = temp;
	}
	/* else, if right pointer has not crossed pivot */
	else if (pivot<right) {

		/* exchange elements at pivot and right */
		temp = a[right];
		a[right--] = a[pivot];
		a[pivot] = temp;
	}
	/* left and right pointers have now crossed; set output bounds */
	*j = right;
	*k = left;
}

static void
qkinss (float a[], int p, int q)
/*****************************************************************************
quicksort insertion sort (FOR INTERNAL USE ONLY):
Sort a subarray bounded by p and q so that
a[p] <= a[p+1] <= ... <= a[q]
******************************************************************************
Input:
a		subarray[p:q] containing elements to be sorted
p		lower bound of subarray; must be less than q
q		upper bound of subarray; must be greater then p

Output:
a		subarray[p:q] sorted
******************************************************************************
Notes:
Adapted from Sedgewick, R., 1983, Algorithms, Addison Wesley, p. 96.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int i,j;
	float ai;

	for (i=p+1; i<=q; i++) {
		for (ai=a[i],j=i; j>p && a[j-1]>ai; j--)
			a[j] = a[j-1];
		a[j] = ai;
	}
}

void
qksort (int n, float a[])
/*****************************************************************************
Sort an array such that a[0] <= a[1] <= ... <= a[n-1]
******************************************************************************
Input:
n		number of elements in array a
a		array[n] containing elements to be sorted

Output:
a		array[n] containing sorted elements
******************************************************************************
Notes:
n must be less than 2^NSTACK, where NSTACK is defined above.

This function is adapted from procedure quicksort by
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321;
the main difference is that recursion is accomplished
explicitly via a stack array for efficiency; also, a simple
insertion sort is used to sort subarrays too small to be
partitioned efficiently.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int pstack[NSTACK],qstack[NSTACK],j,k,p,q,top=0;

	/* initialize subarray lower and upper bounds to entire array */
	pstack[top] = 0;
	qstack[top++] = n-1;

	/* while subarrays remain to be sorted */
	while(top!=0) {

		/* get a subarray off the stack */
		p = pstack[--top];
		q = qstack[top];

		/* while subarray can be partitioned efficiently */
		while(q-p>NSMALL) {

			/* partition subarray into two subarrays */
			qkpart(a,p,q,&j,&k);

			/* save larger of the two subarrays on stack */
			if (j-p<q-k) {
				pstack[top] = k;
				qstack[top++] = q;
				q = j;
			} else {
				pstack[top] = p;
				qstack[top++] = j;
				p = k;
			}
		}
		/* use insertion sort to finish sorting small subarray */
		qkinss(a,p,q);
	}
}

void
qkfind (int m, int n, float a[])
/*****************************************************************************
Partially sort an array so that the element a[m] has the value it
would have if the entire array were sorted such that 
a[0] <= a[1] <= ... <= a[n-1]
******************************************************************************
Input:
m		index of element to be found
n		number of elements in array a
a		array[n] to be partially sorted

Output:
a		array[n] partially sorted
******************************************************************************
Notes:
This function is adapted from procedure find by
Hoare, C.A.R., 1961, Communications of the ACM, v. 4, p. 321.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int j,k,p,q;

	/* initialize subarray lower and upper bounds to entire array */
	p = 0;  q = n-1;

	/* while subarray can be partitioned efficiently */
	while(q-p>NSMALL) {

		/* partition subarray into two subarrays */
		qkpart(a,p,q,&j,&k);

		/* if desired value is in lower subarray, then */
		if (m<=j)
			q = j;

		/* else, if desired value is in upper subarray, then */
		else if (m>=k)
			p = k;
		
		/* else, desired value is between j and k */
		else
			return;
	}
			
	/* completely sort the small subarray with insertion sort */
	qkinss(a,p,q);
}
