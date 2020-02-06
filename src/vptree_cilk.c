#include <stdio.h>
#include <stdlib.h>
#include "../inc/vptree.h"
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>

#define DISTANCE_THRESHOLD 10000
#define SPLIT_THRESHOLD 20000

// Functions for vptree.h
vptree * buildvp(double *X, int n, int d);
vptree * getInner(vptree * T);
vptree * getOuter(vptree * T);
double getMD(vptree * T);
double * getVP(vptree * T);
int getIDX(vptree * T);
//newnode function
vptree * newnode(double *data, int *index, int n, int d);
void  euclidean_distance(double *data,int *index, double *dist, int n, int d);
//functions for quickselect to find median value
void swap(double *x, double *y);
void swap_int(int *x, int *y);
int partition (double *arr,int *index, int l, int r);
double quickselect(double *arr,int *index, int l,int r, int k);












/*              VANTAGE POINT TREE IMPLEMENTATION

    Newnode function will create a new node following the structure:
    T.vp    : vantage point 
    T.md    : median distance of vantage point to the other
    T.idx   : the index of the vantage point in the original set
    T.inner : inner subtree
    T.outer : outer subtree
    The structure is declared in vptree.h
    Newnode will fill the vp, md , idx variables and then recursively call itself
    to calculate T.inner and T.outer

*/

vptree * newnode(double *data, int *index, int n, int d)
{
    //allocate memory for our new node
    vptree *node = malloc(sizeof(vptree));

    //0 case
    if (n<=0){
        free(node);
        return NULL;
    }
    // Recursion guard(1 point case)
    if (n == 1)
    {
        node->vp = &data[index[n-1]*d];
        node->md = 0;
        node->idx = index[0];
        node->inner = NULL;
        node->outer = NULL;

        return node;
    }

    //For n > 1 :

    //STEP 1: Choose last point as Vantage point and calculte euclidean distances
    node->vp = &data[index[n-1]*d];
    node->idx = index[n-1];

    //allocate array to save our distances(n-1 because last point is vantage point)
    double *dist=(double *)malloc((n-1)*sizeof(double));
    //call euclidean distance function
    euclidean_distance(data, index, dist, n, d);

    //STEP 2: Find Median
    //We always choose (n+1)/2 as our median
    //In our case n/2 because dist array is n-1 size
    int median = floor(n/2);


    //n-1 because our last point is vantage point, n-2 because our qselect works for 
    //an array of [l..r] size
    node->md = quickselect(data,index,0, n-2, median);
    //free memory
    free(dist);

    //STEP 3: Find the arrays left and right(inner / outer) of median
    // We avoid to add a cilk_spawn attribute to the second recursive call to vpt()
    // because that will create an empty continuation
    if (n>=SPLIT_THRESHOLD) 
    {
        node->inner =cilk_spawn newnode(data, index, median, d);
        node->outer = newnode(data, &index[median], n-median-1, d); 
    }
    else 
    {
        node->inner = newnode(data, index, median, d);
        node->outer = newnode(data, &index[median], n-median-1, d); 
    }
    //cilc sync
    cilk_sync;
    return node;

   

}

/*             EUCLIDEAN_DISTANCE IMPLEMENTATION
    We implement a function that will return an array of euclidean distances of our 
    data points from a point of choice.(in VPT case, the vantage point)

*/

void  euclidean_distance(double *data,int *index, double *dist, int n, int d){

    //cilk parallelism
    if (n>=DISTANCE_THRESHOLD)
    {
        cilk_for(int i = 0; i < n-1; i++){
            dist[i] = 0;
            for (int j = 0; j < d; j++){
                dist[i] += pow(data[index[n-1]*d+j] - data[index[i]*d+j], 2);
        }
        dist[i] = sqrt(dist[i]);  
    }
    }
    else
    {
    //n-1 because the last point is our vantage point, we dont want to compare it with itself.
        for (int i = 0; i < n-1; i++){
            dist[i] = 0;
            for (int j = 0; j < d; j++){
                dist[i] += pow(data[index[n-1]*d+j] - data[index[i]*d+j], 2);
            }
            dist[i] = sqrt(dist[i]);
        }

    }
}


/* QUICKSELECT IMPLEMENTATION
	-KTHSMALLEST
    -SWAP
    -PARTITION
    
slighty modified version of c++ implementation found here:
https://www.geeksforgeeks.org/quickselect-algorithm/
https://en.wikipedia.org/wiki/Quickselect(pseudocode)

*/

// This function returns k'th smallest element in arr[l..r] using QuickSort.
double quickselect(double *arr,int *index, int l,int r, int k)
{

     // If k is smaller than number of  
    // elements in array 
    if (k > 0 && k <= r - l + 1) { 
  
          
        int pivot = partition(arr,index, l, r); 
  
        // If position is same as k 
        if (pivot - l == k - 1) 
            return arr[pivot]; 
  
        // If position is more, recur  
        // for left subarray 
        if (pivot - l > k - 1) 
            return quickselect(arr,index, l, pivot - 1, k); 
  
        // Else recur for right subarray 
        return quickselect(arr,index, pivot + 1, r, k - pivot + l - 1); 


   
    } 
    return -1;

}

// Standard partition process of QuickSort().
int partition (double *arr,int *index, int l, int r)
{
    double x = arr[r];
    int i = l;

    for (int j = l; j <= r - 1; j++)
    {
        if (arr[j] <= x)
        {
            swap(&arr[i], &arr[j]);
            swap_int(&index[i], &index[j]);
            i++;            
        }
    }
    swap(&arr[i], &arr[r]);
    swap_int(&index[i], &index[r]);
    return i;
}

// Swaps two elements(double)
void swap(double *x, double *y)
{
    double temp = *x;
    *x = *y;
    *y = temp;
}
//swap two elements(int) for index array
void swap_int(int *x, int *y)
{
    int temp = *x;
    *x = *y;
    *y = temp;
}





// ACCESSORS

vptree * buildvp(double *X, int n, int d)
{
    //index array to keep the original indexes during tree building
    int *index = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++)
        index[i] = i;
    vptree *root = newnode(X, index, n, d);
    return root;
}


vptree * getInner(vptree * T){
    return T->inner;
}


vptree * getOuter(vptree * T){
    return T->outer;
}


double getMD(vptree * T){
    return T->md;
}


double * getVP(vptree * T){
    return T->vp;
}

int getIDX(vptree * T){
    return T->idx;
}


