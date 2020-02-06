#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include<string.h>
#include "../inc/vptree.h"

void deltree(vptree* node,int d) ;

void deleteTree(vptree** node_ref,int d) 
{ 
    
  deltree(*node_ref,d); 
  
  *node_ref = NULL; 
} 


void deltree(vptree* node,int d) 
{ 
    if (node == NULL) return; 
    
    /* first delete both subtrees */
    deltree(getInner(node),d); 
    deltree(getOuter(node),d); 
    free(node); 
} 

int main()
{

    int n = 10;
    int d = 1;
    struct timeval start, end;
    double cpu_time_used;
    double hold_temp = 1;
    double time_increase;

 
    FILE *f = fopen("xronoi.csv", "w+");
    fprintf(f, "\nDimensions, Time(sec), Increase");

    for (int i=10; i<=200000; i*=5){
    	int b = 2;
		double *X = malloc(i * b * sizeof(double));
    	for(int k = 0; k <i * b; k++)
        	X[k] = rand();    			 
    			

    	gettimeofday (&start, NULL);
    	vptree *tree = buildvp(X, i, b);
    	gettimeofday (&end, NULL);
    	deleteTree(&tree,d);
    	free(X);
    	cpu_time_used = (double)((end.tv_usec - start.tv_usec)/1.0e6 + end.tv_sec - start.tv_sec);
    	time_increase = (cpu_time_used / hold_temp);

    	fprintf(f,"\n %dx%d, %f, %f ",i,b,cpu_time_used, time_increase);
    	printf("\nDimensions: %dx%d Time: %f sec Increase: %f ",i,b,cpu_time_used, time_increase);
    	hold_temp = cpu_time_used;

	}
	hold_temp = 1;


	for (int i=2;i<=2500; i*=2){
    	int b = 100000;
		double *X = malloc(b * i * sizeof(double));
    	for(int k = 0; k <b * i; k++)
        	X[k] = rand();    			 
    			

    	gettimeofday (&start, NULL);
    	vptree *tree = buildvp(X, i, b);
    	gettimeofday (&end, NULL);
    	deleteTree(&tree,d);
    	free(X);
    	cpu_time_used = (double)((end.tv_usec - start.tv_usec)/1.0e6 + end.tv_sec - start.tv_sec);
    	time_increase = (cpu_time_used / hold_temp);
    	fprintf(f,"\n %dx%d, %f, %f ",b,i,cpu_time_used, time_increase);
    	printf("\nDimensions: %dx%d,Time: %f sec Increase: %f ",b,i,cpu_time_used, time_increase);
    	hold_temp = cpu_time_used;

	}

	fclose(f);
    
			

 

    return 0;

}