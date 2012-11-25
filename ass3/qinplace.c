/**
 * Quick-sort implementation in OpenMP using the task directive.
 *
 * Author: Paolo Boschini
 * Author: Marcus Ihlar
 */

#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

// prototypes
void serial_quicksort(double*,int,int);
void* paral_quicksort(double*,int,int,int);
int partition(double*,int,int,int);
void swap(double*,double*);
void print_array(double*,int);
int is_sorted(double*,int);
int choose_pivot(double*,int,int,int);
int timer();

int max_level;

int main (int argc, char const *argv[]) {
  int i;

  if(argc < 4) {
    printf("Usage: ./qinplace size_array seed max_level\n");
    exit(-1);
  }

  int size = atoi(argv[1]);
  int seed = atoi(argv[2]);
  max_level = atoi(argv[3]);

  double *a = (double*)malloc(size * sizeof(double));
  if (a == NULL) {
    printf("Something went wrong with memory allocation...\n");
    exit(-1);
  }

  srand(seed);
  for(i = 0; i < size; ++i)
    a[i] = ((double)rand()*1000/(double)(RAND_MAX));

  int ttime=timer();

  #pragma omp parallel {
    #pragma omp single nowait {
      paral_quicksort(a,0,size-1,0);
    }
  }

  ttime=timer()-ttime;

  printf("Time: %f\n",ttime/1000000.0);
  printf(is_sorted(a,size) ? "Success!!!\n" : "Fail!!!\n");

  return 0;
}


void* paral_quicksort(double *a, int left, int right, int depth) {
  if(left < right) {
    int pivotNewIndex = partition(a,left,right,choose_pivot(a,left,(left+right)/2,right));

    if(depth < max_level) {
      paral_quicksort(a,left,pivotNewIndex-1, depth+1);
      paral_quicksort(a,pivotNewIndex+1,right,depth+1);
    }
    else {
#pragma omp tasks
      {
#pragma omp task
	{
	  serial_quicksort(a,left,pivotNewIndex-1);
	}
#pragma omp task
	{
	  serial_quicksort(a,pivotNewIndex+1,right);
	}
      }
    }
  }
}

void serial_quicksort(double *a,int left, int right) {
  if(left<right) {
    // Get lists of bigger and smaller items and final position of pivot
    int pivotNewIndex = partition(a,left,right,choose_pivot(a,left,(left+right)/2,right));

    // Recursively sort elements smaller than the pivot
    serial_quicksort(a,left,pivotNewIndex-1);

    // Recursively sort elements at least as big as the pivot
    serial_quicksort(a,pivotNewIndex+1,right);
  }
}

int partition(double *a,int left,int right,int pivot_index) {
  int i;

  double pivotValue = a[pivot_index];
  // Move pivot to the end
  swap(&a[pivot_index],&a[right]);

  int storeIndex = left;

  for(i = left; i < right; ++i)
    if(a[i] < pivotValue)
      swap(&a[i],&a[storeIndex++]);

  // Move pivot to its final place
  swap(&a[storeIndex],&a[right]);

  return storeIndex;
}

void swap(double *x, double *y) {
  double temp = *x;
  *x = *y;
  *y = temp;
}

void print_array(double *a, int size) {
  int i;
  for(i = 0; i < size; ++i)
    printf(i==size-1 ? "%.2f" : "%.2f, ", a[i]);
  printf("\n\n");
}

int is_sorted(double *a, int size) {
  int i;
  for(i = 0; i < size-1; ++i)
    if(a[i] > a[i+1]) {
      printf("%f - %f\n", a[i],a[i+1]);
      return 0;
    }
  return 1;
}

int choose_pivot(double *a, int left,int mid,int right) {
  if(a[left] <= a[mid] && a[mid] <= a[right])
    return mid;
  if(a[left] <= a[right])
    return left;
  return right;
}

int timer() {
  struct timeval tv;
  gettimeofday(&tv, (struct timezone*)0);
  return (tv.tv_sec*1000000+tv.tv_usec);
}
