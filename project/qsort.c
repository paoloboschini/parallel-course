/* Paolo Boschini, Programming of Parallel Computers, VT12 */
/* MPI Quicksort */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

// prototypes
void print_array(double*,int);
void fill_array(double*,int);
void swap(double*,double*);
int partition(double*,int,int,int);
void serial_quicksort(double*,int,int);
void paral_quicksort(int,int,MPI_Comm,double*,int,double*);
int choose_pivot(double*,int,int,int);
void merge(double*,double*,double*,int,int);
void countPartitions(int*,int*,double,double*,int);
int is_sorted(double*,int);
int timer();

int main (int argc, char *argv[]) {

  int wrank, nproc, size, i, level = 0, local_size;
  double *array, *local_array, *local_left, *local_right, median;

  MPI_Init(&argc, &argv);                /* initialize mpi               */
  MPI_Comm_size(MPI_COMM_WORLD, &nproc); /* get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &wrank); /* Get my number                */

  if(argc < 2 && wrank == 0) {
    printf(wrank == 0 ? "Usage: mpirun -np 4 ./qsort n (where n is the array size)\n" : "");
    MPI_Finalize(); 
    return -1;
  }

  size = atoi(argv[1]);
  
  // allocate whole array if on proc 0
  if (wrank == 0) {
    array = (double *)calloc(size,sizeof(double));
    fill_array(array,size);
  }

  // allocate local array for each processor
  local_size = size/nproc;
  local_array = (double*)calloc(local_size,sizeof(double));

  int ttime=timer();

  MPI_Scatter(array, local_size, MPI_DOUBLE, local_array, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
  serial_quicksort(local_array,0,local_size-1);
  paral_quicksort(0, nproc, MPI_COMM_WORLD, local_array, local_size, array);

  if(wrank == 0) {
    ttime=timer()-ttime;
    printf("procs: %d, Time: %f\n",nproc,ttime/1000000.0);
    //    printf ("\n");
    printf(is_sorted(array,size) ? "Success!!!\n" : "Fail!!!\n");
    print_array(array,size);
  }

  if (wrank == 0) free(array);
  MPI_Finalize(); 
  return 0;
}

void paral_quicksort(int level, int nproc, MPI_Comm comm, double *local_array, int local_size, double *global) {
  int i, rank, other_size, local_nproc;
  double median, *received_array, *new_local;
  MPI_Request request;
  MPI_Status status;

  // extract local rank and local nproc
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &local_nproc);

  // basecase
  if(level < log2(nproc)) {
    
    // choose pivot
    if(rank == 0)
      median = local_array[local_size/2];

    // broadcast pivot
    MPI_Bcast(&median, 1, MPI_DOUBLE, 0, comm);
    
    int lt_pivot = 0, gt_pivot = 0;
    countPartitions(&lt_pivot, &gt_pivot, median, local_array, local_size);

    if(rank < local_nproc/2) {
      // send data
      MPI_Isend(local_array+lt_pivot, gt_pivot, MPI_DOUBLE, rank+(local_nproc/2), 111, comm, &request);

      // receive data
      MPI_Probe(rank+(local_nproc/2), 222, comm, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &other_size);
      received_array = (double *)calloc(other_size,sizeof(double));
      MPI_Recv(received_array, other_size, MPI_DOUBLE, rank+(local_nproc/2), 222, comm, &status);
      MPI_Wait(&request, &status);

      // recalculate the new size for the merged array
      local_size = lt_pivot + other_size;
      new_local = (double *) calloc(local_size, sizeof(double));
    
      merge(local_array, received_array, new_local, lt_pivot, other_size);
    }
    if(rank >= local_nproc/2) {
      // send data
      MPI_Isend(local_array, lt_pivot, MPI_DOUBLE, rank-(local_nproc/2), 222, comm, &request);

      // receive data
      MPI_Probe(rank-(local_nproc/2), 111, comm, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &other_size);
      received_array = (double *)calloc(other_size,sizeof(double));    
      MPI_Recv(received_array, other_size, MPI_DOUBLE, rank-(local_nproc/2), 111, comm, &status);
      MPI_Wait(&request, &status);

      // recalculate the new size for the merged array 
      local_size = other_size + gt_pivot;
      new_local = (double *) calloc(local_size, sizeof(double));
    
      merge(received_array, local_array+lt_pivot, new_local, other_size, gt_pivot);      
    }

    free(local_array);
    free(received_array);

    // create new communicators for the recursive calls
    MPI_Comm new_comm;
    MPI_Comm_split(comm, rank < local_nproc/2, 0, &new_comm);
    paral_quicksort(level+1, nproc, new_comm, new_local, local_size, global);
  }
  else { // basecase reached
    int *displacements = (int*)calloc(nproc,sizeof(int));
    int *rcounts = (int*)calloc(nproc,sizeof(int));

    // send the sizes to proc 0
    MPI_Gather(&local_size, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // fix displacements
    displacements[0] = 0;
    for (i = 1; i < nproc; ++i)
      displacements[i] = displacements[i-1] + rcounts[i-1];

    // gather into global array
    MPI_Gatherv(local_array, local_size, MPI_DOUBLE, global, rcounts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    free(local_array);
  }
}

void countPartitions(int *lt_pivot, int *gt_pivot, double median, double *local_array, int local_size) {
  int i;
  // less than 4 elements, do a linear search
  if(local_size < 4) {
    for (i = 0; i < local_size; ++i) {
      if(local_array[i] <= median)
	++(*lt_pivot);
      else
	++(*gt_pivot);
    }
  }
  // more than 3 elements, start searching in the middle
  else {
    int el = local_size / 2;
    if(local_array[el] <= median) {
      for(i = el; i < local_size; ++i) {
	if(local_array[i] > median)
	  break;
      }
      *lt_pivot = i; *gt_pivot = local_size-i;
    }
    else if(local_array[el] > median) {
      for(i = el; i >= 0; --i) {
	if(local_array[i] <= median)
	  break;
      }
      *lt_pivot = i+1; *gt_pivot = local_size-i-1;
    }
  }
}

void merge(double *left, double *right, double* output, int n_left, int n_right) {
  int m = 0, n = 0, i;
  
  for(i = 0; i < n_left+n_right; ++i) {
    if(m < n_left && n < n_right) {
      if(left[m] <= right[n])
	output[i] = left[m++];
      else
	output[i] = right[n++];
    }
    else if(m < n_left)
      output[i] = left[m++];
    else if(n < n_right)
      output[i] = right[n++];
  }
}

void fill_array(double *array, int size) {
  int i;
  srand(29);
  for (i = 0; i < size; i++)
    array[i] = ((double)rand()*1000/(double)(RAND_MAX));
}

void print_array(double *array, int size) {
  int i;
  //  printf ("rank: %d\n",rank);
  for (i = 0; i < size; i++)
    printf (i == size-1 ? "%f\n" : "%f, ",array[i]);
}

void serial_quicksort(double *a,int left, int right) {
  if(left<right) {
    // Get lists of bigger and smaller items and final position of pivot
    //    int pivotNewIndex = partition(a,left,right,choose_pivot(a,left,(left+right)/2,right));
    int pivotNewIndex = partition(a,left,right,(left+right)/2);

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

int choose_pivot(double *a, int left,int mid,int right) {
  if(a[left] <= a[mid] && a[mid] <= a[right])
    return mid;
  if(a[left] <= a[right])
    return left;
  return right;
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

int timer() {
  struct timeval tv;
  gettimeofday(&tv, (struct timezone*)0);
  return (tv.tv_sec*1000000+tv.tv_usec);
}
