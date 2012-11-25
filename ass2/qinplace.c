/**
 * Parallelization of quicksort for shared
 * address space platforms using POSIX threads.
 *
 * Author: Paolo Boschini
 * Author: Marcus Ihlar
 */

#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

// prototypes
void serial_quicksort(double*,int,int);
void* paral_quicksort(void*);
int partition(double*,int,int,int);
void swap(double*,double*);
void print_array(double*,int);
int is_sorted(double*,int);
int choose_pivot(double*,int,int,int);
int timer();

struct thread_data {
    int depth, left, right;
    double * array;
};

pthread_attr_t attr;
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

    
    int n_thread = 0;
    for(i = 1; i <= max_level; ++i)
        n_thread += pow(2,i);
    printf("\nMax level is %d, n_threads is %d\n\n",max_level,n_thread);


    double *a = (double*)malloc(size * sizeof(double));
    if (a == NULL) {
        printf("Something went wrong with memory allocation...\n");
        exit(-1);
    }

    srand(seed);
    for(i = 0; i < size; ++i)
        a[i] = ((double)rand()*1000/(double)(RAND_MAX));
    printf("Generated array was:\n");
    print_array(a, size);


    // init thread
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    void *status;

    struct thread_data data;
    data.depth = 0;
    data.array = a;
    data.left = 0;
    data.right = size-1;

    pthread_t t;

    int ttime=timer();
    pthread_create(&t, &attr, paral_quicksort, (void *)&data);
    pthread_join(t, &status);
    ttime=timer()-ttime;

    
    printf("Sorted array was:\n");
    print_array(a, size);
    printf("Time: %f\n",ttime/1000000.0);
    printf(is_sorted(a,size) ? "Success!!!\n" : "Fail!!!\n");

    //clear threads
    pthread_attr_destroy(&attr);
    pthread_exit(NULL);

    return 0;
}

void* paral_quicksort(void *arg) {
    pthread_t tL, tR;
    struct thread_data leftData;
    struct thread_data rightData;
    void *status;

    struct thread_data *myData = (struct thread_data*)arg;
    int depth = myData -> depth;
    double *array = myData -> array;
    int left = myData -> left;
    int right = myData -> right;

    // quicksort
    if(left < right) {
        int pivotNewIndex = partition(array,left,right,choose_pivot(array,left,(left+right)/2,right));

        // choose either serial or parallel quicksort
        if(depth < max_level) {
            leftData.depth = ++depth;
            leftData.array = array;
            leftData.left = left;
            leftData.right = pivotNewIndex-1;

            rightData.depth = depth;
            rightData.array = array;
            rightData.left = pivotNewIndex+1;
            rightData.right = right;

            int l=pthread_create(&tL, &attr, paral_quicksort, (void *)&leftData);
            int r=pthread_create(&tR, &attr, paral_quicksort, (void *)&rightData);
            // printf("l = %d, r = %d\n",l,r);
            
            pthread_join(tL, &status);
            pthread_join(tR, &status);
        }
        else {
            serial_quicksort(array,left,pivotNewIndex-1);
            serial_quicksort(array,pivotNewIndex+1,right);
        }
    }
    pthread_exit(NULL);
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