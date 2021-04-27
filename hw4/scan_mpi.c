#include <stdio.h>
#include <math.h>
#include <mpi.h>

// Scan A array and write result into prefix_sum array;
// use long data type to avoid overflow
// void scan_seq(long* prefix_sum, const long* A, long n) {
//   if (n == 0) return;
//   prefix_sum[0] = 0;
//   for (long i = 1; i < n; i++) {
//     prefix_sum[i] = prefix_sum[i-1] + A[i-1];
//   }
// }

// void scan_omp(long* prefix_sum, const long* A, long n) {
//   // TODO: implement multi-threaded OpenMP scan
//   int nthreads, block_len;
//   #pragma omp parallel
//   {
//     int tid = omp_get_thread_num();
//     if (tid == 0) {
//       nthreads = omp_get_num_threads();
//       block_len = n / nthreads;
//       //printf("Number of threads = %d\n", nthreads);
//     }

//     // Q: why adding a barrier here would cause compiling error?
//     // #pragma omp barrier
//     // printf("Thread %d is starting...\n", tid)
//     #pragma omp barrier
   
//     // 1. split A into chuncks and do local scan
//     #pragma omp for
//     for (long i = 0; i < n; ++i) {
//       if (i % block_len == 0) prefix_sum[i] = 0;
//       else prefix_sum[i] = prefix_sum[i - 1] + A[i - 1];   
//     }
//   }

//   // 2. calculate offsets
//   long offsets[nthreads];
//   long sum = 0;
//   for (long i = 0; i < n; ++i) {
//     if (i % block_len == 0) offsets[i / block_len] = sum;
//     sum += A[i];
//   }
  

//   #pragma omp parallel
//   {
//     int tid = omp_get_thread_num();
//     // 3. update prefixSum with offset
//     #pragma omp for
//     for (long i = 0; i < n; ++i) {
//       prefix_sum[i] += offsets[i / block_len];
//     }
//     //printf("Thread %d is done!\n", tid);
//   }
// }

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int rank, size;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
  if (argc < 2) {
	  printf("usage: ./scan-mpi <N>\n");
	  exit(1);
  }

  long N = atoi(argv[1]);
  long* A;
  int elements_per_proc = N / size;

    if (rank == 0) {
        // 1. initialize array with random numbers
        A = (long*) malloc(N * sizeof(long));
        time_t t;
        srand((int) time(&t));
        for (long i = 0; i < N; ++i) A[i] = rand();
    }

    // 2. split up vector and scatter to other processes
    long* sub_array = (long*) malloc(elements_per_proc * sizeof(long));
    
    MPI_Scatter(A, elements_per_proc, MPI_LONG, sub_array, 
                elements_per_proc, MPI_LONG, 0, MPI_COMM_WORLD);

    // 3. local scan on each process
    long* sub_prefix_array = (long*) malloc(elements_per_proc * sizeof(long));
    for (int i = 0; i < elements_per_proc; ++i) {
        if (i == 0) sub_prefix_array[i] = sub_array[i];
        else sub_prefix_array[i] = sub_prefix_array[i - 1] + sub_array[i];
    }

    // 4. share their offset with everybody else
    long* offsets = (long*) malloc(size * sizeof(long));
    MPI_Allgather(&sub_prefix_array[elements_per_proc - 1], 1, MPI_LONG, 
                  offsets, 1, MPI_LONG, MPI_COMM_WORLD);

    int offset = 0;
    for (int i = 1; i < rank; ++i) {
        offset += offsets[i - 1];
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    // 5. update vector entries
    printf("rank %d\n", rank);
    for (int i = 0; i < elements_per_proc; ++i) {
      sub_prefix_array[i] += offset;
      printf("%ld ", sub_prefix_array[i]);
    }
    printf("\n");

    MPI_Finalize();
    return 0;
}

// int main() {
//   long N = 100000000;
//   long* A = (long*) malloc(N * sizeof(long));
//   long* B0 = (long*) malloc(N * sizeof(long));
//   long* B1 = (long*) malloc(N * sizeof(long));
//   for (long i = 0; i < N; i++) A[i] = rand();

//   double tt = omp_get_wtime();
//   scan_seq(B0, A, N);
//   printf("sequential-scan = %fs\n", omp_get_wtime() - tt);

//   tt = omp_get_wtime();
//   scan_omp(B1, A, N);
//   printf("parallel-scan   = %fs\n", omp_get_wtime() - tt);

//   long err = 0;
//   for (long i = 0; i < N; i++) err = std::max(err, std::abs(B0[i] - B1[i]));
//   printf("error = %ld\n", err);

//   free(A);
//   free(B0);
//   free(B1);
//   return 0;
// }
