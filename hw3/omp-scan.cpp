#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>

// Scan A array and write result into prefix_sum array;
// use long data type to avoid overflow
void scan_seq(long* prefix_sum, const long* A, long n) {
  if (n == 0) return;
  prefix_sum[0] = 0;
  for (long i = 1; i < n; i++) {
    prefix_sum[i] = prefix_sum[i-1] + A[i-1];
  }
}

void scan_omp(long* prefix_sum, const long* A, long n) {
  // TODO: implement multi-threaded OpenMP scan
  int nthreads, block_len;
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    if (tid == 0) {
      nthreads = omp_get_num_threads();
      block_len = n / nthreads;
      //printf("Number of threads = %d\n", nthreads);
    }

    // Q: why adding a barrier here would cause compiling error?
    // #pragma omp barrier
    // printf("Thread %d is starting...\n", tid)
    #pragma omp barrier
   
    // 1. split A into chuncks and do local scan
    #pragma omp for
    for (long i = 0; i < n; ++i) {
      if (i % block_len == 0) prefix_sum[i] = 0;
      else prefix_sum[i] = prefix_sum[i - 1] + A[i - 1];   
    }
  }

  // 2. calculate offsets
  long offsets[nthreads];
  long sum = 0;
  for (long i = 0; i < n; ++i) {
    if (i % block_len == 0) offsets[i / block_len] = sum;
    sum += A[i];
  }
  

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    // 3. update prefixSum with offset
    #pragma omp for
    for (long i = 0; i < n; ++i) {
      prefix_sum[i] += offsets[i / block_len];
    }
    //printf("Thread %d is done!\n", tid);
  }
}

int main() {
  long N = 100000000;
  long* A = (long*) malloc(N * sizeof(long));
  long* B0 = (long*) malloc(N * sizeof(long));
  long* B1 = (long*) malloc(N * sizeof(long));
  for (long i = 0; i < N; i++) A[i] = rand();

  double tt = omp_get_wtime();
  scan_seq(B0, A, N);
  printf("sequential-scan = %fs\n", omp_get_wtime() - tt);

  tt = omp_get_wtime();
  scan_omp(B1, A, N);
  printf("parallel-scan   = %fs\n", omp_get_wtime() - tt);

  long err = 0;
  for (long i = 0; i < N; i++) err = std::max(err, std::abs(B0[i] - B1[i]));
  printf("error = %ld\n", err);

  free(A);
  free(B0);
  free(B1);
  return 0;
}
