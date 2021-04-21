#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

int time_int_ring(long Nrepeat, long Nsize, MPI_Comm comm, double* time) {
	int rank, size; // process's rank, # of total processes
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	int sum = 0; // local variable to keep track of sum

	MPI_Barrier(comm);

	double tt = MPI_Wtime();

	for (long repeat = 0; repeat < Nrepeat; ++repeat) {
		// printf("process %d started repeat %d\n", rank, repeat);

		if (rank > 0) {
			MPI_Recv(&sum, Nsize, MPI_INT, rank - 1, repeat, comm, MPI_STATUS_IGNORE);
			// printf("process %d receives %d at repeat %d\n", rank, sum, repeat);
			sum += rank;
		}

		MPI_Send(&sum, Nsize, MPI_INT, (rank + 1) % size, repeat, comm);
		// printf("process %d sends sum = %d at repeat %d\n", rank, sum, repeat);

		if (rank == 0) {
			MPI_Recv(&sum, Nsize, MPI_INT, size - 1, repeat, comm, MPI_STATUS_IGNORE);
			// printf("process %d receives %d at repeat %d\n", rank, sum, repeat);
		}

		//printf("process %d finished repeat %d\n", rank, repeat);
	}

	tt = MPI_Wtime() - tt;
	*time = tt;
	return sum;
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	
	if (argc < 2) {
		printf("Usage: ./int_ring <# of repeat>\n");
		exit(1);
	}

	int Nrepeat = atoi(argv[1]);

	int rank;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);

	double tt;
	int result = time_int_ring(Nrepeat, 1, comm, &tt);
	if (rank == 0) {
		printf("int ring latency: %e ms\n", tt/Nrepeat * 1000);
		printf("result: %d\n", result);
	}

	MPI_Finalize();
}

