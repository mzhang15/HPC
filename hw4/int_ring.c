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

void time_array_ring(long Nrepeat, long Nsize, MPI_Comm comm, double* time) {
	int rank, size; // process's rank, # of total processes
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	char* msg = (char*) malloc(Nsize);
	for (int i = 0; i < Nsize; ++i) msg[i] = 42;

	MPI_Barrier(comm);

	double tt = MPI_Wtime();

	for (long repeat = 0; repeat < Nrepeat; ++repeat) {
		// printf("process %d started repeat %d\n", rank, repeat);

		if (rank > 0) {
			MPI_Recv(msg, Nsize, MPI_CHAR, rank - 1, repeat, comm, MPI_STATUS_IGNORE);
			// printf("process %d receives %d at repeat %d\n", rank, sum, repeat);
		}

		MPI_Send(msg, Nsize, MPI_CHAR, (rank + 1) % size, repeat, comm);
		// printf("process %d sends sum = %d at repeat %d\n", rank, sum, repeat);

		if (rank == 0) {
			MPI_Recv(msg, Nsize, MPI_CHAR, size - 1, repeat, comm, MPI_STATUS_IGNORE);
			// printf("process %d receives %d at repeat %d\n", rank, sum, repeat);
		}

		//printf("process %d finished repeat %d\n", rank, repeat);
	}

	tt = MPI_Wtime() - tt;
	*time = tt;
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	
	if (argc < 2) {
		printf("Usage: ./int_ring <# of repeat>\n");
		exit(1);
	}

	int Nrepeat = atoi(argv[1]);

	int rank, size;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	double tt;
	int result = time_int_ring(Nrepeat, 1, comm, &tt);
	if (rank == 0) {
		printf("int ring latency: %e ms\n", tt/Nrepeat * 1000);
		printf("result: %d\n", result);
	}

	time_array_ring(Nrepeat, 2000000, comm, &tt);
	if (rank == 0) {
		printf("ring bandwidth: %e GB/s\n", (Nsize * Nrepeat *size)/tt/1e9);
	}

	MPI_Finalize();
}

