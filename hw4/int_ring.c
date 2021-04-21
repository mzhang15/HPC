#include <stdio.h>
#include <mpi.h>
#include <string.h>

unsigned int time_int_ring(long Nrepeat, long Nsize, MPI_Comm comm, double* time) {
	int rank;
	MPI_Comm_rank(comm, &rank);
	char* msg = (char*) malloc(Nsize);
	unsigned int sum = 0;

	MPI_Barrier(comm);

	double tt = MPI_Wtime();

	for (long repeat = 0; repeat < Nrepeat; ++repeat) {
		MPI_Status status;

		switch(rank) {
			case 0:
				if (repeat == 0) memcpy(msg, (char*)&sum,sizeof(unsigned int)); // put sum = 0 into msg
				else MPI_Recv(msg, Nsize, MPI_CHAR, 2, comm, &status);

				MPI_Send(msg, Nsize, MPI_CHAR, 1, comm);
				break;
			case 1:
				MPI_Recv(msg, Nsize, MPI_CHAR, 0, comm, &status);

				sum = *(unsigned int*)(msg);
				sum += 1;
				memcpy(msg, (char*)&sum,sizeof(unsigned int));

				MPI_Send(msg, Nsize, MPI_CHAR, 2, comm);
				break;
			case 2:
				MPI_Recv(msg, Nsize, MPI_CHAR, 1, comm, &status);

				sum = *(unsigned int*)(msg);
				sum += 2;
				memcpy(msg, (char*)&sum,sizeof(unsigned int));

				MPI_Send(msg, Nsize, MPI_CHAR, 0, comm);
				break
		}
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
	unsigned int result = time_int_ring(Nrepeat, 8, comm, &tt);
	if (rank == 0) printf("int ring latency: %e ms\n", *tt/Nrepeat * 1000);
	printf("result: %d\n", result);

	MPI_Finalize();
}

