

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "omp.h"

typedef struct DiArray {
  int n1;
  int n2;
 } Darr ;

Darr display() {
  printf("here I am\n");
  Darr d1 ;
  d1.n1 = 5 ;
  d1.n2 = 10 ;
  return d1 ;
}

int Compute( MPI_Comm comm, int rank, int numProcs ) {
  int Sum ;
  int send_partner , recv_partner ;
  MPI_Status stat ;
  if(rank == numProcs -1) send_partner = 0 ; //--Loop around for the last one--
  else                    send_partner = rank+1 ;
  if(rank == 0) recv_partner = numProcs - 1 ;
  else          recv_partner = rank - 1 ;
  
  if(rank==0) {
     Sum = rank ;
     MPI_Send(&Sum, 1, MPI_INT, send_partner, 1, comm);
	 MPI_Recv(&Sum, 1, MPI_INT, recv_partner, 1, comm, &stat);
  }
  else {
     MPI_Recv(&Sum, 1, MPI_INT, recv_partner, 1, comm, &stat);
	 Sum += rank ;
	 MPI_Send(&Sum, 1, MPI_INT, send_partner, 1, comm);
  }

  return Sum ;

}

int main(int argc, char* argv[])
{
  int numProcs, rank , Final_Sum;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm comm = MPI_COMM_WORLD ;
  Darr d1 = display();

  if(d1.n1 < d1.n2)
     printf("SUCCESS\n");

  Final_Sum = Compute( comm, rank, numProcs) ; // Sum at each individual stages
  if(rank==0)
  	printf("Final Sum on a Ring = %d\n", Final_Sum);

  

  MPI_Finalize() ;

}
