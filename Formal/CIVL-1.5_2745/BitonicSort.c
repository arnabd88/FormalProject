
//#include <iostream>
//#include <fstream>
//#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>


#define NUM_THREADS 8
#define ARRAYSIZE 5

//using namespace std;
//using std::string ;

typedef struct DiArrays {
   //int *b1 ;
   //int *b2 ;
   int N1  ;
   int N2  ;
} Darr ;

int numProcs, rank ;
MPI_Comm comm ;

Darr reverse_split( int* b, size_t N)
{
  Darr b12 ;
  b12.N1 = N-N/2 ;
  b12.N2 = N/2 ;
  int temp ;
  int residue = N%2 ;
  //#pragma omp parallel for
  for(int i=0; i<N/2; i++)
  {
    if(b[i] < b[N/2 + residue + i])
	{
	  temp = b[i] ;
	  b[i] = b[N/2 + residue + i] ;
	  b[N/2 + residue + i] = temp ;
	}
	//#pragma omp critical
	//{
	if(residue==1 && b[N/2] < b[N/2+residue+i])
	{
	  temp = b[N/2] ;
	  b[N/2] = b[N/2 + residue+i];
	  b[N/2 + residue+i] = temp ;
	}
	//}
  }
  if(residue==1 && b[N/2] < b[N/2+1])
  {
     temp = b[N/2];
	 b[N/2] = b[N/2+1] ;
	 b[N/2+1] = temp ;
  }
  return b12 ;
}

Darr split( int* b, size_t N)
{
  Darr b12 ;
  b12.N1 = N/2 ;
  b12.N2 = N - N/2 ;
  //b12.b1 = b ;
  //b12.b2 = &b[N/2] ;
  int temp ;
  int residue = N%2 ;
 // #pragma omp parallel for
  for(int i=0 ; i<N/2; i++)
  {
    if(b[i] > b[N/2 + residue+ i])
	{
	   temp = b[i] ;
	   b[i] = b[N/2 + residue + i] ;
	   b[N/2 + residue + i] = temp ;
	 //  #pragma omp critical
	   //{
	   if(residue==1 && b[i] >  b[N/2])
	   {
	     temp = b[i] ;
		 b[i] = b[N/2] ;
		 b[N/2] = temp ;
	   }
	   //}
	}
  }
  if(residue==1 && (b[N/2-1] > b[N/2]))
  {
     temp = b[N/2-1];
	 b[N/2-1] = b[N/2] ;
	 b[N/2] = temp ;
  }

  return b12 ;
}

int ReverseSortBitonic( int* b, int N)
{
  if(N==1) return 1;
  Darr d_split = reverse_split(b,N);
 // #pragma omp parallel
  //{
  ReverseSortBitonic(b, d_split.N1);
  ReverseSortBitonic(&b[d_split.N1], d_split.N2);
  //}
  //#pragma omp barrier

  return 1 ;
}

int SortBitonic( int* b, int N)
{
  if (N==1) return 1 ;
  Darr d_split = split(b, N);
 // #pragma omp parallel
 // {
  SortBitonic(b, d_split.N1);
  SortBitonic(&b[d_split.N1], d_split.N2);
 // }
 // #pragma omp barrier

  return 1 ;

}

void MergeBitonic( int* b1, int*b2, int N1, int N2)
{
  SortBitonic(b1, N1);
  ReverseSortBitonic(b2, N2);

  //for(int i=0; i<N1; i++)  b1[i] = b_1[i] ;
  for(int i=0; i<N2; i++)  b1[N1+i] = b2[i] ;
}

int* MakeBitonic( int* seq, int N)
{
  int* seq1 = seq ;
  int* seq2 = &seq[N/2];
  if(N > 4)
  {
  // #pragma omp parallel
  // {
    MakeBitonic(seq1, N/2);
	MakeBitonic(seq2, N-N/2);
  // }
  // #pragma omp barrier
  }
   MergeBitonic(seq1, seq2, N/2, N-N/2);
   return seq1 ;


}

void bitonicSearch(int* buffer, size_t N, size_t size)
 {
    //omp_set_num_threads(NUM_THREADS) ;
	//--- Make bitonic sequences for the sequential parts with processors
	int r = N/numProcs ; //--- Number of elements per process
	if(r > 1)
	MakeBitonic(&buffer[rank*r], r);
	 // if(rank==1)
	 // {
	 //   cout << "Make Bitonic" << endl ;
	 //   for(int k=0; k<r; k++) cout << buffer[rank*r +k] << " " ;
	 // }

	 //--- Start distributed transaction and exchange ---
	 int dim = ceil(log2(numProcs));
	 int comm_partner ;
	 int mysize = r;
	 int partnersize = r ;
	 //cout << "Dim = " << dim << endl ;
	 //-- for d==0
     for(int j=0; j<dim; j++)
	 {
	    comm_partner = rank ^ (1 << j) ;
	//    if(rank==0) cout << "Rank = " << rank << "; RankPArtner = " << comm_partner << endl ;
	    //--- Send a piece of the array to partner
	    MPI_Status status ;
	    if( rank < numProcs && comm_partner < numProcs)
	    {
	      if(rank < comm_partner)
	      {
	         MPI_Send(&mysize, 1, MPI_INT, comm_partner, 1, MPI_COMM_WORLD);
	         MPI_Recv(&partnersize, 1, MPI_INT, comm_partner, 1, MPI_COMM_WORLD, &status);
             MPI_Send(&buffer[rank*r], mysize, MPI_INT, comm_partner, 1, MPI_COMM_WORLD);
	         MPI_Recv(&buffer[comm_partner*r], partnersize, MPI_INT, comm_partner, 1 , MPI_COMM_WORLD, &status);
	      }
	      else if(rank > comm_partner)
	      {
	        MPI_Recv(&partnersize, 1, MPI_INT, comm_partner, 1, MPI_COMM_WORLD, &status);
	        MPI_Send(&mysize, 1, MPI_INT, comm_partner, 1, MPI_COMM_WORLD);
	        MPI_Recv(&buffer[comm_partner*r], partnersize, MPI_INT, comm_partner, 1, MPI_COMM_WORLD, &status);
	        MPI_Send(&buffer[rank*r], mysize, MPI_INT, comm_partner, 1, MPI_COMM_WORLD);
	      }
	       //--- Do the compare and exchange ---
	       if(rank < comm_partner)
	       {
	         int* seq1 = &buffer[rank*r] ;
	         int* seq2 = &buffer[comm_partner*r];
	         MergeBitonic(seq1, seq2, mysize, partnersize);
			 mysize += partnersize ;
	         // if(rank==0)
	         // {
	         //   cout << "After first level of merge" << endl ;
	         //   for(int k=0; k<mysize; k++)  cout << seq1[rank*r + k] << " " ;
	         // }
	       }
	    }
	}
	if(rank==0)
	{
	  SortBitonic(buffer, N);
	        //   cout << "\n\nLast level of merge" << endl ;
			printf("\n\n Last level of merge \n");
	        for(int k=0; k<N; k++) // cout << buffer[rank*r + k] << " " ;
			    printf(" %d ", buffer[rank*r + k]) ;
	}
	  MPI_Bcast(buffer, N, MPI_INT, 0, comm);
	  MPI_Barrier(comm);
	 

	
 }


void BitonicDriver( int N )
 {

    int* buffer = (int*)malloc(sizeof(int)*N) ;
	int r = N / numProcs ;
	int last = 29 ;
	//--- Initialize with random numbers -----
	if(rank==0)
	{
	printf("--- Unsorted Input-Array ------\n");
	#pragma parallel for
	//for(int k=0; k<N; k++)  buffer[k] = rand()%243 ;
	for(int k=0; k<N; k++) {
	//	buffer[k] = (last%(numProcs%(k+1)+1)) << k%numProcs ; 
	    buffer[k] = (k*71)%(k+7) ;
	}

    
	for(int k=0; k<N; k++)  //cout << buffer[k] << " " ;
	   printf("%d  ", buffer[k]);
	}
	MPI_Bcast(buffer, N, MPI_INT, 0, comm);
	MPI_Barrier(comm);
	bitonicSearch(buffer, N, sizeof(int));
	if(rank==0)
	   printf("Merging Completed Successfully......");
 }


 int main(int argc, char* argv[])
 {
   //int numProcs, rank ;
   int N=0;
   //---- Initializations ----//
   MPI_Init(&argc, &argv);
   comm = MPI_COMM_WORLD ;
   MPI_Comm_size(comm, &numProcs);
   MPI_Comm_rank(comm, &rank);
   double start = MPI_Wtime();
   double end, WCT ;

   N = ARRAYSIZE ;

   if(argc !=2 || argv[1]=="-h")
   {
     if(rank==0)
	 {
	   printf("Incorrect options- See below help menu\n");
	   printf("Format: mpirun ./a.out <N/P(num keys per processor ratio>\n");
	 }
   }
   else {
    // sscanf(argv[1], "%d", &N);
	 if(N==0)
	 {
	   if(rank==0)
	   	printf("Invalid problem size .... Exiting\n");
	 }
	 else
	 {
	   BitonicDriver(N*numProcs) ;
	   WCT = MPI_Wtime() - start ;
	   if(rank==0)
	   {
	       printf("Wall Clock Time = %f\n", WCT);
	   }
	 }
   }


  //BitonicDriver(ARRAYSIZE) ;
  MPI_Finalize();
 }
