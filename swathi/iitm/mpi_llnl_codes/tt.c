#include<stdio.h>
#include<mpi.h>

#define ROOT 0

main(int argc,char *argv[])	
{
int num,rank;
int a,b;

MPI_Status status;
MPI_Request handle;

MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&num);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);


if(rank==0)
{
a=10;
MPI_Bcast(&a,1,MPI_INT,ROOT,MPI_COMM_WORLD);
}
MPI_Barrier(MPI_COMM_WORLD);

printf("a:%d from process%d\n",a,rank);



MPI_Finalize();

}