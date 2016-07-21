#include <stdio.h>
#include <mpi.h>

#define  MAX_ITERATIONS 10000
#define Root 0

double converge_criteria(double *X_Old, double *X_New, int n);

main(int argc, char** argv) 
{

  int n, NoofRows_Bloc;
  int Numprocs, MyRank;
  int irow, icol, index, Iteration, GlobalRowNo;

  MPI_Init(&argc, &argv); 
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &Numprocs);
  
  n=4;
  NoofRows_Bloc = n/Numprocs;

  double Matrix_A[n][n],Input_B[n],ARecv[NoofRows_Bloc][n],BRecv[NoofRows_Bloc];
  double X_New[n], X_Old[n], Bloc_X[NoofRows_Bloc];

  if(MyRank == Root)
  {

  for(irow = 0; irow < n; irow++)
  {

  for(icol = 0; icol < n; icol++)
    //  scanf("%lf", &Matrix_A[irow][icol]);
    if(irow==icol) {
      Matrix_A[irow][icol] = 4;
      Input_B[irow] = 1; }
    else 
      Matrix_A[irow][icol] = -1;
  }

  //for (irow = 0; irow<n; irow++)
  //scanf("%lf",&Input_B[irow]);

  }

  if(n % Numprocs != 0) 
  {
	  MPI_Finalize();
	  if(MyRank == 0)
 	  {
	  printf("No. of processors should be a factor of \n");
	  }
	  exit(-1);
  }  	

  MPI_Scatter (Matrix_A, NoofRows_Bloc * n, MPI_DOUBLE, ARecv, NoofRows_Bloc * n,MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Scatter (Input_B, NoofRows_Bloc, MPI_DOUBLE, BRecv, NoofRows_Bloc, MPI_DOUBLE, 0,MPI_COMM_WORLD);

  for(irow=0; irow<NoofRows_Bloc; irow++)
  Bloc_X[irow] = BRecv[irow];

  MPI_Allgather(Bloc_X, NoofRows_Bloc, MPI_DOUBLE, X_New, NoofRows_Bloc,MPI_DOUBLE, MPI_COMM_WORLD);

  for (irow = 0; irow<n; irow++)
  X_New[irow]=0;

  Iteration = 0;
  do{

   for(irow=0; irow<n; irow++)
   X_Old[irow] = X_New[irow];

   for(irow=0; irow<NoofRows_Bloc; irow++)
   {

   GlobalRowNo = (MyRank * NoofRows_Bloc) + irow;
   Bloc_X[irow] = BRecv[irow];
   index = irow * n;

   for(icol=0; icol<n; icol++)
   {
   if(icol!=GlobalRowNo)
   Bloc_X[irow] -= X_Old[icol] * ARecv[irow][icol];
   }

   Bloc_X[irow] = Bloc_X[irow] / ARecv[irow][GlobalRowNo];
   }

  MPI_Allgather(Bloc_X, NoofRows_Bloc, MPI_DOUBLE, X_New,NoofRows_Bloc, MPI_DOUBLE, MPI_COMM_WORLD);
   Iteration++;

  }while( (Iteration < MAX_ITERATIONS) && (converge_criteria(X_Old, X_New, n) >= 1.0E-24)); 
  
  if (MyRank == 0) 
  {

     printf ("\n");

     printf("Matrix Input_A \n");
     printf ("\n");
     for (irow = 0; irow < n; irow++)
     {
     for (icol = 0; icol < n; icol++)
     printf("%.3lf  ", Matrix_A[irow][icol]);
     printf("\n");
     }
     printf ("\n");
     printf("Matrix Input_B \n");
     printf("\n");
     for (irow = 0; irow < n; irow++)
     printf("%.3lf\n", Input_B[irow]);
    
     printf ("\n");
     printf("Solution vector \n");
     printf("Number of iterations = %d\n",Iteration);
     printf ("\n");
     for(irow = 0; irow < n; irow++)
     printf("%.12lf\n",X_New[irow]);
  }
  
  MPI_Finalize(); 
}

double converge_criteria(double *X_Old, double *X_New, int n)
{
   int  index;
   double Sum;

   Sum = 0.0;
   for(index=0; index<n; index++)
   Sum += (X_New[index] - X_Old[index])*(X_New[index]-X_Old[index]);

   return(Sum);
}





