#include<algorithm>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <vector>

using namespace std;

double x[10]={1,2,3,4,5,6,7,8,9,0};
double y[10]={1,2,3,4,5,6,7,8,9,0};
double z[10]={1,2,3,4,5,6,7,8,9,0};
double Xmax=9, Xmin=0;
double Ymax=9, Ymin=0;
double Zmax=9, Zmin=0;
double re=1, DELTA=0;
int NUM=100;


//declarations of global variables such as x, y, z, re, DELTA, xmax, etc/ 

//CUDA kernel for calcHash
__global__ void calcHash(double *d_x, double *d_y, double *d_z, uint2 *d_particleHash,\
	int *d_NUM, double *d_Xmax, double *d_Xmin, double *d_re, double *d_DELTA, double *d_Ymin, double *d_Ymax, double *d_Zmax, double *d_Zmin){

  int ncx, ncy, ncz;
  ncx = int((*d_Xmax - *d_Xmin) / (*d_re + *d_DELTA)) + 1;     // Number of cells in x direction
  ncy = int((*d_Ymax - *d_Ymin) / (*d_re + *d_DELTA)) + 1;     // Number of cells in y direction
  ncz = int((*d_Zmax - *d_Zmin) / (*d_re + *d_DELTA)) + 1;     // Number of cells in z direction

  int k = threadIdx.x + blockIdx.x * blockDim.x;
  int m = *d_NUM;
  int *icell, *jcell, *kcell;

  icell[k] = int((d_x[k] - *d_Xmin) / (*d_re + *d_DELTA)) + 1;
  jcell[k] = int((d_y[k] - *d_Ymin) / (*d_re + *d_DELTA)) + 1;
  kcell[k] = int((d_z[k] - *d_Zmin) / (*d_re + *d_DELTA)) + 1;

  int cellNum = icell[k] + (jcell[k] - 1)*ncx + (kcell[k] - 1)*ncx*ncy;

  d_particleHash[k].x = cellNum;
  d_particleHash[k].y = k;

}


/*&CUDA kernel for findCellStart

__global__ void findCellStart(int *particleHash, int *cellStart, int *cellEnd){

  int k = 1 + threadIdx.x + blockIdx.x * blockDim.x; //take care of index shift in # of threads
  if (particleHash[k] != particleHash[k+1]){
    cellEnd[particleHash[k]] = k;
    cellStart[particleHash[k+1]] = k+1;
  }

}
//CUDA kernel for neighb[][] creation

__global__ void createNeighbourArrays(int **neighb, int *cellStart, int *cellEnd, int *particleHash){

  //check which type of cell it is and loop surrounding cells
  loop{
    add the particles in those cells as neighbours in the 2-d array
  }

}
*/


void neighbour_cuda(){

  int *d_cellEnd, *d_cellStart, *d_NUM;
  uint2 *particleHash, *d_particleHash;
  double *d_x, *d_y, *d_z, *d_Xmax, *d_Xmin, *d_Ymax, *d_Ymin, *d_Zmax, *d_Zmin, *d_re, *d_DELTA;

  int arrsizeint = NUM * sizeof(int);
  int sizeint = sizeof(int);
  int arrsizedouble = NUM * sizeof(double);
  int sizedouble = sizeof(double);
  int arruint2size = NUM * sizeof(uint2);

  particleHash = (uint2 *)malloc(arruint2size);


  cudaMalloc((void **)&d_particleHash, arruint2size);
  cudaMalloc((void **)&d_cellStart, arrsizeint);
  cudaMalloc((void **)&d_cellEnd, arrsizeint);
  cudaMalloc((void **)&d_x, arrsizedouble);
  cudaMalloc((void **)&d_y, arrsizedouble);
  cudaMalloc((void **)&d_z, arrsizedouble);
  cudaMalloc((void **)&d_Xmin, sizedouble);
  cudaMalloc((void **)&d_Xmax, sizedouble);
  cudaMalloc((void **)&d_Ymin, sizedouble);
  cudaMalloc((void **)&d_Ymax, sizedouble);
  cudaMalloc((void **)&d_Zmin, sizedouble);
  cudaMalloc((void **)&d_Zmax, sizedouble);
  cudaMalloc((void **)&d_re, sizedouble);
  cudaMalloc((void **)&d_DELTA, sizedouble);
  cudaMalloc((void **)&d_NUM, sizeint);

  //cudaMemcpy(d_particleHash, particleHash, arrsizeint, cudaMemcpyHostToDevice);
  //cudaMemcpy(d_cellStart, cellStart, arrsizeint, cudaMemcpyHostToDevice);
  //cudaMemcpy(d_cellEnd, cellEnd, arrsizeint, cudaMemcpyHostToDevice);
  cudaMemcpy(d_x, &x, arrsizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, &y, arrsizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_z, &z, arrsizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Xmin, &Xmin, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Xmax, &Xmax, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ymin, &Ymin, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ymax, &Ymax, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Zmin, &Zmin, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Zmax, &Zmax, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_re, &re, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_DELTA, &DELTA, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_NUM, &NUM, sizeint, cudaMemcpyHostToDevice);




  //copy x, y, z onto the device
  //run calc hash on the global array
  calcHash<<<5,5>>>(d_x, d_y, d_z, d_particleHash, d_NUM, d_Xmax, d_Xmin, d_re, d_DELTA, d_Ymin, d_Ymax, d_Zmax, d_Zmin);
  cudaMemcpy(particleHash, d_particleHash, arruint2size, cudaMemcpyDeviceToHost);
  for(int k=0; k<10; k++){
  cout<<particleHash[k].x<<"  "<<particleHash[k].y<<endl;
  }

  //radixsort - the particleHash array 
  //run findcellstart on the global array - Ista and Iend - check collision
  //findcellstart<<<NUM/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(particleHash);
  //create the neighbour arrays for each particle
  //createNeighbourArrays<<<NUM/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(cellstart, cellend);

  cudaFree(d_particleHash);
  cudaFree(d_cellStart);
  cudaFree(d_cellEnd);
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);
  cudaFree(d_Xmin);
  cudaFree(d_Xmax);
  cudaFree(d_Ymin);
  cudaFree(d_Ymax);
  cudaFree(d_Zmin);
  cudaFree(d_Zmax);
  cudaFree(d_re);
  cudaFree(d_NUM);

}
/*

void NEIGHBOR()
{

	// ------------------PARAMETERS DEFENTION -------------------------------------
	int ncx = int((Xmax - Xmin) / (re + DELTA)) + 1;     // Number of cells in x direction
	int ncy = int((Ymax - Ymin) / (re + DELTA)) + 1;     // Number of cells in y direction
	int ncz = int((Zmax - Zmin) / (re + DELTA)) + 1;     // Number of cells in z direction

	int tnc = ncx*ncy*ncz;							   // Total number of cells   
	int m, k, kmax, Cnum;


	int *Ista, *Iend, *nc, *icell, *jcell, *kcell;
	int *ip;                             // I is sorted number of ip[I] th paricle
	Ista = new int[tnc + 1]; //this points to the index of the first element in a cell in the array ip
	Iend = new int[tnc + 1]; //index of the last element in a cell in the array ip
	nc = new int[tnc + 1];
	icell = new int[TP + 1];
	jcell = new int[TP + 1];
	kcell = new int[TP + 1];
	ip = new int[TP + 1]; //this is the main array that we are looking for, it is sorted 
  // according to cell numbers and it contains particle indices 



	//----------------- ALLOCATING PRTICLES IN CELLS --------------------------


	for (k = 1; k <= tnc; k++) //cell loop 
	{
		Ista[k] = 1;
		Iend[k] = 0;
		nc[k] = 0;
	}

	for (k = 1; k <= NUM; k++) //particle loop
	{
		icell[k] = int((x[k] - Xmin) / (re + DELTA)) + 1;
		jcell[k] = int((y[k] - Ymin) / (re + DELTA)) + 1;
		kcell[k] = int((z[k] - Zmin) / (re + DELTA)) + 1;

		Cnum = icell[k] + (jcell[k] - 1)*ncx + (kcell[k] - 1)*ncx*ncy;     // Cell number in which particle k located

		nc[Cnum]++;						            // Number of particle in cell Cnum
		Iend[Cnum]++;						        // Number of particle in cell Cnum 

		for (m = Iend[tnc]; m >= Iend[Cnum]; m--)
		{
			if (m>0) ip[m + 1] = ip[m];
		} //this block is there to create space at the end as and when new particles are added

		for (m = Cnum + 1; m <= tnc; m++)
		{
			Ista[m]++;
			Iend[m]++;
		}

		ip[Iend[Cnum]] = k;
	}


	//--------------- FINDIND NEIGHBORS ----------------------------------
	int JJ;
	for (I = 1; I <= NUM; I++)
	{
		k = 2;
		int row, colu, elev, m1, m2, m3, m4, m5, m6;
		if (icell[I] == 1)m1 = 0; else m1 = -1;
		if (icell[I] == ncx)m2 = 0; else m2 = +1;
		if (jcell[I] == 1)m3 = 0; else m3 = -1;
		if (jcell[I] == ncy)m4 = 0; else m4 = +1;
		if (kcell[I] == 1)m5 = 0; else m5 = -1;
		if (kcell[I] == ncz)m6 = 0; else m6 = +1;

		for (row = m1; row <= m2; row++) //could be -1 to 1 , the triple loop is basically there to find all the 9 cells around that particle, including the one in which it itself is
		{
			for (colu = m3; colu <= m4; colu++) 
			{
				for (elev = m5; elev <= m6; elev++)
				{

					Cnum = icell[I] + row + (jcell[I] - 1 + colu)*ncx + (kcell[I] - 1 + elev)*ncx*ncy;

					for (JJ = Ista[Cnum]; JJ <= Iend[Cnum]; JJ++)
					{
						J = ip[JJ]; //J is tha ACTUAL particle index 
						neighb[I][k] = J;
						k++;
					}
				}
			}
		}
		kmax = k - 2;
		neighb[I][1] = kmax; //this is the total number of neighbours, which is stored at the beginning 
		//if( neighb[I][1]>1098 ||neighb[I][1]*0!=0) printf("ERROR, the neighbors of particles %d is %d", I, neighb[I][1]);
	}
	//--------------------Clearing dynamic arrays ----------------------------

	delete[]Ista;
	delete[]Iend;
	delete[]nc;
	delete[]icell;
	delete[]jcell;
	delete[]kcell;
	delete[]ip;
	Ista = NULL; Iend = NULL; nc = NULL; icell = NULL; jcell = NULL; kcell = NULL, ip = NULL;
}

*/


int main(){

	neighbour_cuda();
	return 0;
}