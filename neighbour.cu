#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <vector>
#include <thrust/sort.h>

using namespace std;

double x[10]={1,1.001,3,4,5,6,7,8,9,0};
double y[10]={1,1.001,3,4,5,6,7,8,9,0};
double z[10]={1,1.001,3,4,5,6,7,8,9,0};
double Xmax=9, Xmin=0;
double Ymax=9, Ymin=0;
double Zmax=9, Zmin=0;
double re=5, DELTA=0;
int NUM=10;
int **neighb;


//declarations of global variables such as x, y, z, re, DELTA, xmax, etc/ 

//CUDA kernel for calcHash
__global__ void calcHash(double *d_x, double *d_y, double *d_z, int *d_particleHash,\
	int *d_NUM, double *d_Xmax, double *d_Xmin, double *d_re, double *d_DELTA, double *d_Ymin, \
  double *d_Ymax, double *d_Zmax, double *d_Zmin, int *d_particleid, int *d_tnc, int *ncx, int *ncy,\
  int *ncz){

  *ncx = int((*d_Xmax - *d_Xmin) / (*d_re + *d_DELTA)) + 1;     // Number of cells in x direction
  *ncy = int((*d_Ymax - *d_Ymin) / (*d_re + *d_DELTA)) + 1;     // Number of cells in y direction
  *ncz = int((*d_Zmax - *d_Zmin) / (*d_re + *d_DELTA)) + 1;     // Number of cells in z direction
  *d_tnc = *ncx * *ncy * *ncz;

  int k = threadIdx.x + blockIdx.x * blockDim.x;
  int *icell, *jcell, *kcell;

  int cellarrsize =  *d_NUM * sizeof(int);
  icell = (int *)malloc(cellarrsize);
  jcell = (int *)malloc(cellarrsize);
  kcell = (int *)malloc(cellarrsize);
  
  icell[k] = int((d_x[k] - *d_Xmin) / (*d_re + *d_DELTA)) + 1;
  jcell[k] = int((d_y[k] - *d_Ymin) / (*d_re + *d_DELTA)) + 1;
  kcell[k] = int((d_z[k] - *d_Zmin) / (*d_re + *d_DELTA)) + 1;

  int cellNum = icell[k] + (jcell[k] - 1)* *ncx + (kcell[k] - 1)* *ncx * *ncy;

  d_particleHash[k] = cellNum;
  d_particleid[k] = k;

}

__global__ void findCellStart(int *particleHash, int *cellStart, int *cellEnd){

  int k = threadIdx.x + blockIdx.x * blockDim.x; // here index value is equal to the cell number which starts with 1 
  if (particleHash[k] != particleHash[k+1] and k!=9){
    cellEnd[particleHash[k] - 1] = k;
    cellStart[particleHash[k+1] - 1] = k+1;
  }
  if(k == 9){
    cellEnd[particleHash[k] - 1] = k;
  }
              
}

//CUDA kernel for neighb[][] creation

__global__ void createNeighbourArrays(int **neighb, int *cellStart, int *cellEnd, int *particleHash, int *particleid, int *ncx, int *ncy, int *ncz){

  int index = threadIdx.x + blockIdx.x * blockDim.x; 
  //neighb[particleid[k]][2, 3, 4, 5, 6, 7....] = ids of other particles in that cell 
  int pid, icell, jcell, kcell, cellNum;

  cellNum = particleHash[index]; pid = particleid[index];
  kcell = (cellNum - 1)/((*ncx) * (*ncy)) + 1;
  jcell = ((cellNum - 1) - ((kcell - 1)* (*ncx) * (*ncy)))/ *ncx + 1;
  icell = cellNum - 1 - *ncx * (jcell - 1) - (*ncx * *ncy)*(kcell - 1) + 1;

  int k = 2;
  int Cnum, J;

  int row, colu, elev, m1, m2, m3, m4, m5, m6;
  if (icell == 1)m1 = 0; else m1 = -1;
  if (icell == *ncx)m2 = 0; else m2 = +1;
  if (jcell == 1)m3 = 0; else m3 = -1;
  if (jcell == *ncy)m4 = 0; else m4 = +1;
  if (kcell == 1)m5 = 0; else m5 = -1;
  if (kcell == *ncz)m6 = 0; else m6 = +1;

  for (row = m1; row <= m2; row++)
  {
    for (colu = m3; colu <= m4; colu++) 
    {
      for (elev = m5; elev <= m6; elev++)
      {

        Cnum = icell + row + (jcell - 1 + colu)* *ncx + (kcell - 1 + elev)* *ncx* *ncy;

        for (int JJ = cellStart[Cnum]; JJ <= cellEnd[Cnum]; JJ++)
        {
          J = particleid[JJ]; //J is tha ACTUAL particle index 
          neighb[pid][k] = J;
          k++;
        }
      }
    }
  }
  int kmax = k - 2;
  neighb[pid][1] = kmax;

}

__global__ void InitializeCellDetails(int *cellStart, int *cellEnd){
  int index = threadIdx.x + blockIdx.x * blockDim.x; 
  cellStart[index] = 0; cellEnd[index] = 0;
}



void neighbour_cuda(){

  // ------------------ variable declarations and initializations ------------------------------

  int *d_cellEnd, *d_cellStart, *d_NUM, *cellStart, *cellEnd, *d_tnc, *tnc, *d_ncx, *d_ncy, *d_ncz;
  int *particleHash, *d_particleHash, *d_particleid, *particleid, **d_neighb;
  double *d_x, *d_y, *d_z, *d_Xmax, *d_Xmin, *d_Ymax, *d_Ymin, *d_Zmax, *d_Zmin, *d_re, *d_DELTA;

  int arrsizeint = NUM * sizeof(int);
  int sizeint = sizeof(int);
  int arrsizedouble = NUM * sizeof(double);
  int sizedouble = sizeof(double);
  int sizeneighb = NUM * sizeof(int*);
  
  size_t *d_pitch, *pitch;

  particleHash = (int *)malloc(arrsizeint);
  particleid = (int *)malloc(arrsizeint);
  tnc = (int *)malloc(sizeint);
  
  

  cudaMalloc((void **)&d_particleHash, arrsizeint);
  cudaMalloc((void **)&d_particleid, arrsizeint); 
  
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
  cudaMalloc((void **)&d_tnc, sizeint);
  cudaMalloc((void **)&d_ncx, sizeint);
  cudaMalloc((void **)&d_ncy, sizeint);
  cudaMalloc((void **)&d_ncz, sizeint);

  cudaMallocPitch(d_neighb, d_pitch, 1500*sizeof(int), NUM);

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

  // --------------- running the calcHash kernel ----------------------------------------

  calcHash<<<5,2>>>(d_x, d_y, d_z, d_particleHash, d_NUM, d_Xmax, d_Xmin, d_re, d_DELTA, d_Ymin, d_Ymax, d_Zmax, d_Zmin, d_particleid, d_tnc, d_ncx, d_ncy, d_ncz);

  // ---------------- sorting the particleHash array -----------------------------

  thrust::device_ptr<int> dev_Hash(d_particleHash);
  thrust::device_ptr<int> dev_id(d_particleid);
  //thrust::device_ptr<int> dev_NUM(d_NUM);
  thrust::sort_by_key(dev_Hash, dev_Hash + 10, dev_id); //need to generalise this 10

  cudaMemcpy(particleHash, d_particleHash, arrsizeint, cudaMemcpyDeviceToHost);
  cudaMemcpy(particleid, d_particleid, arrsizeint, cudaMemcpyDeviceToHost);


  for(int k=0; k<10; k++){
  cout<<particleHash[k]<<"  "<<particleid[k]<<endl;
  }

  // --------------------- finding cell start and cell end for each cell -----------------------------

  cudaMemcpy(tnc, d_tnc, sizeint, cudaMemcpyDeviceToHost);
  int cellarrsize = *tnc * sizeof(int);
  cellStart = (int *)malloc(cellarrsize);
  cellEnd = (int *)malloc(cellarrsize);
  cudaMalloc((void **)&d_cellStart, cellarrsize); //need to initialise this with all zeros
  cudaMalloc((void **)&d_cellEnd, cellarrsize); //need to initialise this with all zeros 

  InitializeCellDetails<<<*tnc,1>>>(d_cellStart, d_cellEnd);
  findCellStart<<<5,2>>>(d_particleHash, d_cellStart, d_cellEnd);

  cudaMemcpy(cellStart, d_cellStart, cellarrsize, cudaMemcpyDeviceToHost);
  cudaMemcpy(cellEnd, d_cellEnd, cellarrsize, cudaMemcpyDeviceToHost);

  for(int i=0; i< *tnc; i++){
    cout<<cellStart[i]<<" "<<cellEnd[i]<<endl;
  }




  // -------------------------- Creating neighbour arrays for each particle ------------------------------

  createNeighbourArrays<<<5,2>>>(d_neighb, d_cellStart, d_cellEnd, d_particleHash, d_particleid, d_ncx, d_ncy, d_ncz);

  neighb = (int **)malloc(sizeof(int *)*NUM);
  for(int i=0; i<NUM; i++){
    neighb[i] = (int *)malloc(1500*sizeof(int));
  }

  cudaMemcpy2D(neighb, *pitch, d_neighb, *d_pitch, 1500*sizeof(int), NUM, cudaMemcpyDeviceToHost);

  cout<<neighb[1][1];


  cudaFree(d_particleHash);
  cudaFree(d_particleid);
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
  cudaFree(d_tnc);
  cudaFree(d_ncx);
  cudaFree(d_ncy);
  cudaFree(d_ncz);
  cudaFree(d_neighb);

  free(particleHash);
  free(particleid);
  free(cellStart);
  free(cellEnd);
  //free(x);
  //free(y);
  //free(z);
  //free(Xmin);
  //free(Xmax);
  //free(Ymin);
  //free(Ymax);
  //free(Zmin);
  //free(Zmax);
  //free(re);
  //free(NUM);
  free(tnc);
}

int main(){

	neighbour_cuda();
	return 0;
}