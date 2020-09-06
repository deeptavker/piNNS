/*
Copyright (C) 2018-2020 Deep Tavker
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef _NEIGHBOUR_KERNELS_
#define _NEIGHBOUR_KERNELS_ 

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

  int cellarrsize =   sizeof(int);
  icell = (int *)malloc(cellarrsize);
  jcell = (int *)malloc(cellarrsize);
  kcell = (int *)malloc(cellarrsize);
  
  *icell = int((d_x[k] - *d_Xmin) / (*d_re + *d_DELTA)) + 1;
  *jcell = int((d_y[k] - *d_Ymin) / (*d_re + *d_DELTA)) + 1;
  *kcell = int((d_z[k] - *d_Zmin) / (*d_re + *d_DELTA)) + 1;

  int cellNum = *icell + (*jcell - 1)* *ncx + (*kcell - 1)* *ncx * *ncy;

  d_particleHash[k] = cellNum;
  d_particleid[k] = k;

}

__global__ void findCellStart(int *particleHash, int *cellStart, int *cellEnd, int *NUM){

  int k = threadIdx.x + blockIdx.x * blockDim.x; // here index value is equal to the cell number which starts with 1 
  if (particleHash[k] != particleHash[k+1] and k!= *NUM - 1){
    cellEnd[particleHash[k] - 1] = k;
    cellStart[particleHash[k+1] - 1] = k+1;
  }
  if(k == *NUM - 1){
    cellEnd[particleHash[k] - 1] = k;
  }
              
}

__global__ void createNeighbourArraysCUDA(int *d_neighb, int *cellStart, int *cellEnd, int *particleHash, int *particleid, int *ncx, int *ncy, int *ncz, int *d_max_neighb, int *test){

  int index = threadIdx.x + blockIdx.x * blockDim.x; 
  int pid, icell, jcell, kcell, cellNum;

  cellNum = particleHash[index]; 
  pid = particleid[index];
  
  int neighb_index = pid * (*d_max_neighb + 1);

  kcell = (cellNum - 1)/((*ncx) * (*ncy)) + 1;
  jcell = ((cellNum - 1) - ((kcell - 1)* (*ncx) * (*ncy)))/ *ncx + 1;
  icell = cellNum - 1 - *ncx * (jcell - 1) - (*ncx * *ncy)*(kcell - 1) + 1;

  int Cnum, J;
  int curr_neighb_num = 0;
  
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

        if (cellEnd[Cnum - 1] != -1){

        for (int JJ = cellStart[Cnum -1]; JJ <= cellEnd[Cnum - 1]; JJ++)
        {
          J = particleid[JJ];
          curr_neighb_num++;
          d_neighb[neighb_index + curr_neighb_num] = J+1; //here the index is shifted by one unit to conform to the original MPS convention
          
        }
      }
      }
    }
  }
  
  
  d_neighb[neighb_index] = curr_neighb_num;
  test[index] = d_neighb[neighb_index];
 
}

__global__ void InitializeCellDetails(int *cellStart, int *cellEnd){
  int index = threadIdx.x + blockIdx.x * blockDim.x; 
  cellStart[index] = 0; cellEnd[index] = -1;
}

__global__ void Template(int *particleHash, int *particleid, int *cellStart, int *cellEnd, int *ncx, int *ncy, int *ncz, int *size_neighbours){
  int index = threadIdx.x + blockDim.x * blockIdx.x;
  int pid, icell, jcell, kcell, cellNum;
  int *neighbours;
  neighbours = (int *)malloc(*size_neighbours);
  cellNum = particleHash[index]; 
  pid = particleid[index];

  kcell = (cellNum - 1)/((*ncx) * (*ncy)) + 1;
  jcell = ((cellNum - 1) - ((kcell - 1)* (*ncx) * (*ncy)))/ *ncx + 1;
  icell = cellNum - 1 - *ncx * (jcell - 1) - (*ncx * *ncy)*(kcell - 1) + 1;

  int Cnum, J;
  int curr_neighb_num = 0;
  
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

        if (cellEnd[Cnum - 1] != -1){

        for (int JJ = cellStart[Cnum -1]; JJ <= cellEnd[Cnum - 1]; JJ++)
        {
          J = particleid[JJ];
          curr_neighb_num++;
          neighbours[curr_neighb_num] = J;
          
        }
      }
      }
    }
  }
  
  
  neighbours[0] = curr_neighb_num;
  //any further operations can be done using this neighbour array

}

#endif
