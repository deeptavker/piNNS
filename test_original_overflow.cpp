#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <time.h>

using namespace std;



double *x, *y, *z;
double Xmax=9, Xmin=0;
double Ymax=9, Ymin=0;
double Zmax=9, Zmin=0;
double re=0.072, DELTA=0;
int NUM=265000;
int MAX_NEIGHB=1500;

void create_particles(int NUM){
    x = (double *)malloc(sizeof(double)*NUM);
    y = (double *)malloc(sizeof(double)*NUM);
    z = (double *)malloc(sizeof(double)*NUM);


    srand((unsigned)time(0)); 
    double lowest=0, highest=8; 
    double range=(highest-lowest)+1; 
    for(int index=0; index<NUM; index++){ 
        x[index] = lowest+double(range*rand()/(RAND_MAX + 1.0)); 
        y[index] = lowest+double(range*rand()/(RAND_MAX + 1.0)); 
        z[index] = lowest+double(range*rand()/(RAND_MAX + 1.0)); 
    } 
}

int **neighb;

void NEIGHBOUR_serial(){

  cout<<endl<<"Time study for NEIGHBOUR_serial()"<<endl;

  // ------------------PARAMETERS DEFENTION -------------------------------------
  int ncx = int((Xmax - Xmin) / (re + DELTA)) + 1;     // Number of cells in x direction
  int ncy = int((Ymax - Ymin) / (re + DELTA)) + 1;     // Number of cells in y direction
  int ncz = int((Zmax - Zmin) / (re + DELTA)) + 1;     // Number of cells in z direction

  int tnc = ncx*ncy*ncz;                 // Total number of cells   
  int m, k, kmax, Cnum;

  neighb = new int*[NUM+1];
  for(int i=0; i<NUM+1; i++){
    neighb[i] = new int[MAX_NEIGHB + 2];
  }


  int *Ista, *Iend, *nc, *icell, *jcell, *kcell;
  int *ip;                             // I is sorted number of ip[I] th paricle
  Ista = new int[tnc + 1]; //this points to the index of the first element in a cell in the array ip
  Iend = new int[tnc + 1]; //index of the last element in a cell in the array ip
  nc = new int[tnc + 1];
  icell = new int[NUM + 1];
  jcell = new int[NUM + 1];
  kcell = new int[NUM + 1];
  ip = new int[NUM + 1]; //this is the main array that we are looking for, it is sorted 
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
    icell[k] = int((x[k-1] - Xmin) / (re + DELTA)) + 1;
    jcell[k] = int((y[k-1] - Ymin) / (re + DELTA)) + 1;
    kcell[k] = int((z[k-1] - Zmin) / (re + DELTA)) + 1;

    Cnum = icell[k] + (jcell[k] - 1)*ncx + (kcell[k] - 1)*ncx*ncy;     // Cell number in which particle k located

    nc[Cnum]++;                       // Number of particle in cell Cnum
    Iend[Cnum]++;                   // Number of particle in cell Cnum 

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
  int JJ, J;
  for (int I = 1; I <= NUM; I++)
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

int main(){
	create_particles(NUM);
	NEIGHBOUR_serial();
}