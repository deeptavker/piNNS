


/*
_____________________________________________________
|													   |
|     MPARS (Mesh-free Particle CFD Simulator)        |
|   Copyright (C)2014 
|_____________________________________________________|



                                      ``...........``                                               
                              ``.-::///////////////////::-..``                                      
                          `.://///////////////////////////////::..`                                 
                       .-/+++++/////////////////////////////////::::-.`                             
                    `-/+++++//////////////////////////////////::::::::::-`                          
                   -/++++++////////////////////////////////::::::::::::::::-.                       
                 `/++++++///////:--....--:////////////////::::::::::::::::::::-`                    
                ./+++++/////:.`            .:////////////:::::::::::::::::::::::-.                  
               `/+++++///:.                 `/////////::::: Copyright :::::::::::-.`               
               :+++////-`                   `:///////:::::::: AHMAD ::::::::::::::-.`             
               /++///:`                     .//////::::::: SHAKIBAEINIA ::::::::::::::-.            
              `/////.                      `://///:::::::::::::::::::::::::::::::::::----`          
              `////.                      `:///:::::::::::::::::::::::::::::::::::::------.         
               ///.                      .://:::::::::::::::::::::::::::::::::::::---------.        
               ./-                     `-://::::::::::::::::::::::::::::::::::::------------.       
                .                    `-://::::::::::::::::::::::::::::-----------------------.      
                                   `-:///::::::::::::::::::::---..`````                  `````      
                                `.:://::::::::::::::::::--.``         .:/-. .://:. .-/:.  `..`      
                             `.-:////:::::::::::::--..``      `-+o+/.ohhyhyoyhyyhy+yhyyho/yyyyo.    
                         `.-::////:::::::::::--.``      ./sss/+hhyyhyyhyyyhshhyyhyohhyyhyhhyyyho    
                   ``.-::://////::::::::--.``     `-:/:-yhyyyhyhhyyho.////:`.:/::.`:+//:.:oooo+.    
           ```..--::///////////:::::-.``    `://:.ohhyhhyyyyyy/.:---` ..``. `.```. ..``.` ....`     
       `.-:::////////////:::::--.``   `.::.-yhyyhyyhyyyho.:::. ..```` .```.``.```` ..`````-`  .`    
           ````.........````    `---`-shhyhyhhyyyh/:///:`-```.``.`````://:-`.+++/. -///:` .---`     
     ./oo+:``---.`  ```  `://:`/yhyyyyhhyyhh///:/. ..``.``...``:ooo+-shhyhhohhyyhh+hhhhhyoyyyyy-    
    .yhyyyhsshhhyo:syyyo/yhyyhyyhyyyhy:+++/-`.```. .````.+ssso/hhyyhhoyhhyy/syyyyo:yhhyyoyhhyhh/    
    `shyyyyohhyyyhyhyyyhyhhyyyh//++++.`.`````.```.`:///:shhyyhyoyyys/ .::-`  ....  `-//::/sys+-     
     `-:--.`:++++-:sssss-.////- `.``.`..` `.`////.ohhyhhyoysyo- ````              `-:::--..-//`     
      -```. ..``.  `...` `.```. .`  ``.:::-.yhhyhhsyyhyy:                       ./:-.-::.`.-:y+`    
      ..``` -`  `` -`  .``.` `. .-----hhhhhysyyyys.`-:.`            ``..::-`   :/.`.sso+o+os+:-`    
     -+ooo/``....  ..... ./+//-+hhhhhhhhhyys`.--.            ```-::/:::--.-//``s.``.///::://:-`     
     hhyyyhoohhhho+yyyys/hhyyhhyhhyhho`-/:-                 o/:/-.--:/:.....y+`o-.```````````-s-    
     /syyyo-hhyyhyhhhyyhsshyyys..:+:-         `..::/:       h..../ys+/s-.../y+ `:///////:.````yo.   
       ```  -/oo/``/o++:` .--.` ```````      ++//:--/o-     h..../yo/:-..:os+. ``.:o/---:s```-yo`   
 .```````     ````.---  `/+//++++oo+++++-   +/::++---/s-    h.....--..../os:` :+:-../+/::.`.-so:    
.hyyyhyyy/` `oyyyssssd:`/y++++ssysy////+h/`o/::/yy+---:y:   h....:hs+:.....+/``+-.`````..-/os+-     
.hsoohysosy/yoosdsoo+ho./y++++sssoo//+ossoo/:::hssy:----s/  h....:h+ :/.....:s. -://+++ooo+:.       
.hsysho ysssooyyhsoo+ho.+s++++yssoosso+:/s:::::/:::/:----s/ h-...:h+  .+/++ooso-   `...``           
 ::::::  ++syyo/ hyssho.+s++++do.`     :o/::/ysoo+//s///+oh/s/oooss/     ...``                      
              `   `----``-:+++o/`      /:oooo+/      .-:---.` ``                                    


                                                            
         `:oydmNmdho:          hddddddddddddhyys+:`         
       `sNMMMMMMMMMMMNo`       NMMMMMMMMMMMMMMMMMMMh/       
      `dMMMMNyo+odMMMMMd`      NMMMMMdhhhhhhhmNMMMMMMd-     
      oNMMMm`     yMMMMMo      NMMMMM-         -sMMMMMN:    
         `..      sMMMMM+      NMMMMM-           -NMMMMN.   
               `-sMMMMMy       NMMMMM-            oMMMMMs   
             :MMMMMMNy-        NMMMMM-            .MMMMMm   
             sMMMMMMmy+`       NMMMMM-            `MMMMMM`  
             :---+dMMMMNo      NMMMMM-            .MMMMMN`  
                   sMMMMMo     NMMMMM-            :MMMMMd   
                   -MMMMMN     NMMMMM-            yMMMMM+   
     .hdmNN+       /MMMMMd     NMMMMM-           sMMMMMh    
      hMMMMN+`   `+NMMMMN:     NMMMMM/-------:/sNMMMMMy`    
      `yMMMMMMmdmMMMMMMm-      NMMMMMMMMMMMMMMMMMMMMm/      
        -yNMMMMMMMMMNh/        NMMMMMMMMMMMMMMMMNho-        
           ./+oso+/-           :::::::::::::--.`            
                                           


*/

#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <vector>
#define STCK_SIZE 500000000

#include <algorithm>
#include <vector>
#include <thrust/sort.h>
#include <chrono>
#define THREADS_PER_BLOCK 512


using namespace std;
using namespace std::chrono;



//_______________________ Global Variables Definition ____________________________________

int I, J, K, NUM, l;

double n0 = 0, counter[21], Rho, MEUi, MEUj, NEU;
double b[10];
double Guess[10];
double A[10][10];
double M[10][10], Minv[10][10];
double *x, *y, *z, *p, *u, *v, *w, *SFX, *SFY, *SFZ;
int *PTYPE;
double *xstar, *ystar, *zstar, *ustar, *vstar, *wstar, *pnew, *phat, *unew, *vnew, *wnew, *NEUt, *TURB1, *TURB2, *n, *nstar, *MEU, *C, *RHO, *RHOSMOOTH, *MEU_Y, *p_rheo, *p_rheo_new, *II, *Inertia;
double *Tau_xx, *Tau_yy, *Tau_zz, *Tau_xy, *Tau_xz, *Tau_yz;

int    **neighb;
int    p_count;
double lambda;                          //         MPS discretization coefficient
double Xmin;						    //         Minimum x of searching grid
double Ymin;						    //         Minimum y of searching grid
double Zmin;						    //         Minimum z of searching grid
double Xmax;						    //         Maximum x of searching grid
double Ymax;						    //         Maximum y of searching grid
double Zmax;						    //         Maximum z of searching grid
int    FP;							    //         Number of Fluid particles: in this model FP calculate in each time step
int    WP;							    //         Number of wall particles
int    GP;								//         Number of ghost particles
double DT = 0.0050;                       //         Time step size
double c0;
int MAX_NEIGHB = 1500;



//**************************************
//***  Assigning Model Parameters ******
//**************************************

//*************** basic model conditions  *****************

double DL = 0.03;						   //         Average particle distance (or particle size)
double re = 0.1;                         //		   Support area radius
double BETA = 0.93;                       //         Coefficient for determination of free surface
double relaxp = 0.5;						//         Relaxation factor for pressure correction
double relaxu = 1.0;						//         Relaxation factor for u-velocity correction
double relaxv = 1.0;						//         Relaxation factor for v-velocity correction
double COURANT = 1;                       //         Courant number
double correction = 0;                    //         Correction factor to modify problem caused by shortage of ghoast particles
double coll = 0.85*DL;                   //       Minimum particle distance to prevent collision
double CC = 0.50;                         //         Collision coefficient
double EXPANSION = 0.0000001;				//         Expansion coefficient
double MAXresi = 0.001;                   //         Maximum acceptable residual for pressure calculation
double Cs = 0.18;                         //         Smogorinsky Constant (For using in SPS-LES turbulence model
double DELTA = DL / 5.0;                    //         Expansion value of background grid cells size.
double Ncorrection = 0.99;                //         Correction of n0 value to improve the incompressibility and initial pressure fluctuation
int    SP = 0;                            //         Number of storage particles
int    TP = 264815;						//         Total number of particles
int    KTYPE = 6;                         //         Kernel type
int    DIM = 3;                           //         Dimension
int    WBC = 0;                           //         Type of wall B.C.  No-Slip:0 & -1, Slip:1
int    TURB = 0;                          //         TURB=0: NO turbulence model, TURB=1: SPS turbulence model
int    Fraction_method = 2;               //         Method of calculation of volume of fraction. 1: Linear dist across the interface, 2: smoothed value

//*************** basic flow conditions  *****************
double NEU1 = 0.000002;                   //         Kinematic Viscosity
double NEU2 = 0.000002;                    //         K inematic Viscosity
double Rho1 = 1000;						    //         Density of phase 1
double Rho2 = 2500;						//         Density of phase 2
double gx = 0.0;                          //         Gravity acceleration in x direction
double gy = 0.0;                          //         Gravity acceleration in y direction
double gz = -9.806;                          //         Gravity acceleration in z direction
double VMAX = 10;	 					//         To avoid jumping particles out of domain (If v>Vmax ---> v=0)

//*************** pressure and pressure gradient Calculation parameters *****************

int    Method = 3;						//         Fully incompressible MPS: Method=1 ; Fully incompressible M-MPS: Method=2 ; Weakly compressible: Method=3 .
int    SOLVER = 1;                        //         SOLVER=1: CGM, SOLVER=2: PCGM (Fully incompressible models)
double c01 = 22;                        //         Numerical sound speed fluid 1. (Weakly compressible model)
double c02 = 22;                        //         Numerical sound speed fluid 2. (Weakly compressible model)
int    KHcorrection = 0;                  //         Khayyer and Gotoh pressure correction(1=yes, 0=no)
int    IterMax = 100;                     //         Maximum number of iterations for pressure calculation in each time step (Fully incompressible models)
double PMAX = 400000;                       //         A limit for the value of calculated pressure to avoid program crashing
double PMIN = -2.0;                        //         Minimum pressure, to avoid a high negative pressures
double GAMA = 7.0;

//*************** inflow & outflow parameters   *****************
int    Inflow_type = 0;                   //         0: No inflow,  1: Horizontal inflow section, 2: Vertical
double Inflow_velo = 0.0;                 //        inflow velocity perpendicular to the boundary
double Inflow_start_x = 0.0;              //        x coordination of inflow start point
double Inflow_start_y = 0.0;              //        y coordination of inflow start point
double Inflow_length = 0.0;               //        Length of the inflow section
int    Inflow_par_type = 1;               //         1: phase 1 particles, 2: phase 2 particles
int    counter_inflow = 1;
double Outflow_location = 0.0;           //        Location of the outflow boundary (i.e., the x coordinates) it has been assumed that the outflow is vertical
double Threshold = 0.0;                  //        Particles whose distance to the inflow are less than the (Threshold*DL), are removed from the domain

//*************** Rheological parameters   *****************

int    Fluid2_type = 1;                   // Newtonian:0  , H-B fluid:1
double N = 1.0;							//flow behaviour (power law) index
double MEU0 = 0.03;						// consistency index
double PHI = 0.70;                           // friction angle (RAD)
double PHI_wall = 0.04;                           // friction angle (RAD)
double PHI_bed = 0.04;                           // friction angle (RAD)
double PHI_2 = 0.95;                            //  second friction angle Values are based on  Minatti &  Paris (2015)
double cohes = 0;							//cohesiveness coefficient
double visc_max = 20;                // maximum viscosity uses to avoid singularity
double dg = 0.013;                    // grain size
double I0 = 0.70;                     // I0 value in Meu9I0 rheology     Values are based on  Minatti &  Paris (2015)
double mm = 200;
int    stress_cal_method = 1;        // Method 1; viscosity is directly used in momentum equation. Method 2: first the stress tensor is calculated then it used in momentum equation
int visc_itr_num = 1;
double visc_error = 0.0, visc_ave = 0.0;
double yeild_stress;





//*************** Surface tension parameters   *****************
int    surface_method = 0;                 // method of calculation of surface tension, 0:No surface tension  1:Koshizuka's model 2: CFS model (concentration based)
double SIGMA = 0.072;                       //surface tention coefficient (N/m)
double sharpness = 0.125;                  // interface sharpness  :   surface tension is added to the particles with C[I]=0.5+-sharpness
int    Fluid_surface = 1;                  // The fluid for which the suface tension is to be applied


//*************** Time parameters  *****************
double  t, T = 4.20;                       //         Simulation time (sec)
double  DT_MAX = 0.005;                   //         Maximum size of time interval allowed
double  output_time1 = 0.05;              //         Output times
double  output_time2 = 0.1;
double	output_time3 = 0.15;
double  output_time4 = 0.2;
double  output_time5 = 0.25;
double	output_time6 = 0.3;
double	output_time7 = 0.35;
double	output_time8 = 0.4;
double	output_time9 = 0.45;
double	output_time10 = 0.5;
double  output_time11 = 0.55;
double  output_time12 = 0.6;
double	output_time13 = 0.65;
double	output_time14 = 0.7;
double  output_time15 = 0.75;
double  output_time16 = 0.8;
double  output_time17 = 0.825;
double	output_time18 = 0.85;
double  output_time19 = 0.875;
double  output_time20 = 0.9;
double	output_time21 = 0.925;
double	output_time22 = 0.95;
double	output_time23 = 0.975;
double	output_time24 = 1;
double	output_time25 = 1.025;
double  output_time26 = 1.05;
double  output_time27 = 1.075;
double  output_time28 = 1.1;
double  output_time29 = 1.125;
double	output_time30 = 1.15;
double  output_time31 = 1.175;
double  output_time32 = 1.2;
double  output_time33 = 1.225;
double	output_time34 = 1.25;
double  output_time35 = 1.275;
double  output_time36 = 1.3;
double  output_time37 = 1.325;
double	output_time38 = 1.35;
double  output_time39 = 1.375;
double  output_time40 = 1.4;


double tec_out_intrval = 0.05;            // time interval for tecplot output




//_________________________________________________________________________________________________
// ######   ##    ##  ######   #######   ##    ## ######### ########  ###   ##  #######  ######
// ######   ##    ##  ######   #######   ##    ## ######### ########  ###   ##  #######  ######
//##    ##  ##    ##  ##   ##  ##    ##  ##    ##    ##        ##     ###   ##  ##      ##    ##
//##    ##  ##    ##  ##   ##  ##    ##  ##    ##    ##        ##     ###   ##  ##      ##    ##
//##        ##    ##  ##   ##  ##    ##  ##    ##    ##        ##     ## #  ##  ##      ##
// #####    ##    ##  ######   ##   ##   ##    ##    ##        ##     ## ## ##  ##       #####
// #####    ##    ##  ######   ##   ##   ##    ##    ##        ##     ## ## ##  #######  #####
//   #####  ##    ##  ##   ##  ######    ##    ##    ##        ##     ##  # ##  #######    #####
//   #####  ##    ##  ##   ##  ######    ##    ##    ##        ##     ##  # ##  ##          #####
//       ## ##    ##  ##    ## ##   ##   ##    ##    ##        ##     ##   ###  ##              ##
// ##    ## ##    ##  ##    ## ##    ##  ##    ##    ##        ##     ##    ##  ##        ##    ##
// ##    ## ##    ##  ##    ## ##    ##  ##    ##    ##        ##     ##    ##  ##        ##    ##
//  ######   ######   #######  ##    ##   ######     ##     ########  ##    ##  #######   ######
//  ######   ######   #######  ##    ##   ######     ##     ########  ##    ##  #######   ######


//===========================================================================================
//=====================     Kernel function       ===========================================
//===========================================================================================

//----------- This subroutine generate weight function or kernel function  --------- ----------
//-- Input variables are particle distance (R), kernel type, Space imension and support radius-----------
//            ________________________________________________
//           |   No. |    Kernel type                         |
//            ================================================
//           |...1...|Second order polynomial function        |
//           |...2...|Proportioal function (Kushizuka, 1998)  |
//           |...3...|Cubic splin func.                       |
//           |...4...|Quartic function                        |
//           |...5...|Quartic splin func.                     |
//           |___6___|________________________________________|

double W(double R, int KTYPE, int dim)
{
	double w;
	double q = R / re;
	//------------------------------------
	if (KTYPE == 1)                      // Second order polynomialfunction (Koshizuka and Oka 1996)
	{
		if (q<0.5)         w = 2.0 - 4.0*pow(q, 2);
		else if (q <= 1.0)   w = (2 * q - 2)*(2 * q - 2);
		else               w = 0;
	}
	//------------------------------------
	if (KTYPE == 2)                    // Rational function (Koshizuka and Oka 1998)
	{
		if (q <= 1.0)        w = (1 / q - 1);
		else               w = 0;
	}
	//------------------------------------
	if (KTYPE == 3)
	{
		double C;
		if (dim == 1) C = 0.6666;
		if (dim == 2) C = 1.43*3.14;
		if (dim == 3) C = 1 / 3.14;
		if (R<re)                     w = C / pow(re, dim)*(1 - 1.5*pow(q, 2.0) + 0.75*pow(q, 3.0));
		else if (R >= re && R<(2.0*re)) w = C / pow(re, dim)*(0.25*pow(2.0 - q, 3.0));
		else                          w = 0;
	}
	//------------------------------------

	if (KTYPE == 5)
	{

		if (q <= 1.0)                w = 1.5*log(1 / (q + 0.000001));
		else                       w = 0;
	}
	//--------------------------------

	if (KTYPE == 6)           // #rd oerder polynomial function (Shakibaeinia and Jin 2010)
	{

		if (q <= 1.0)                w = pow((1 - q), 3);
		else                       w = 0;

	}


	return (w);
}
//______________________________  END OF SUBROUTINE (W)  ____________________________________________


//=================================================================================================
//=====================     XIJ Calculation    ====================================================
//=================================================================================================
double DX(int i, int j)
{
	return (x[j] - x[i]);
}
//______________________________  END OF SUBROUTINE (DX)  __________________________________________

//================================================================================================
//=====================     XIJ* Calculation    ==================================================
//================================================================================================
double DXSTAR(int i, int j)
{
	return (xstar[j] - xstar[i]);
}
//______________________________  END OF SUBROUTINE (DXSTAR)  _______________________________________

//=================================================================================================
//=====================     YIJ Calculation     ===================================================
//=================================================================================================
double DY(int i, int j)
{
	return (y[j] - y[i]);
}
//______________________________  END OF SUBROUTINE (DY)  ___________________________________________


//=================================================================================================
//=====================     YIJ* Calculation     ==================================================
//=================================================================================================
double DYSTAR(int i, int j)
{
	return (ystar[j] - ystar[i]);
}
//______________________________  END OF SUBROUTINE (DYSTAR)  _______________________________________



//=================================================================================================
//=====================     ZIJ Calculation     ===================================================
//=================================================================================================
double DZ(int i, int j)
{
	return (z[j] - z[i]);
}
//______________________________  END OF SUBROUTINE (DZ)  ___________________________________________


//=================================================================================================
//=====================     ZIJ* Calculation     ==================================================
//=================================================================================================
double DZSTAR(int i, int j)
{
	return (zstar[j] - zstar[i]);
}
//______________________________  END OF SUBROUTINE (DZSTAR)  _______________________________________


//=================================================================================================
//=====================     Distance Calculation     ==============================================
//=================================================================================================
double DIST(int i, int j)
{
	double R;
	R = sqrt(pow(DX(i, j), 2.0) + pow(DY(i, j), 2.0) + pow(DZ(i, j), 2.0));
	return (R);
}
//______________________________  END OF SUBROUTINE (DIST)  _________________________________________


//=================================================================================================
//=====================     Distance Calculation     ==============================================
//=================================================================================================
double DISTSTAR(int i, int j)
{
	double R;
	R = sqrt(pow(DXSTAR(i, j), 2.0) + pow(DYSTAR(i, j), 2.0) + pow(DZSTAR(i, j), 2.0));
	return (R);
}
//______________________________  END OF SUBROUTINE (DISTSTAR)  _____________________________________



//=================================================================================================
//=====================     Particle Number Density Calc.  ========================================
//=================================================================================================
double PNUM(int I, int ktype, int dim, int num)
{
	double sum = 0.0;
	double d = 0;
	for (l = 2; l <= neighb[I][1]; l++)
	{
		J = neighb[I][l];
		d = DIST(I, J);
		if (I != J) sum = sum + W(d, ktype, dim);
	}
	if (KTYPE != 2)sum = sum*Ncorrection;
	return(sum);
}
//______________________________  END OF SUBROUTINE PNUM  ___________________________________________


//=================================================================================================
//=====================     Particle Number Density Calc.  ========================================
//=================================================================================================
double PNUMSTAR(int I, int ktype, int dim, int num)
{
	double sum = 0.0;
	double d = 0;
	for (l = 2; l <= neighb[I][1]; l++)
	{
		J = neighb[I][l];
		d = DISTSTAR(I, J);
		if (I != J) sum = sum + W(d, ktype, dim);
	}

	return(sum);
}
//______________________________  END OF SUBROUTINE PNUMSTAR  _______________________________________



//=================================================================================================
//==============   Conjugate Gradient Method For Linear Eq.s   ====================================
//=================================================================================================
//   ______________________________________________________________________________________
//  | This SUBROUTINE solve linear symmetric posetive definite system of equation (Ax=b)   |
//  | using conjugate gradient method. This Subroutine only change initial guess matrix      |
//  | to answer matrix and does not return any value.                                      |
//  |   Notice: A,b and Guess matrixes are global and should defined out of this function. |
//  |______________________________________________________________________________________|

//     NOTE: This subroutine has not been update since the initial code development in 2007 and 2008, Therefore it may need some modification to be used in the present code

double CGM(int size)
{
	double r1[13000], r2[13000], x1[13000], x2[13000], p1[13000], p2[13000];
	double alpha, beta, sum, sum1, sum2, redidual = 10;
	int k, m;

	//---------------------------- Initiation --------------------------
	for (k = 1; k <= size; k++)
	{
		x1[k] = Guess[k];        //replasing guess of initial answer matrix
	}

	for (k = 1; k <= size; k++)
	{
		sum = 0;
		for (m = 1; m <= size; m++)
		{
			sum = sum + A[k][m] * x1[m];  // calculation of initial errors and p1
		}
		r1[k] = b[k] - sum;
		p1[k] = r1[k];
	}
	//------------------------- repeatation loop -------------------------
	for (int repeat = 1; repeat<IterMax; repeat++)
	{
		sum1 = 0;
		sum2 = 0;
		for (k = 1; k <= size; k++)              //Calculation of alpha
		{
			sum = 0;
			for (m = 1; m <= size; m++)
			{
				sum = sum + A[k][m] * p1[m];
			}
			sum1 = sum1 + sum*p1[k];
			sum2 = sum2 + r1[k] * r1[k];
		}
		alpha = sum2 / sum1;

		for (k = 1; k <= size; k++)         //Calculation of new answers & errors
		{

			sum = 0;
			for (m = 1; m <= size; m++)
			{
				sum = sum + A[k][m] * p1[m];
			}
			r2[k] = r1[k] - alpha*sum;
			x2[k] = x1[k] + alpha*p1[k];
		}

		sum1 = 0;
		sum2 = 0;
		for (k = 1; k <= size; k++)         //Calculation of Beta
		{
			sum1 = sum1 + r1[k] * r1[k];
			sum2 = sum2 + r2[k] * r2[k];
		}
		beta = sum2 / sum1;

		sum = 0;
		for (k = 1; k <= size; k++)        //replacement of values for new repeat
		{
			p2[k] = r2[k] + beta*p1[k];
			p1[k] = p2[k];
			r1[k] = r2[k];
			x1[k] = x2[k];
			sum = sum + fabs(r2[k]);
		}
		redidual = sum / size;                      //calculation of average error to stop the loop
		if (redidual<MAXresi)repeat = IterMax;    // End of calculation if residual is less than maximum acceptable residual
	}

	for (k = 1; k <= size; k++)        //replacement of values for new reoeat
	{
		Guess[k] = x1[k];
	}
	return (0);
}
//______________________________  END OF SUBROUTINE (CGM) ___________________________________________



//=================================================================================================
//=================  preparing MATRIXES (A,b and guess for x) for CGM  ============================
//=================================================================================================

//     NOTE: This subroutine has not been update since the initial code development in 2007 and 2008, Therefore it may need some modification to be used in the present code
double MATRIX()
{
	int K = 0;
	double sum = 0, sum1 = 0, d = 0;
	//----- preparing matrixes (A,b and guess for x) for CGM by Original MPS-------
	if (Method == 1)
	{
		for (I = 1; I <= (FP + WP); I++)
		{
			sum1 = 0;
			for (J = 1; J <= NUM; J++)
			{
				d = DIST(I, J);
				if (I != J && d>coll / 10 && d<2 * re) sum1 = sum1 + W(d, KTYPE, 2);
			}


			if (PTYPE[I] == 0) correction = 1.0; else correction = 1.0;
			if (nstar[I]<BETA*n0)
			{
				b[I] = 0.00;
			}
			else
			{
				sum = 0;
				for (J = 1 + FP + WP; J <= NUM; J++)
				{
					d = DIST(I, J);
					if (d<(2.0*re) && d>coll / 10)sum = sum + pnew[J] * W(d, KTYPE, 2);
				}
				b[I] = correction*(-Rho*relaxp*lambda / (4 * pow(DT, 2)))*(nstar[I] - n0) - sum;

			}
			Guess[I] = p[I];
			for (J = 1; J <= (FP + WP); J++)
			{
				d = DIST(I, J);
				if (I == J)    { A[I][J] = -n[I] * correction; }

				//		if (I==J)A[I][J]=-sum1-EXPANSION*Rho*lambda*n0/(4*pow(DT,2));

				else if (d<(2.0*re) && d>coll / 10) A[I][J] = W(d, KTYPE, 2);
				else A[I][J] = 0;
			}
		}
	}



	//----- preparing matrixes (A,b and guess for x) for CGM by Modified MPS-------
	if (Method == 2)
	{
		for (I = 1; I <= NUM - GP; I++)
		{
			if (nstar[I]<BETA*n0)
			{
				b[I] = 0.00;
			}
			else
			{
				sum = 0;
				for (J = 3088; J <= NUM; J++)
				{
					sum = sum + pnew[J] * W(DIST(I, J), KTYPE, 2);
				}
				b[I] = (-Rho / (4 * pow(DT, 2)))*(nstar[I] - n0) - sum;
			}

			Guess[I] = p[I];
			for (J = 1; J <= NUM - GP; J++)
			{
				if (I == J)
				{
					sum = 0;
					for (K = 1; K <= NUM; K++)
					{
						if (I != K && DIST(I, K)<(2 * re) && DIST(I, K)>0.001)
							sum = sum + ((1.0 / DIST(I, K))*W(DIST(I, K), KTYPE, 2));
					}
					A[I][J] = -sum;
				}
				else if (DIST(I, J)<(2.0*re)) A[I][J] = W(DIST(I, J), KTYPE, 2) / DIST(I, J);
				else A[I][J] = 0;

			}
		}
	}
	//-------------------------------------------------------------------------------


	if (Method == 3)
	{

	}
	//-------------------------------------------------------------------------------
	return(0);
}
//______________________________  END OF Subroutine (MATRIX) ________________________________________




//==================================================================================================
//========================        Boundary condition (BC)      =====================================
//==================================================================================================

//----This subroutine assigne the values to boundariy particles (i.e. ghost and wall particles)

double BC(int slip)
{
	double MINIMUM = 100;
	int k1 = TP + 1, k2 = TP + 1;

	//-----------------------------------------------------

	if (PTYPE[I] == -2)                 // bottom
	{
		for (l = 2; l <= neighb[I][1]; l++)
		{
			J = neighb[I][l];
			/*	if ((PTYPE[J]>0) && z[J] <= 0.0 + DL)
			{
			MINIMUM = 100;
			if (fabs(DX(I, J))<MINIMUM)
			{
			k1 = J;
			MINIMUM = fabs(DX(I, J));
			}
			}*/

			if (PTYPE[J] == 0 && x[J] <= x[I] + DL / 2.0 && x[J] >= x[I] - DL / 2.0 && z[J]<0.0 + DL)
			{
				k2 = J;
			}

		}
		w[I] = 0.0, wstar[I] = 0.0, wnew[I] = 0.0;
		u[I] = 0.0, ustar[I] = 0.0, unew[I] = 0.0;
		v[I] = 0.0, vstar[I] = 0.0, vnew[I] = 0.0;

		p[I] = p[k2], pnew[I] = pnew[k2];
	}
	//-----------------------------------------------------
	//-----------------------------------------------------

	else if (PTYPE[I] == -1)                 // right walls
	{
		for (l = 2; l <= neighb[I][1]; l++)
		{
			J = neighb[I][l];
			/*	if ((PTYPE[J]>0) && x[J] >= 2.52 - DL)
			{
			MINIMUM = 100;
			if (fabs(DY(I, J))<MINIMUM)
			{
			k1 = J;
			MINIMUM = fabs(DY(I, J));
			}
			}*/
			if (PTYPE[J] == 0 && z[J] <= z[I] + DL / 2 && z[J] >= z[I] - DL / 2 && x[J] >= 2.52 - DL)
			{
				k2 = J;
			}

		}

		w[I] = 0.0, wstar[I] = 0.0, wnew[I] = 0.0;
		u[I] = 0.0, ustar[I] = 0.0, unew[I] = 0.0;
		v[I] = 0.0, vstar[I] = 0.0, vnew[I] = 0.0;

		p[I] = p[k2], pnew[I] = pnew[k2];
	}
	//-----------------------------------------------------

	/*else if (PTYPE[I]==-3  )                 // left walls
	{
	for (l=2;l<=neighb[I][1];l++)
	{
	J=neighb[I][l];
	if ((PTYPE[J]>0) && x[J]<=0.0+DL)
	{
	MINIMUM=100;
	if (fabs(DY(I,J))<MINIMUM)
	{
	k1=J;
	MINIMUM=fabs(DY(I,J));
	}
	}
	if (PTYPE[J]==0  && y[J]<=y[I]+DL/2 && y[J]>=y[I]-DL/2 && x[J]<=0.0+DL)
	{
	k2=J;
	}

	}

	v[I]=slip*v[k1],vstar[I]=slip*vstar[k1],vnew[I]=slip*vnew[k1];
	u[I]=0.0,ustar[I]=0.0,unew[I]=0.0;
	p[I]=p[k2],pnew[I]=pnew[k2];
	}*/
	//-----------------------------------------------------
	else if (PTYPE[I] == 0)       // wall particles
	{
		u[I] = 0, ustar[I] = 0, unew[I] = 0;
		v[I] = 0, vstar[I] = 0, vnew[I] = 0;

	}


	return(0.0);
}
//______________________________  END OF SUBROTINE (BC) ___________________________________________



//=================================================================================================
//===============  Finding particles in viciniy of given particle  ================================
//=================================================================================================
int NEIGHBOR1(int i, int num)
{
	double min = 100.00;
	int neighbor = 13001;
	for (int j = 1; j <= num; j++)
	{
		if (DIST(i, j)<min && PTYPE[j] == 0)
		{
			min = DIST(i, j);
			neighbor = j;
		}
	}
	return (neighbor);
}
//______________________________  END OF FUNCTION (NEIGHBOR1) _____________________________________

__global__ void calcHash(double *d_x, double *d_y, double *d_z, int *d_particleHash,\
  int *d_NUM, double *d_Xmax, double *d_Xmin, double *d_re, double *d_DELTA, double *d_Ymin, \
  double *d_Ymax, double *d_Zmax, double *d_Zmin, int *d_particleid, int *d_tnc, int *ncx, int *ncy,\
  int *ncz){

  int k = threadIdx.x + blockIdx.x * blockDim.x;
  if(k < *d_NUM){


  *ncx = int((*d_Xmax - *d_Xmin) / (*d_re + *d_DELTA)) + 1;     // Number of cells in x direction
  *ncy = int((*d_Ymax - *d_Ymin) / (*d_re + *d_DELTA)) + 1;     // Number of cells in y direction
  *ncz = int((*d_Zmax - *d_Zmin) / (*d_re + *d_DELTA)) + 1;     // Number of cells in z direction
  *d_tnc = *ncx * *ncy * *ncz;

  
  int *icell, *jcell, *kcell, *cellNum;

  int sizeint = sizeof(int);
  icell = (int *)malloc(sizeint);
  jcell = (int *)malloc(sizeint);
  kcell = (int *)malloc(sizeint);
  cellNum = (int *)malloc(sizeint);
  
  *icell = int((d_x[k] - *d_Xmin) / (*d_re + *d_DELTA)) + 1;
  *jcell = int((d_y[k] - *d_Ymin) / (*d_re + *d_DELTA)) + 1;
  *kcell = int((d_z[k] - *d_Zmin) / (*d_re + *d_DELTA)) + 1;

  *cellNum = *icell + (*jcell - 1)* *ncx + (*kcell - 1)* *ncx * *ncy;

  d_particleHash[k] = *cellNum;
  d_particleid[k] = k;

  
  free(icell);
  free(jcell);
  free(kcell);
  free(cellNum);
}

}

__global__ void findCellStart(int *particleHash, int *cellStart, int *cellEnd, int *NUM){

  int k = threadIdx.x + blockIdx.x * blockDim.x; // here index value is equal to the cell number which starts with 1 
  if(k < *NUM){
  if (particleHash[k] != particleHash[k+1] and k!= *NUM - 1){
    cellEnd[particleHash[k] - 1] = k;
    cellStart[particleHash[k+1] - 1] = k+1;
  }
  if(k == *NUM - 1){
    cellEnd[particleHash[k] - 1] = k;
  }
    }

  free(&k);            
}

__global__ void createNeighbourArraysCUDA(int *d_neighb, int *cellStart, int *cellEnd, int *particleHash, int *particleid, int *ncx, int *ncy, int *ncz, int *d_max_neighb,  int *d_NUM){

  int index = threadIdx.x + blockIdx.x * blockDim.x; 

  if(index < *d_NUM){
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
 }
}

__global__ void InitializeCellDetails(int *cellStart, int *cellEnd, int *d_tnc){
  int index = threadIdx.x + blockIdx.x * blockDim.x; 
  if(index < *d_tnc){
  cellStart[index] = 0; cellEnd[index] = -1;
}
free(&index);
}





void neighbour_cuda_1(){

  //cout<<endl<<"Time study for neighbour_cuda_1()"<<endl;

  // ------------------ variable declarations and initializations ------------------------------

  int *d_cellEnd, *d_cellStart, *d_NUM, *d_tnc, *tnc, *d_ncx, *d_ncy, *d_ncz, *d_max_neighb;
  int *d_particleHash, *d_particleid, *d_neighb, *h_neighb, *d_sizeof_neighbours;
  double *d_x, *d_y, *d_z, *d_Xmax, *d_Xmin, *d_Ymax, *d_Ymin, *d_Zmax, *d_Zmin, *d_re, *d_DELTA;

  int arrsizeint = NUM * sizeof(int);
  int sizeint = sizeof(int);
  int arrsizedouble = NUM * sizeof(double);
  int sizedouble = sizeof(double);
  int sizeneighb = NUM * (MAX_NEIGHB + 1) * sizeof(int);
  int sizeof_neighbours = (MAX_NEIGHB + 1) * sizeof(int);

  tnc = (int *)malloc(sizeint);
  h_neighb = (int *)malloc(sizeneighb);



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
  cudaMalloc((void **)&d_neighb, sizeneighb);
  cudaMalloc((void **)&d_max_neighb, sizeint);
  cudaMalloc((void **)&d_sizeof_neighbours, sizeof_neighbours);

  cudaMemcpy(d_x, x, arrsizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, arrsizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_z, z, arrsizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Xmin, &Xmin, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Xmax, &Xmax, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ymin, &Ymin, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Ymax, &Ymax, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Zmin, &Zmin, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_Zmax, &Zmax, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_re, &re, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_DELTA, &DELTA, sizedouble, cudaMemcpyHostToDevice);
  cudaMemcpy(d_NUM, &NUM, sizeint, cudaMemcpyHostToDevice);
  cudaMemcpy(d_max_neighb, &MAX_NEIGHB, sizeint, cudaMemcpyHostToDevice);
  cudaMemcpy(d_sizeof_neighbours, &sizeof_neighbours, sizeint, cudaMemcpyHostToDevice);


  

  // --------------- running the calcHash kernel ----------------------------------------
 
  calcHash<<<NUM/THREADS_PER_BLOCK + 1,THREADS_PER_BLOCK>>>(d_x, d_y, d_z, d_particleHash, d_NUM, d_Xmax, d_Xmin, d_re, d_DELTA, d_Ymin, d_Ymax, d_Zmax, d_Zmin, d_particleid, d_tnc, d_ncx, d_ncy, d_ncz);
 
  // ---------------- sorting the particleHash array -----------------------------

 
  thrust::device_ptr<int> dev_Hash(d_particleHash);
  thrust::device_ptr<int> dev_id(d_particleid);
  thrust::sort_by_key(dev_Hash, dev_Hash + NUM, dev_id); //need to generalise this 10
 
  
  // --------------------- finding cell start and cell end for each cell -----------------------------

  cudaMemcpy(tnc, d_tnc, sizeint, cudaMemcpyDeviceToHost);
  int cellarrsize = *tnc * sizeof(int);
  cudaMalloc((void **)&d_cellStart, cellarrsize); 
  cudaMalloc((void **)&d_cellEnd, cellarrsize); 

 
  InitializeCellDetails<<<*tnc/THREADS_PER_BLOCK + 1,THREADS_PER_BLOCK>>>(d_cellStart, d_cellEnd, d_tnc);
  findCellStart<<<NUM/THREADS_PER_BLOCK + 1,THREADS_PER_BLOCK>>>(d_particleHash, d_cellStart, d_cellEnd, d_NUM);
 
  
  // -------------------------- Creating neighbour arrays for each particle ------------------------------

 
  createNeighbourArraysCUDA<<<NUM/THREADS_PER_BLOCK + 1,THREADS_PER_BLOCK>>>(d_neighb, d_cellStart, d_cellEnd, d_particleHash, d_particleid, d_ncx, d_ncy, d_ncz, d_max_neighb, d_NUM);

  
  
 
  cudaMemcpy(h_neighb, d_neighb, sizeneighb, cudaMemcpyDeviceToHost);



  // ---------------------------- Populating neighb array ----------------------
  
  
  
  for(int j=0; j<NUM; j++){
    for(int i=0; i<h_neighb[j*(MAX_NEIGHB + 1)]; i++){
      neighb[j+1][i+2] = h_neighb[j*(MAX_NEIGHB + 1) + i + 1];
    }
    neighb[j+1][1] = h_neighb[j*(MAX_NEIGHB + 1)];
  }
  
 
  
  
  // -------------------------- Deallocating memory ---------------------------

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
  cudaFree(d_max_neighb);
  cudaFree(d_sizeof_neighbours);

  free(h_neighb);
  free(tnc);
}

//=================================================================================================
//===============  Finding particles in viciniy of given particle  ================================
//=================================================================================================
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

//______________________________  END OF SUBROTINE (NEIGHBOR) _____________________________________


//______________________________  END OF SUBROTINE (NEIGHBOR_CUDA) _____________________________________


//================================================================================================
//====================  Collision of Particles computation =======================================
//===========================================================================================
double COLLISION(int i, double MINdistance)
{
	double cc;
	double ug, vg, wg, um, vm, wm, ur, vr, wr, vabs, d;
	double m1;
	double m2;

	if (PTYPE[i] == 1)m1 = Rho1; else m1 = Rho2;


	for (l = 2; l <= neighb[I][1]; l++)
	{
		int j = neighb[I][l];
		if (PTYPE[j] == 1) m2 = Rho1; else m2 = Rho2;

		d = DISTSTAR(i, j);
		if (i != j && d<MINdistance)
		{
			cc = CC*(MINdistance - d) / MINdistance;
			//cc=CC;
			ug = (m1*ustar[i] + m2*ustar[j]) / (m1 + m2);
			vg = (m1*vstar[i] + m2*vstar[j]) / (m1 + m2);
			wg = (m1*wstar[i] + m2*wstar[j]) / (m1 + m2);

			ur = m1*(ustar[i] - ug);
			vr = m1*(vstar[i] - vg);
			wr = m1*(wstar[i] - wg);
			vabs = (ur*DXSTAR(i, j) + vr*DYSTAR(i, j) + wr*DZSTAR(i, j)) / d;

			um = (1.0 + cc)*vabs*DXSTAR(i, j) / d;
			vm = (1.0 + cc)*vabs*DYSTAR(i, j) / d;
			wm = (1.0 + cc)*vabs*DZSTAR(i, j) / d;
			if (vabs>0)
			{
				if (PTYPE[i]>0)
				{
					ustar[i] = ustar[i] - um / m1;
					vstar[i] = vstar[i] - vm / m1;
					wstar[i] = wstar[i] - wm / m1;
					xstar[i] = xstar[i] - DT*um / m1;
					ystar[i] = ystar[i] - DT*vm / m1;
					zstar[i] = zstar[i] - DT*wm / m1;
				}

				if (PTYPE[j]>0)
				{
					ustar[j] = ustar[j] + um / m2;
					vstar[j] = vstar[j] + vm / m2;
					wstar[j] = wstar[j] + wm / m2;
					xstar[j] = xstar[j] + DT*um / m2;
					ystar[j] = ystar[j] + DT*vm / m2;
					zstar[j] = zstar[j] + DT*wm / m2;


				}
			}
		}
	}
	return (0);
}
//______________________________  END OF SUBROTINE (COLLISION) _________________________________


//===========================================================================================
//====================  Collision of Particles computation ==================================
//===========================================================================================
double COLLISION2(int i, double MINdistance)
{
	double cc;
	double ug, vg, wg, um, vm, wm, ur, vr, wr, vabs, d;
	double m1;
	double m2;

	if (PTYPE[i] == 1)m1 = Rho1; else m1 = Rho2;

	for (l = 2; l <= neighb[I][1]; l++)
	{
		int j = neighb[I][l];
		if (PTYPE[j] == 1)m2 = Rho1; else m2 = Rho2;
		d = DIST(i, j);
		if (i != j && d<MINdistance)
		{
			cc = CC*(MINdistance - d) / MINdistance;
			//cc=CC;
			ug = (m1*unew[i] + m2*unew[j]) / (m1 + m2);
			vg = (m1*vnew[i] + m2*vnew[j]) / (m1 + m2);
			wg = (m1*wnew[i] + m2*wnew[j]) / (m1 + m2);


			ur = m1*(unew[i] - ug);
			vr = m1*(vnew[i] - vg);
			wr = m1*(wnew[i] - wg);

			vabs = (ur*DX(i, j) + vr*DY(i, j) + wr*DZ(i, j)) / d;

			um = (1.0 + cc)*vabs*DX(i, j) / d;
			vm = (1.0 + cc)*vabs*DY(i, j) / d;
			wm = (1.0 + cc)*vabs*DZ(i, j) / d;


			if (vabs>0)
			{
				if (PTYPE[i]>0)
				{
					unew[i] = unew[i] - um / m1;
					vnew[i] = vnew[i] - vm / m1;
					wnew[i] = wnew[i] - wm / m1;
					x[i] = x[i] - DT*um / m1;
					y[i] = y[i] - DT*vm / m1;
					z[i] = z[i] - DT*wm / m1;
				}

				if (PTYPE[j]>0)
				{
					unew[j] = unew[j] + um / m2;
					vnew[j] = vnew[j] + vm / m2;
					wnew[j] = wnew[j] + wm / m2;

					x[j] = x[j] + DT*um / m2;
					y[j] = y[j] + DT*vm / m2;
					z[j] = z[j] + DT*wm / m2;
				}
			}
		}
	}
	return (0);
}
//______________________________  END OF SUBROTINE (COLLISION) _________________________________


//===========================================================================================
//====================  Calculation of DT to satisfy Courant no. ============================
//===========================================================================================

void DTcalculation()
{
	double max = 0;

	//	for (I=GP+WP+1; I<=NUM; I++)
	//	{
	//		if (fabs(unew[I])>max) {max=fabs(unew[I]);}
	//		if (fabs(vnew[I])>max) {max=fabs(vnew[I]);}
	//	}
	if (c01>c02)c0 = c01;
	else c0 = c02;

	max = c0;
	double courant = COURANT*DL / max;
	if (DT>courant) DT = courant;
	else if (courant<DT_MAX)DT = courant;
	else DT = DT_MAX;

	return;
}
//______________________________  END OF SUBROUTIN(DTcalculation) ____________________________



//===========================================================================================
//================ Particle reflection after collision with walls============================
//======================== Not applicable in this model======================================

void REFLECTION(double ni, double nj, double nk)
{
	// V'=V-2(V.n)n     : n:unit nornal vector

	double vn = unew[I] * ni + vnew[I] * nj + wnew[I] * nk;                                        //V.n
	unew[I] = unew[I] - 2 * ni*vn;
	vnew[I] = vnew[I] - 2 * nj*vn;
	wnew[I] = wnew[I] - 2 * nk*vn;
	return;
}
//______________________________  END OF Subroutine (REFLECTION) ____________________________


//===========================================================================================
//================================  SPS LES Turbulence Model ================================
//===========================================================================================

void SPS()
{


	double  *S12, *S11, *S22, S, d;


	double Uxx, Uxy, Uyx, Uyy;
	double w;
	int i, j;

	S12 = new double[TP + 1];
	S11 = new double[TP + 1];
	S22 = new double[TP + 1];

	for (i = GP + 1; i <= NUM; i++)
	{
		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
		for (l = 2; l <= neighb[i][1]; l++)
		{
			j = neighb[i][l];
			d = DIST(i, j);
			if (i != j && d>coll / 10)
			{
				w = W(d, KTYPE, 2);
				sum1 = sum1 + (unew[j] - unew[i])*DX(i, j)*w / d / d;
				sum2 = sum2 + (vnew[j] - vnew[i])*DY(i, j)*w / d / d;
				sum3 = sum3 + (unew[j] - unew[i])*DY(i, j)*w / d / d;
				sum4 = sum4 + (vnew[j] - vnew[i])*DX(i, j)*w / d / d;
			}
		}
		Uxx = (2.0 / n0)*sum1;
		Uyy = (2.0 / n0)*sum2;
		Uxy = (2.0 / n0)*sum3;
		Uyx = (2.0 / n0)*sum4;

		S12[i] = 0.5*(Uxy + Uyx);
		S11[i] = 0.5*(Uxx + Uxx);
		S22[i] = 0.5*(Uyy + Uyy);


		// NEUt=pow(Cs*DL,2)*(pow((Uxx+Uxx)*Uxx,0.5)+pow((Uxy+Uyx)*Uxy,0.5)+pow((Uyx+Uxy)*Uyx,0.5)+pow((Uyy+Uyy)*Uyy,0.5);
		S = 2 * (S11[i] * S11[i] + S22[i] * S22[i] + 2 * S12[i] * S12[i]);
		if (S<0) S = 0;
		S = sqrt(S);

		NEUt[i] = pow(Cs*DL, 2)*S;

		if (NEUt[i] * 0 != 0)  NEUt[i] = 0;
		if (NEUt[i]>1)     NEUt[i] = 1;


	}
	for (i = GP + 1; i <= NUM; i++)
	{
		double sum1 = 0, sum2 = 0;
		for (l = 2; l <= neighb[i][1]; l++)
		{
			j = neighb[i][l];
			d = DIST(i, j);
			if (i != j && d>coll / 10)
			{
				w = W(d, KTYPE, 2);

				sum1 = sum1 + (NEUt[j] - NEUt[i])*DX(i, j)*w / d / d;    //d(neo)/dx
				sum2 = sum2 + (NEUt[j] - NEUt[i])*DY(i, j)*w / d / d;	   //d(neo)/dy
			}
		}
		TURB1[i] = (2.0 / n0)*(2 * S11[i] * sum1 + 2 * S12[i] * sum2);       //additional terms in turbulence model (Report #8)
		TURB2[i] = (2.0 / n0)*(2 * S12[i] * sum1 + 2 * S22[i] * sum2);
	}


	delete[]S11; delete[]S12; delete[]S22;
	S11 = NULL; S12 = NULL; S22 = NULL;

	return;

}
//______________________________  END OF Subroutine (SPS) ____________________________


//===========================================================================================
//============================= Fluid particles Re-arangement================================
//===========================================================================================

void SORT()
{
	FP = 0;
	for (I = GP + WP + 1; I <= TP; I++)
		//	for (I=TP;I>=GP+WP+1;I--)
	{

		if (PTYPE[I]>0)                             // Rearange fluid particles
		{
			FP++;
			x[FP + GP + WP] = x[I];
			y[FP + GP + WP] = y[I];
			p[FP + GP + WP] = p[I];
			u[FP + GP + WP] = u[I];
			v[FP + GP + WP] = v[I];
			C[FP + GP + WP] = C[I];
			PTYPE[FP + GP + WP] = PTYPE[I];
			xstar[FP + GP + WP] = xstar[I];
			ystar[FP + GP + WP] = ystar[I];
			ustar[FP + GP + WP] = ustar[I];
			vstar[FP + GP + WP] = vstar[I];
			pnew[FP + GP + WP] = pnew[I];
			ustar[FP + GP + WP] = ustar[I];
			ustar[FP + GP + WP] = ustar[I];
			unew[FP + GP + WP] = unew[I];
			vnew[FP + GP + WP] = vnew[I];
			NEUt[FP + GP + WP] = NEUt[I];
			TURB1[FP + GP + WP] = TURB1[I];
			TURB2[FP + GP + WP] = TURB2[I];
			nstar[FP + GP + WP] = nstar[I];
		}
	}

	for (I = GP + WP + FP + 1; I <= TP; I++)                   // Rearange null particles
	{

		x[I] = 0;
		y[I] = 0;
		p[I] = 0;
		u[I] = 0;
		v[I] = 0;
		C[I] = 0;
		PTYPE[I] = -100;
		xstar[I] = 0;
		ystar[I] = 0;
		ustar[I] = 0;
		vstar[I] = 0;
		pnew[I] = 0;
		ustar[I] = 0;
		ustar[I] = 0;
		unew[I] = 0;
		vnew[I] = 0;
		NEUt[I] = 0;
		TURB1[I] = 0;
		TURB2[I] = 0;
		nstar[I] = 0;
	}

	return;
}
//______________________________  END OF Subroutine (SORT) ____________________________

//===========================================================================================
//============ Add new particle to entrance base on inflow velocity ==========================
//===========================================================================================


void INFLOW()
{

	I = FP + WP + GP + 1;
	int mm = Inflow_length / DL;
	int m = 1;

	//************* Horizontal inflow *************
	if (Inflow_type == 1)
	{
		double vv = Inflow_velo;
		double tt = DL / fabs(vv);
		if (t >= tt*counter_inflow)
		{
			for (m = 1; m <= mm; m++)
			{

				double X = Inflow_start_x + m*DL;

				x[I] = X;
				y[I] = Inflow_start_y;
				p[I] = 0;
				pnew[I] = 0;
				u[I] = 0;
				unew[I] = 0;
				v[I] = vv;
				vnew[I] = vv;
				PTYPE[I] = Inflow_par_type;
				if (PTYPE[I] == 1) RHO[I] = Rho1;
				else             RHO[I] = Rho2;

				I++;
			}
			counter_inflow++;
		}
	}
	//************* Vertical inflow *************
	if (Inflow_type == 2)
	{
		double uu = Inflow_velo;
		double tt = DL / fabs(uu);
		if (t >= tt*counter_inflow)
		{
			for (m = 1; m <= mm; m++)
			{

				double Y = Inflow_start_y + m*DL;

				x[I] = Inflow_start_x;
				y[I] = Y;
				p[I] = 0;
				pnew[I] = 0;
				u[I] = uu;
				unew[I] = uu;
				v[I] = 0.0;
				vnew[I] = 0.0;
				PTYPE[I] = Inflow_par_type;
				if (PTYPE[I] == 1) RHO[I] = Rho1;
				else             RHO[I] = Rho2;
				I++;
			}
			counter_inflow++;
		}
	}


}
//______________________________  END OF Subroutine (INFLOW) ____________________________



//===========================================================================================
//============ AVOIDING PARTICLE TO PENETRATE  BOUNDARIES  ==========================
//===========================================================================================


void BOUNDARIE()
{
	//if (x[I]<-2.00 )															{unew[I]=+fabs(unew[I]);x[I]=-2.00;}           //left wall
	if (x[I]>+2.52 - DL / 2.0)												        { unew[I] = -fabs(unew[I]); x[I] = 2.52 - DL / 2.0; }           //right wall


	if (y[I]<-1.02 + DL / 2)														  { vnew[I] = fabs(vnew[I]); y[I] = -1.02 + DL / 2.0; }            // front wall
	if (y[I]>+1.02 - DL / 2)                                                         { vnew[I] = -fabs(vnew[I]); y[I] = +1.02 - DL / 2.0; }            //back wall

	if (z[I]<0.0 + DL / 2)														      { wnew[I] = fabs(wnew[I]); z[I] = 0.0 + DL / 2; }      // bootom of the tank
	if (z[I]>Zmax)                                                                   { wnew[I] = -fabs(wnew[I]); z[I] = Zmax; }            //Top of domain





	if (x[I] < 0.5151*z[I] - 4.111)										            { unew[I] = fabs(unew[I]); x[I] = 0.5151*z[I] - 4.111 + DL / 2.0; }                     //box left wall
	if (y[I]<-0.6 && x[I] <= 0.5151*z[I] - 1.453)                                    { vnew[I] = +fabs(vnew[I]); y[I] = -0.6 + DL / 2.0; }            // box front wall
	if (y[I]> 0.6  && x[I] <= 0.5151*z[I] - 1.453)                                      { vnew[I] = -fabs(vnew[I]); y[I] = 0.6 - DL / 2.0; }            // box back wall
	if (x[I] > 0.5151*z[I] - 1.453 && t<0.25 && PTYPE[I]>1)										    { unew[I] = -fabs(unew[I]); x[I] = 0.5151*z[I] - 1.453 ; }                     //box rigth wall
	




	//double h;

	//if (t>.1)  h=0.9*(t-0.1)+0.08;
	//else h=0;

	//if (PTYPE[I]>1  && y[I]>h && x[I]>0.0605 && x[I]<0.0605+5.0*DL)      {unew[I]=-fabs(unew[I])*CC;x[I]=0.0605;}            // gate



	//if(y[I]>=Ymax-DL/2)                                   {vnew[I]=0.0;y[I]=Ymax-DL/2;}                                 //top

	double  COEFF, nx, ny, nz;

	if (x[I] >= -3.5 && x[I] <= 0.0 && z[I] <= -0.515*x[I] + DL / 2)
	{

		COEFF = 1 / sqrt(1 + 0.515*0.515);
		nx = COEFF*0.515;
		ny = 0;
		nz = COEFF;

		REFLECTION(nx, ny, nz);
		z[I] = -0.515*x[I] + DL / 2;
	}


	return;
}
//______________________________  END OF Subroutine (BOUNDARIE) ____________________________

//===========================================================================================
//============ Calculation of dynamic viscosity =================================
//===========================================================================================
void VISCOSITY(double *x_vel, double *y_vel, double *z_vel)
{
	double  *S12, *S13, *S23, *S11, *S22, *S33, d, phi = 0, phi2 = 0, grain_VF, meu0, normal_stress, *p_smooth;
	double **BL, **WL, **PS;
	int aa = 2, kx, ky;
	int kx_max = int((Xmax - Xmin) / aa / DL) + 1;
	int ky_max = int((Ymax - Ymin) / aa / DL) + 1;



	double Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz;
	double w;

	int i, j;


	S11 = new double[TP + 1];
	S22 = new double[TP + 1];
	S33 = new double[TP + 1];
	S12 = new double[TP + 1];
	S13 = new double[TP + 1];
	S23 = new double[TP + 1];
	p_smooth = new double[TP + 1];
	BL = new double*[kx_max + 1];
	WL = new double*[kx_max + 1];
	PS = new double*[kx_max + 1];

	for (int m = 1; m <= kx_max; m++)
	{
		BL[m] = new double[ky_max + 1];
		WL[m] = new double[ky_max + 1];
		PS[m] = new double[ky_max + 1];
	}


	//-------------------------------determining the bed level
	for (kx = 1; kx <= kx_max; kx++)
	{
		for (ky = 1; ky <= ky_max; ky++)
		{
			BL[kx][ky] = Zmin;
			WL[kx][ky] = Zmin;
		}

	}

	for (i = 1; i <= NUM; i++)
	{
		kx = int((x[i] - Xmin) / aa / DL) + 1;
		ky = int((y[i] - Ymin) / aa / DL) + 1;

		if (z[i]>BL[kx][ky] && C[i]>0.5) { BL[kx][ky] = z[i]; PS[kx][ky] = pnew[i]; }
		if (z[i]>WL[kx][ky] && PTYPE[i] == 1) { WL[kx][ky] = z[i]; }
	}


	//---------------------------------------- Strain rate calculation --------------------------------------------

	for (i = 1; i <= NUM; i++)
	{

		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;
		for (l = 2; l <= neighb[i][1]; l++)
		{
			j = neighb[i][l];
			d = DIST(i, j);
			if (i != j && d <= re)
			{
				w = W(d, KTYPE, 2);
				sum1 = sum1 + (x_vel[j] - x_vel[i])*DX(i, j)*w / d / d;
				sum2 = sum2 + (x_vel[j] - x_vel[i])*DY(i, j)*w / d / d;
				sum3 = sum3 + (x_vel[j] - x_vel[i])*DZ(i, j)*w / d / d;

				sum4 = sum4 + (y_vel[j] - y_vel[i])*DX(i, j)*w / d / d;
				sum5 = sum5 + (y_vel[j] - y_vel[i])*DY(i, j)*w / d / d;
				sum6 = sum6 + (y_vel[j] - y_vel[i])*DZ(i, j)*w / d / d;

				sum7 = sum7 + (z_vel[j] - z_vel[i])*DX(i, j)*w / d / d;
				sum8 = sum8 + (z_vel[j] - z_vel[i])*DY(i, j)*w / d / d;
				sum9 = sum9 + (z_vel[j] - z_vel[i])*DZ(i, j)*w / d / d;

				sum10 = sum10 + pnew[j] * w;



			}
		}
		Uxx = (DIM / n0)*sum1;
		Uxy = (DIM / n0)*sum2;
		Uxz = (DIM / n0)*sum3;

		Uyx = (DIM / n0)*sum4;
		Uyy = (DIM / n0)*sum5;
		Uyz = (DIM / n0)*sum6;

		Uzx = (DIM / n0)*sum7;
		Uzy = (DIM / n0)*sum8;
		Uzz = (DIM / n0)*sum9;






		p_smooth[i] = sum10 / n0;
		if (p_smooth[i]<0) p_smooth[i] = 0;

		S11[i] = 0.5*(Uxx + Uxx);
		S12[i] = 0.5*(Uxy + Uyx);
		S13[i] = 0.5*(Uxz + Uzx);
		S22[i] = 0.5*(Uyy + Uyy);
		S23[i] = 0.5*(Uyz + Uzy);
		S33[i] = 0.5*(Uzz + Uzz);


		//II[i] = 0.5*Uxx*Uxx + 0.5*Uyy*Uyy + 0.25*(Uxy + Uyx)*(Uxy + Uyx);
		II[i] = 0.5*(S11[i] * S11[i] + S12[i] * S12[i] + S13[i] * S13[i] + S12[i] * S12[i] + S22[i] * S22[i] + S23[i] * S23[i] + S13[i] * S13[i] + S23[i] * S23[i] + S33[i] * S33[i]);
		//	II[i]= S11[i]*S22[i] +S22[i]*S33[i]+ S11[i]*S33[i] - S12[i]*S12[i] -S13[i]*S13[i]- S23[i]*S23[i] ;
		if (II[i]<0 || II[i] * 0 != 0) II[i] = 0;

		//II=fabs(S11[i]*S22[i]-S12[i]*S12[i]);


	}




	//--------------  Newtonian visc ---------------------
	if (Fluid2_type == 0)
	{
		for (i = 1; i <= NUM; i++)
		{
			if (PTYPE[i] <= 1)MEU[i] = NEU1*Rho1;
			if (PTYPE[i] != 1)MEU[i] = NEU2*Rho2;
		}


		if (TURB>0)
		{
			NEUt[i] = Cs*DL*Cs*DL * 2 * sqrt(II[i]);

			if (NEUt[i] * 0 != 0)  NEUt[i] = 0;
			if (NEUt[i]>1)     NEUt[i] = 1;
		}
	}
	//--------------------- Granular Fluid  -------------------------
	if (Fluid2_type == 1)
	{
		for (i = 1; i <= NUM; i++)
		{


			if (TURB>0)
			{
				NEUt[i] = Cs*DL*Cs*DL * 2 * sqrt(II[i]);

				if (NEUt[i] * 0 != 0)  NEUt[i] = 0;
				if (NEUt[i]>1)     NEUt[i] = 1;
			}

			if (PTYPE[i] == 1)MEU[i] = NEU1*Rho1;
			else
			{

				phi = (C[i] - 0.25)*PHI / (1 - 0.25);
				phi2 = (C[i] - 0.25)*PHI_2 / (1 - 0.25);
				if (C[i] <= 0.25) { phi = 0.00001; phi2 = 0.00001; }
				if (PTYPE[i] <= 0) phi = PHI_bed;


				// --------------------- normal stress calculation ----------------------------
				p_rheo_new[i] = p_smooth[i];

				kx = int((x[i] - Xmin) / aa / DL) + 1;
				ky = int((y[i] - Ymin) / aa / DL) + 1;

				//		normal_stress=(BL[k]-y[i]+DL/2)*(Rho2)*9.81;                                       // normal_stress= Gama.H
				normal_stress = (BL[kx][ky] - z[i] + DL / 2)*(Rho2 - Rho1)*9.81 - (x_vel[i] * x_vel[i] + y_vel[i] * y_vel[i] + z_vel[i] * z_vel[i])*(Rho2 - Rho1) / 2.0;                                       // normal_stress= Gama.H




				if (p_smooth[i] - (WL[kx][ky] - z[i])*Rho1*9.8<0) p_smooth[i] = (WL[kx][ky] - z[i])*Rho1*9.8;
				if (t <= 1) normal_stress = 1.0*(1 - t)*(p_smooth[i] - (WL[kx][ky] - z[i])*Rho1*9.8) + 1.0*(t)*normal_stress;



				//		normal_stress=normal_stress*0.61*1500/Rho2;
				if (normal_stress < 1 || C[i] < 0.5) normal_stress = 1;

				p_rheo_new[i] = normal_stress;


				// --------------------- yeild stress calculation ----------------------------
				//	Inertia[i] = sqrt(II[i])* dg / sqrt(normal_stress / Rho2);        // Free-fall regime
				Inertia[i] = sqrt(II[i])* dg / sqrt(normal_stress / (Rho1*0.47));        // Grain inertia regime
				//	Inertia[i] = sqrt(II[i])* (NEU1*Rho1) / normal_stress ;                 //viscous regime

				grain_VF = 0.65 - (0.65 - 0.25)*Inertia[i];
				phi = phi*grain_VF / 0.65;


				yeild_stress = cohes*cos(phi) + normal_stress*sin(phi);

				if (yeild_stress < 0)yeild_stress = 0;

				visc_max = (yeild_stress*mm + MEU0);

				if (II[i]>0) MEU_Y[i] = yeild_stress*(1 - exp(-mm*sqrt(II[i]))) / 2.0 / sqrt(II[i]);
				else MEU_Y[i] = visc_max;






				// ---------------H-B rhology--------------------------------------

				//meu0 = MEU0;


				// ---------------Non-linear Meu(I) rhology--------------------------------------
				//		meu0=  0.5*    0.36               *normal_stress* dg/ (I0* sqrt(normal_stress/ Rho2)+ sqrt(II[i])*dg);               //free fall
				meu0 = 0.5*(sin(phi2) - sin(phi)) *normal_stress* dg / (I0* sqrt(normal_stress / (Rho1*0.47)) + sqrt(II[i])*dg);               //grain inertia
				//		meu0=0.5*0.36*normal_stress* (NEU1*Rho1)/ (I0* normal_stress+ sqrt(II[i])*(NEU1*Rho1));               //viscous


				// --------------linear Meu(I) rhology ----------------------------------
				//		meu0 = 0.5*(tan(phi2) - tan(phi)) * dg* sqrt(normal_stress * Rho2)      / I0;               //free fall
				//        meu0 = 0.5*(tan(phi2) - tan(phi)) * dg* sqrt(normal_stress * Rho1*0.47) / I0;               //grain inertia
				//		meu0 = 0.5*(tan(phi2) - tan(phi)) * (NEU1*Rho1)                         / I0;               //viscous



				if (II[i] <= 0 || (meu0 * 0) != 0) meu0 = MEU0;

				visc_max = (yeild_stress*mm + meu0);


				MEU[i] = MEU_Y[i] + MEU0*pow(4 * II[i], (N - 1) / 2);

				if (II[i] == 0 || MEU[i]>visc_max) MEU[i] = visc_max;
				if (PTYPE[i] <= 0) MEU[i] = MEU[i] * C[i] + 0.001*(1 - C[i]);

			}
		}

		//---------------------------------- Direct stress calculation method -----------------------------------------

		if (stress_cal_method == 2)
		{
			for (i = 1; i <= NUM; i++)
			{

				double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0;
				for (l = 2; l <= neighb[i][1]; l++)
				{
					j = neighb[i][l];
					d = DIST(i, j);
					if (i != j && d <= re)
					{

						w = W(d, KTYPE, 2);

						double meuij = 2 * MEU[i] * MEU[j] / (MEU[i] + MEU[j]);
						if ((NEUt[i] + NEUt[j])>0) meuij = meuij + 2 * NEUt[i] * RHO[i] * NEUt[j] * RHO[j] / (NEUt[i] * RHO[i] + NEUt[j] * RHO[j]);



						sum1 = sum1 + meuij*(x_vel[j] - x_vel[i])*DX(i, j)*w / d / d;
						sum2 = sum2 + meuij*(x_vel[j] - x_vel[i])*DY(i, j)*w / d / d;
						sum3 = sum3 + meuij*(x_vel[j] - x_vel[i])*DZ(i, j)*w / d / d;

						sum4 = sum4 + meuij*(y_vel[j] - y_vel[i])*DX(i, j)*w / d / d;
						sum5 = sum5 + meuij*(y_vel[j] - y_vel[i])*DY(i, j)*w / d / d;
						sum6 = sum6 + meuij*(y_vel[j] - y_vel[i])*DZ(i, j)*w / d / d;

						sum7 = sum7 + meuij*(z_vel[j] - z_vel[i])*DX(i, j)*w / d / d;
						sum8 = sum8 + meuij*(z_vel[j] - z_vel[i])*DY(i, j)*w / d / d;
						sum9 = sum9 + meuij*(z_vel[j] - z_vel[i])*DZ(i, j)*w / d / d;

					}
				}

				Tau_xx[i] = (DIM / n0) * 2 * sum1;
				Tau_yy[i] = (DIM / n0) * 2 * sum5;
				Tau_zz[i] = (DIM / n0) * 2 * sum9;

				Tau_xy[i] = (DIM / n0)*(sum2 + sum4);
				Tau_xz[i] = (DIM / n0)*(sum3 + sum7);
				Tau_yz[i] = (DIM / n0)*(sum6 + sum8);





			}
		}



	}
	//---------------------------------------------------------------

	
	for (int m = 1; m <= kx_max; m++)
	{
		delete[] BL[m];
		delete[] WL [m];
		delete[] PS [m];
	}
	
	delete[]S11; delete[]S12; delete[]S13; delete[]S22; delete[]S23; delete[]S33; delete[]BL; delete[]WL; delete[]PS; delete[]p_smooth;




	S11 = NULL; S12 = NULL; S13 = NULL; S22 = NULL; S23 = NULL; S33 = NULL; BL = NULL; WL = NULL; PS = NULL; p_smooth = NULL;



}







//===========================================================================================
//========================== Calculation of pressure  ========================================
//===========================================================================================


void PRESSURECALC()
{

	//-------------- Using Equation of State instead of Poisson eq. -------------------

	if (Method == 3)
	{

		for (I = GP + 1; I <= NUM; I++)
		{

			if (PTYPE[I] >= 2)c0 = c02;
			else c0 = c01;




			//	Rho=Rho2*C[I]+Rho1*(1-C[I]);

			if (PTYPE[I] >= 2)Rho = Rho2;
			else Rho = Rho1;


			if (nstar[I]<BETA*n0)

			{
				pnew[I] = 0.0;

			}

			else
			{
				pnew[I] = (c0*c0*Rho / GAMA)*(pow(nstar[I] / n0, GAMA) - 1);    //P=B((rho/rho0)^7-1), B=(10*vmax)^2*Rho0/7
			}
			if (pnew[I]<PMIN)
			{
				pnew[I] = PMIN;
			}
			if (pnew[I]>PMAX)
			{
				pnew[I] = PMAX;
			}

		}

	}
	else
	{
		//---------------------- Using Conjugate gradient method -------------------

		MATRIX();
		CGM(FP + WP);
		//--------replacing result of conjugate gradient method-----------------------
		for (I = 1; I <= NUM; I++)
		{
			if (PTYPE[I]>0)
			{

				pnew[I] = Guess[GP + I];
			}

			if (nstar[I]<BETA*n0)
			{
				pnew[I] = 0.0;

			}

			if (pnew[I]<PMIN)
			{
				pnew[I] = PMIN;
			}
		}
	}


	return;
}
//______________________________  END OF Subroutine (PRESSURECALC) ____________________________

//===========================================================================================
//========================== Surface Tension model ==========================================
//===========================================================================================

//                  THIS SUBROUTINE HAS BEEN DISABLED IN THIS VESION


void Surfacetension()
{

	return;
}
//______________________________  END OF Subroutine (Surfacetension) ____________________________

//===========================================================================================
//=============== Calculation of the volume of fraction if phase II in the mixture =========
//===========================================================================================
void V_FRACTION()
{
	double sum1, sum2, d;

	if (Fraction_method == 1)   //Linear distribution
	{
		for (I = 1; I <= NUM; I++)
		{
			sum1 = 0;
			sum2 = 0;
			for (l = 2; l <= neighb[I][1]; l++)
			{
				J = neighb[I][l];

				if (I != J && PTYPE[J]>0)
				{
					sum1 = sum1 + 1;
					if (PTYPE[J] >= 2)sum2 = sum2 + 1;
				}
			}
			C[I] = sum2 / sum1;
			if (sum1 == 0)C[I] = 0;


		}
	}

	if (Fraction_method == 2)   //Non linear :  Smoothed using the weight funtion
	{
		for (I = 1; I <= NUM; I++)
		{
			sum1 = 0;
			sum2 = 0;
			for (l = 2; l <= neighb[I][1]; l++)
			{
				J = neighb[I][l];
				d = DIST(I, J);
				if (I != J && PTYPE[J]>0)
				{
					sum1 = sum1 + W(d, KTYPE, 2);
					if (PTYPE[J] >= 2)sum2 = sum2 + W(d, KTYPE, 2);
				}
			}
			C[I] = sum2 / sum1;
			if (sum1 == 0)C[I] = 0;


		}
	}

}
//______________________________  END OF FUNCTION FRACTION ____________________________




//********************************************************************************************
//           #       #    ####    ########  ###   ##
//           #       #    ####    ########  ###   ##
//           ##     ##   ##  ##      ##     ###   ##
//           ##     ##   ##  ##      ##     ###   ##
//           ###   ###  ##    ##     ##     ## #  ##
//           ###   ###  ##    ##     ##     ## #  ##
//           ## # # ##  ##    ##     ##     ## ## ##
//           ## # # ##  ##    ##     ##     ## ## ##
//           ##  #  ##  ########     ##     ##  # ##
//           ##  #  ##  ########     ##     ##  # ##
//           ##  #  ##  ##    ##     ##     ##   ###
//           ##  #  ##  ##    ##     ##     ##   ###
//           ##     ##  ##    ##  ########  ##    ##
//           ##     ##  ##    ##  ########  ##    ##
//********************************************************************************************


int main()
{
	printf(" M@@@;   s@@@@  @@@@@@@.      @@,     @@@@@@@@.    .@@@@@@   \n");
	printf(" P@@@@   @@@@@  @@@@@@@@@   .@@@@     @@@B##@@@@  #@@2sX@@@  \n");
	printf(" A@@#@, :@@@@@  @@    :@@,  @@@@@@    @@,    @@@  @@@    :@  \n");
	printf(" R@# @@ M@ @@@  @@@@@@@@@  @@@  @@@   @@@@@##@@r   @@@@@@r   \n");
	printf(" S@@ @@,@@ @@@  @@@@@@@.   @@@;;@@@   @@#9@@@5        @@@@@; \n");
	printf(" @@@ :@@@: @@@  @@        @@@@@@@@@@  @@:  .@@M   @:,    @@@ \n");
	printf(" @@@  @@@  @@@  @@       @@@     r@@, @@:    @@@. @@@@MM@@@. \n");
	printf(" @@@  :@:  @@@  @@       @@:     :@@. @@:      @@   @@@@@@   \n");

	printf("   ____________________________________________________\n");
	printf("  |          MPARS (Mesh-free Particle CFD Simulator)  |\n");
	printf("  |                        Ahmad Shakibaeinia          |\n");
	printf("  |____________________________________________________|\n");
	printf("   ____________________________________________________\n");
	printf("  |             TEST CASE: MPARS EXAMPLE 1             |\n");
	printf("  |              By: AHMAD SHAKIBAEINIA, PhD           |\n");
	printf("  |                   Date: Sept., 2015                |\n");
	printf("  |____________________________________________________|\n");
	printf("\n");



	//------------------------- Defining Dynamic Matrices -------------------------
	x = new double[TP + 1];
	y = new double[TP + 1];
	z = new double[TP + 1];
	p = new double[TP + 1];
	u = new double[TP + 1];
	v = new double[TP + 1];
	w = new double[TP + 1];
	PTYPE = new int[TP + 1];
	xstar = new double[TP + 1];
	ystar = new double[TP + 1];
	zstar = new double[TP + 1];
	ustar = new double[TP + 1];
	vstar = new double[TP + 1];
	wstar = new double[TP + 1];
	pnew = new double[TP + 1];
	phat = new double[TP + 1];
	unew = new double[TP + 1];
	vnew = new double[TP + 1];
	wnew = new double[TP + 1];
	NEUt = new double[TP + 1];
	TURB1 = new double[TP + 1];
	TURB2 = new double[TP + 1];
	p_rheo = new double[TP + 1];
	p_rheo_new = new double[TP + 1];
	C = new double[TP + 1];
	MEU = new double[TP + 1];
	RHO = new double[TP + 1];
	RHOSMOOTH = new double[TP + 1];
	SFX = new double[TP + 1];
	SFY = new double[TP + 1];
	SFZ = new double[TP + 1];
	n = new double[TP + 1];
	nstar = new double[TP + 1];
	neighb = new int*[TP + 1];
	MEU_Y = new double[TP + 1];
	II = new double[TP + 1];
	Inertia = new double[TP + 1];

	Tau_xx = new double[TP + 1];
	Tau_yy = new double[TP + 1];
	Tau_zz = new double[TP + 1];
	Tau_xy = new double[TP + 1];
	Tau_xx = new double[TP + 1];
	Tau_yz = new double[TP + 1];


	for (int m = 1; m <= TP; m++)
		neighb[m] = new int[1500];


	//---------------------- Openning INPUT and OUTPUT files ---------------------------------

	using namespace std;
	fstream in, out, history, tecplot_out;	                    // Openning  input and  output files
	fstream  out1, out2, out3, out4, out5, out6, out7, out8, out9, out10;
	fstream  out11, out12, out13, out14, out15, out16, out17, out18, out19, out20;
	fstream  out21, out22, out23, out24, out25, out26, out27, out28, out29, out30;
	fstream  out31, out32, out33, out34, out35, out36, out37, out38, out39, out40;
	in.open("input.txt", ios::in);
	out.open("output.txt", ios::out);
	history.open("history.txt", ios::out);
	tecplot_out.open("TECPLOT.txt", ios::out);


	out1.open("output1.txt", ios::out);
	out2.open("output2.txt", ios::out);
	out3.open("output3.txt", ios::out);
	out4.open("output4.txt", ios::out);
	out5.open("output5.txt", ios::out);
	out6.open("output6.txt", ios::out);
	out7.open("output7.txt", ios::out);
	out8.open("output8.txt", ios::out);
	out9.open("output9.txt", ios::out);
	out10.open("output10.txt", ios::out);
	out11.open("output11.txt", ios::out);
	out12.open("output12.txt", ios::out);
	out13.open("output13.txt", ios::out);
	out14.open("output14.txt", ios::out);
	out15.open("output15.txt", ios::out);
	out16.open("output16.txt", ios::out);
	out17.open("output17.txt", ios::out);
	out18.open("output18.txt", ios::out);
	out19.open("output19.txt", ios::out);
	out20.open("output20.txt", ios::out);
	out21.open("output21.txt", ios::out);
	out22.open("output22.txt", ios::out);
	out23.open("output23.txt", ios::out);
	out24.open("output24.txt", ios::out);
	out25.open("output25.txt", ios::out);
	out26.open("output26.txt", ios::out);
	out27.open("output27.txt", ios::out);
	out28.open("output28.txt", ios::out);
	out29.open("output29.txt", ios::out);
	out30.open("output30.txt", ios::out);
	out31.open("output31.txt", ios::out);
	out32.open("output32.txt", ios::out);
	out33.open("output33.txt", ios::out);
	out34.open("output34.txt", ios::out);
	out35.open("output35.txt", ios::out);
	out36.open("output36.txt", ios::out);
	out37.open("output37.txt", ios::out);
	out38.open("output38.txt", ios::out);
	out39.open("output39.txt", ios::out);
	out40.open("output40.txt", ios::out);



	double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13, D, MAX = 0, weigth;

	//tecplot_out<<"VARIABLES ='x','y','ptype','u','v','p','c'\n";

	//------------------------INITIALIZATION ---------------------------------------------------
	I = 0, J = 0, K = 0, NUM = 0, l = 0;
	p_count = 1;
	int tec_counter = 0;
	for (I = 1; I <= TP + 1; I++)
	{
		x[I] = 0, y[I] = 0, z[I] = 0, p[I] = 0, u[I] = 0, v[I] = 0, w[I] = 0, PTYPE[I] = 0;
		xstar[I] = 0, ystar[I] = 0, zstar[I] = 0, ustar[I] = 0, vstar[I] = 0, wstar[I] = 0;
		pnew[I] = 0, phat[I] = 0, unew[I] = 0, vnew[I] = 0, wnew[I] = 0, NEUt[I] = 0, TURB1[I] = 0, TURB2[I] = 0, n[I] = 0, nstar[I] = 0, SFX[I] = 0, SFY[I] = 0, SFZ[I] = 0, RHO[I] = 0, C[I] = 0, RHOSMOOTH[I] = 0, MEU_Y[I] = 0, p_rheo[I] = 0, p_rheo_new[I] = 0;
	}

	//-------------------------- defining LAMDA -----------------------------------------

	if (KTYPE == 1)   lambda = 0.22143*pow(re, 2);
	if (KTYPE == 2)   lambda = 0.16832*pow(re, 2);
	if (KTYPE == 5)   lambda = 0.250*pow(re, 2);
	if (KTYPE == 6)   lambda = (3.0 / 14.0)*pow(re, 2);
	//----------------------- Reading Input file ---------------------------------------
	printf("     READING INPUT FILE . . .\n");

	Xmin = 99999999;
	Ymin = 99999999;
	Zmin = 99999999;
	Xmax = -99999999;
	Ymax = -99999999;
	Zmax = -99999999;
	FP = 0;
	WP = 0;
	GP = 0;
	TP = 0;

	in >> NUM;
	for (I = 1; I <= NUM; I++)
	{
		in >> x[I];
		in >> y[I];
		in >> z[I];
		in >> PTYPE[I];
		in >> u[I];
		in >> v[I];
		in >> w[I];
		in >> p[I];

		if (PTYPE[I]<0)  GP++;
		if (PTYPE[I] == 0) WP++;
		if (PTYPE[I]>0)  FP++;

		if (x[I]<Xmin) Xmin = x[I];
		if (x[I]>Xmax) Xmax = x[I];
		if (y[I]<Ymin) Ymin = y[I];
		if (y[I]>Ymax) Ymax = y[I];
		if (z[I]<Zmin) Zmin = z[I];
		if (z[I]>Zmax) Zmax = z[I];
	}
	TP = FP + GP + WP + SP;

	Xmin = Xmin - DL;
	Ymin = Ymin - DL;
	Zmin = Zmin - DL;

	Xmax = Xmax + DL;
	Ymax = Ymax + DL;
	Zmax = Zmax + 6 * DL;
	printf("     PRE- ITERATION CALCULATIONS... \n");
	//----------------------------------- Density assigning -----------------------------------------------------

	for (I = 1; I <= NUM; I++)
	{

		if (PTYPE[I] >= 2)RHO[I] = Rho2;
		else            RHO[I] = Rho1;

	}
	//--------------------------- Calculation of initial Particle number---------------------------------------

	neighbour_cuda_1();

	for (I = 1; I <= NUM; I++)
	{

		n[I] = PNUM(I, KTYPE, DIM, NUM);

		pnew[I] = p[I];
		if (n[I]>MAX)   MAX = n[I];
	}
	n0 = MAX;
	//printf("(%f)   ", n[144564]); getchar();

	printf("     The initial particle number density = %f\n", n0);


	//-------------------------------------------------------------------------------------
	//------------------------------- Time Iteration -------------------------------------------
	//-------------------------------------------------------------------------------------
	printf("\n");
	printf("     ITERATION STARTED  (each bar represents 1 time step) \n");
	t = 0;

	// -------------------------Calculation of CPU time	----------------------
	time_t start, end;
	double dif;
	time(&start);
	//------------------------------------------------------------
	for (int Tstep = 1; t <= T; Tstep++)
	{
		DTcalculation();
		if (((Tstep - 1) / 10.00) == int((Tstep - 1) / 10.00))
		{
			printf("     t=%f Sec     ", t);
			printf("\n");
			printf("                  ", t);
		}


		NUM = FP + WP + GP;


		if (Method != 3)
		{
			for (I = 1; I <= NUM; I++)
			{
				n[I] = PNUM(I, KTYPE, DIM, NUM);
			}
		}

		//	if (((Tstep-1)/4.00)==int((Tstep-1)/4.00))   // Setting nighboring list each 4 time step
		//	{
		neighbour_cuda_1();

		//	}
		//	if (TURB==1) SPS();
		V_FRACTION();


		//----------------------------- calculation of surface tention force --------------------------
		if (surface_method>0)Surfacetension();



		//-------------------------------------Prediction--------------------------------------

		//-------------------------------------Prediction--------------------------------------


		//----------------------------- calculation of dynamic visc --------------------------

		visc_error = 0.0, visc_ave = 0.0;

		VISCOSITY(u, v, w);

		//	 for(int visc_itr=1; visc_itr<=visc_itr_num; visc_itr++) // iteration for calculation of stress dependent viscosity
		//	 {

		for (I = 1; I <= NUM; I++)

		{
			if (PTYPE[I] >= 2)
			{
				if (C[I] > 0.5) RHO[I] = Rho2;
				else RHO[I] = C[I] * Rho2 + (1 - C[I]) * 2650;
			}

			if (PTYPE[I] <= 0)

			{

				xstar[I] = x[I];

				ystar[I] = y[I];

				zstar[I] = z[I];

				ustar[I] = u[I];

				vstar[I] = v[I];

				wstar[I] = w[I];
			}

			else

			{

				sum1 = 0.0;
				sum2 = 0.0;
				sum3 = 0.0;
				sum4 = 0.0;
				sum5 = 0.0;
				sum6 = 0.0;
				sum7 = 0.0;
				sum8 = 0.0;
				sum9 = 0.0;
				sum10 = 0.0;
				sum11 = 0.0;
				sum12 = 0.0;
				sum13 = 0.0;


				for (l = 2; l <= neighb[I][1]; l++)

				{

					J = neighb[I][l];

					D = DIST(I, J);



					if (I != J && D>0 && D <= re)

					{


						weigth = W(D, KTYPE, 2);


						if (PTYPE[I] == 1) NEU = 2 * MEU[I] * MEU[J] / (MEU[I] + MEU[J]) / Rho1;

						else NEU = 2 * MEU[I] * MEU[J] / (MEU[I] + MEU[J]) / Rho2;
						if ((NEUt[I] + NEUt[J])>0) NEU = NEU + (2 * NEUt[I] * RHO[I] * NEUt[J] * RHO[J] / (NEUt[I] * RHO[I] + NEUt[J] * RHO[J])) / RHO[I];



						sum1 = sum1 + (pnew[J] - phat[I])*DX(I, J)*weigth / D / D;
						sum2 = sum2 + (pnew[J] - phat[I])*DY(I, J)*weigth / D / D;
						sum11 = sum11 + (pnew[J] - phat[I])*DZ(I, J)*weigth / D / D;

						sum3 = sum3 + weigth*(u[J] - u[I])*NEU;
						sum4 = sum4 + weigth*(v[J] - v[I])*NEU;
						sum12 = sum12 + weigth*(w[J] - w[I])*NEU;
						if (sum3 * 0 != 0) {
							//	printf("prediction,%f,%f,%f, %f,%f,%f,%f, %f,%d,%d\n", weigth, D, u[I], u[J], NEU, MEU[I], MEU[J], NEUt[I], I, J); getchar();
						}


						//		sum5 = sum5 + weigth*(u[J] - u[I]);
						//		sum6 = sum6 + weigth*(v[J] - v[I]);
						//		sum13 = sum13 + weigth*(w[J] - w[I]);



						/*		sum7 = sum7 + (Tau_xx[J] - Tau_xx[I])*DX(I, J)*weigth / D / D;
						sum8 = sum8 + (Tau_yy[J] - Tau_yy[I])*DY(I, J)*weigth / D / D;


						sum9 = sum9 + (Tau_xy[J] - Tau_xy[I])*DX(I, J)*weigth / D / D;
						sum10 = sum10 + (Tau_xy[J] - Tau_xy[I])*DY(I, J)*weigth / D / D;*/



					}

				}

				//	if (stress_cal_method == 1)
				//	{

				ustar[I] = u[I] + gx*DT + 2 * DIM * DT*(sum3) / (lambda*n0) - (1 - relaxp)*(DIM*DT / n0 / RHO[I])*sum1 + DT*SFX[I] / RHO[I];
				vstar[I] = v[I] + gy*DT + 2 * DIM  * DT*(sum4) / (lambda*n0) - (1 - relaxp)*(DIM*DT / n0 / RHO[I])*sum2 + DT*SFY[I] / RHO[I];
				wstar[I] = w[I] + gz*DT + 2 * DIM  * DT*(sum12) / (lambda*n0) - (1 - relaxp)*(DIM*DT / n0 / RHO[I])*sum11 + DT*SFZ[I] / RHO[I];

				//	}
				//	else
				//	{
				// this part is not 3D
				//		ustar[I] = u[I] + gx*DT - (1 - relaxp)*(DIM*DT / n0 / RHO[I])*sum1 + DT*SFX[I] / RHO[I] + (2.0*DT / n0 / RHO[I])* (sum7 + sum10);
				//		vstar[I] = v[I] + gy*DT - (1 - relaxp)*(DIM*DT / n0 / RHO[I])*sum2 + DT*SFY[I] / RHO[I] + (2.0*DT / n0 / RHO[I])* (sum8 + sum9);

				//	}

				xstar[I] = x[I] + DT*ustar[I];
				ystar[I] = y[I] + DT*vstar[I];
				zstar[I] = z[I] + DT*wstar[I];




			}

		}

		//	 VISCOSITY(ustar,vstar);
		//	 if (visc_error<0.01*visc_ave) visc_itr=visc_itr_num;
		//	 }




		//------------ Particle collision and calculation of new particle number density-------

		if (int(Tstep / 2) == 1.0 / 2 * Tstep)
		{
			for (I = 1; I <= NUM; I++)
			{
				COLLISION(I, coll);
			}
		}
		else
		{
			for (I = NUM; I >= 1; I--)
			{
				COLLISION(I, coll);
			}

		}


		for (I = 1; I <= NUM; I++)
		{

			p_rheo[I] = p_rheo_new[I];
			nstar[I] = PNUMSTAR(I, KTYPE, 2, NUM);


		}
		//-------------------------------------Pressure calculation--------------------------------------

		PRESSURECALC();
		//------------------------------------Correction----------------------------------------

		for (I = 1; I <= GP + WP; I++)
		{
			BC(WBC);
		}

		//******************** Calculation of Phat ****************
		for (I = 1; I <= NUM; I++)
		{
			double min = 999999999999999;
			for (l = 2; l <= neighb[I][1]; l++)
			{
				J = neighb[I][l];

				if (pnew[J]<min)min = pnew[J];
			}
			phat[I] = min;
		}
		//**************** Calculation of pressure gradient *****************

		for (I = GP + WP + 1; I <= NUM; I++)
		{

			sum1 = 0;
			sum2 = 0;
			sum3 = 0;
			sum4 = 0;
			sum5 = 0;
			sum6 = 0;


			if (KHcorrection == 1)
			{
				for (l = 2; l <= neighb[I][1]; l++)
				{
					J = neighb[I][l];
					D = DISTSTAR(I, J);
					if (I != J && D <= re)
					{
						sum1 = sum1 + (pnew[J] + pnew[I] - 2.0*phat[I])*DXSTAR(I, J)*W(D, KTYPE, 2) / D / D;
						sum2 = sum2 + (pnew[J] + pnew[I] - 2.0*phat[I])*DYSTAR(I, J)*W(D, KTYPE, 2) / D / D;
						sum3 = sum3 + (pnew[J] + pnew[I] - 2.0*phat[I])*DZSTAR(I, J)*W(D, KTYPE, 2) / D / D;
					}
				}
			}
			else
			{
				for (l = 2; l <= neighb[I][1]; l++)
				{
					J = neighb[I][l];
					D = DISTSTAR(I, J);
					if (I != J && D <= re)
					{

						sum1 = sum1 + (pnew[J] - phat[I])*DXSTAR(I, J)*W(D, KTYPE, 2) / D / D;
						sum2 = sum2 + (pnew[J] - phat[I])*DYSTAR(I, J)*W(D, KTYPE, 2) / D / D;
						sum3 = sum3 + (pnew[J] - phat[I])*DZSTAR(I, J)*W(D, KTYPE, 2) / D / D;
					}
				}
			}
			Rho = RHO[I];
			unew[I] = ustar[I] - relaxp*(DIM*DT / n0 / Rho)*sum1;
			vnew[I] = vstar[I] - relaxp*(DIM*DT / n0 / Rho)*sum2;
			wnew[I] = wstar[I] - relaxp*(DIM*DT / n0 / Rho)*sum3;


			//----------- Damper ----------------------------

			if (fabs(unew[I])>2.0*VMAX) unew[I] = VMAX;
			if (vnew[I]>2.0*VMAX) vnew[I] = 2.0*VMAX;
			if (vnew[I]<-2.0*VMAX) vnew[I] = -2.0*VMAX;

			//------------------------------------------------

		}
		//---------------------------------- Moving particles -----------------------------------
		for (I = GP + WP + 1; I <= NUM; I++)
		{

			x[I] = x[I] + unew[I] * DT;
			y[I] = y[I] + vnew[I] * DT;
			z[I] = z[I] + wnew[I] * DT;
			BOUNDARIE();                       // Avoid penetration of particle to boundaries

		}

		//---------------------------------- Inflow ---------------------------------------------------


		if (Inflow_type>0)
		{
			INFLOW();
			SORT();
		}
		//-----------------------------------Aplying the pair-wise Collision ----------------------------------------------

		if (int(Tstep / 2) == 1.0 / 2 * Tstep)
		{
			for (I = 1; I <= NUM; I++)
			{
				COLLISION2(I, coll);
			}
		}
		else
		{
			for (I = NUM; I >= 1; I--)
			{
				COLLISION2(I, coll);
			}

		}


		//--------------------------------- print results ---------------------------------------------
		/*		if (((Tstep-1)/100.000)==int((Tstep-1)/100.000))
		{


		//	out<<" t= "<<t<<"\n";
		//	out<<"------------------------------------------------------------------\n";
		//	out<<"  I     x      y       u       v       p\n";
		//	out<<"-------------------------------------------------------------------\n";
		out<<t<<"\n";
		out<<NUM<<"\n";

		for (I=1;I<=NUM;I++)
		{

		out<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"\n";
		}
		}


		*/



		if (int(t * 100) == int(tec_out_intrval*tec_counter * 100))
		{

			tecplot_out << "VARIABLES  =  x,y,z,ptype,u,v,w,p,prheo,c, log_Meu, Meu_y\n";
			tecplot_out << " ZONE T=\"" << tec_out_intrval*tec_counter << "\"\n";

			for (I = 1; I <= NUM; I++)
			{


				if (PTYPE[I] >= 0)  tecplot_out << x[I] << "  " << y[I] << "	" << z[I] << "	" << PTYPE[I] << "  " << unew[I] << "  " << vnew[I] << " " << wnew[I] << " " << pnew[I] << "	" << p_rheo[I] << "  " << C[I] << "	" << log10(MEU[I]) << "	" << MEU_Y[I] << "\n";
			}
			tecplot_out << "\n";
			tec_counter++;
		}


		/*

		if (int(t*100)==int(output_time1*100) && p_count==1 )
		{
		out1<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out1<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out1.close();
		}
		if (int(t*100)==int(output_time2*100) && p_count==2 )
		{
		out2<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out2<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out2.close();
		}
		if (int(t*100)==int(output_time3*100) && p_count==3 )
		{
		out3<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out3<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out3.close();
		}
		if (int(t*100)==int(output_time4*100) && p_count==4 )
		{
		out4<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out4<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out4.close();
		}
		if (int(t*100)==int(output_time5*100) && p_count==5 )
		{
		out5<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out5<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out5.close();
		}
		if (int(t*100)==int(output_time6*100) && p_count==6 )
		{
		out6<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out6<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out6.close();
		}
		if (int(t*100)==int(output_time7*100) && p_count==7 )
		{
		out7<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out7<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out7.close();
		}
		if (int(t*100)==int(output_time8*100) && p_count==8 )
		{
		out8<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out8<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out8.close();
		}
		if (int(t*100)==int(output_time9*100) && p_count==9 )
		{
		out9<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out9<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out9.close();
		}

		if (int(t*100)==int(output_time10*100) && p_count==10 )
		{
		out10<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out10<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out10.close();
		}
		if (int(t*100)==int(output_time11*100) && p_count==11 )
		{
		out11<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out11<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out11.close();
		}
		if (int(t*100)==int(output_time12*100) && p_count==12 )
		{
		out12<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out12<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out12.close();
		}
		if (int(t*100)==int(output_time13*100) && p_count==13 )
		{
		out13<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out13<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out13.close();
		}
		if (int(t*100)==int(output_time14*100) && p_count==14 )
		{
		out14<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out14<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out14.close();
		}
		if (int(t*100)==int(output_time15*100) && p_count==15 )
		{
		out15<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out15<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out15.close();
		}
		if (int(t*100)==int(output_time16*100) && p_count==16 )
		{
		out16<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out16<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out16.close();
		}
		if (int(t*100)==int(output_time17*100) && p_count==17 )
		{
		out17<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out17<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out17.close();
		}
		if (int(t*100)==int(output_time18*100) && p_count==18 )
		{
		out18<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out18<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out18.close();
		}
		if (int(t*100)==int(output_time19*100) && p_count==19 )
		{
		out19<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out19<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out19.close();
		}

		if (int(t*100)==int(output_time20*100) && p_count==20 )
		{
		out20<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out20<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out20.close();
		}
		if (int(t*100)==int(output_time21*100) && p_count==21 )
		{
		out21<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out21<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out21.close();
		}
		if (int(t*100)==int(output_time22*100) && p_count==22 )
		{
		out22<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out22<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out22.close();
		}
		if (int(t*100)==int(output_time23*100) && p_count==23 )
		{
		out23<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out23<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out23.close();
		}
		if (int(t*100)==int(output_time24*100) && p_count==24 )
		{
		out24<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out24<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out24.close();
		}
		if (int(t*100)==int(output_time25*100) && p_count==25 )
		{
		out25<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out25<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out25.close();
		}
		if (int(t*100)==int(output_time26*100) && p_count==26 )
		{
		out26<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out26<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out26.close();
		}
		if (int(t*100)==int(output_time27*100) && p_count==27 )
		{
		out27<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out27<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out27.close();
		}
		if (int(t*100)==int(output_time28*100) && p_count==28 )
		{
		out28<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out28<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out28.close();
		}
		if (int(t*100)==int(output_time29*100) && p_count==29 )
		{
		out29<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out29<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out29.close();
		}

		if (int(t*100)==int(output_time30*100) && p_count==30 )
		{
		out30<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out30<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out30.close();
		}
		if (int(t*100)==int(output_time31*100) && p_count==31 )
		{
		out31<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out31<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out31.close();
		}
		if (int(t*100)==int(output_time32*100) && p_count==32 )
		{
		out32<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out32<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out32.close();
		}
		if (int(t*100)==int(output_time33*100) && p_count==33 )
		{
		out33<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out33<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out33.close();
		}

		if (int(t*100)==int(output_time34*100) && p_count==34 )
		{
		out34<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out34<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out34.close();
		}
		if (int(t*100)==int(output_time35*100) && p_count==35 )
		{
		out35<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out35<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out35.close();
		}
		if (int(t*100)==int(output_time36*100) && p_count==36 )
		{
		out36<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out36<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out36.close();
		}
		if (int(t*100)==int(output_time37*100) && p_count==37 )
		{
		out37<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out37<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out37.close();
		}
		if (int(t*100)==int(output_time38*100) && p_count==38 )
		{
		out38<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out38<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out38.close();
		}
		if (int(t*100)==int(output_time39*100) && p_count==39 )
		{
		out39<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out39<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out39.close();
		}

		if (int(t*100)==int(output_time40*100) && p_count==40 )
		{
		out40<<t<<"\n";p_count++;
		for (I=1;I<=NUM;I++)
		out40<<x[I]<<"  "<<y[I]<<"  "<<PTYPE[I]<<"  "<<unew[I]<<"  "<<vnew[I]<<"  "<<pnew[I]<<"  "<<C[I]<<"\n";
		out40.close();
		}*/


		//--------------------- Prepare data for new time step----------------------------------------------------------

		for (I = 1; I <= NUM; I++)
		{
			p[I] = pnew[I];
			u[I] = unew[I];
			v[I] = vnew[I];
			w[I] = wnew[I];


			if (p[I] * 0 != 0 || u[I] * 0 != 0 || v[I] * 0 != 0)
			{
				//cout<<"Error Check"<<p[I]<<" "<<u[I]<< " "<<v[I]<<endl;
				//int a;
				//cin>>a;
				printf("ERROR#1: ERROR in particle %d , x=%f, y=%f, p=%f\n", I, x[I], y[I], p[I]);
				getchar();
			}

			if (x[I] >= Xmax || y[I] >= Ymax || z[I] >= Zmax || x[I] <= Xmin || y[I] <= Ymin || z[I] <= Zmin)
			{
				printf("ERROR#2: ERROR in particle %d , x=%f, y=%f\n", I, x[I], y[I], z[I]);
			}
		}


		t = t + DT;

		if (Tstep == 1)
		{
			time(&end);
			dif = difftime(end, start);

			printf("Estimated running time:%.2fsec (per time step)\n", dif);
			printf("                                      OR %.2f hr (per 1sec simulation)\n", dif / DT / 3600);
			printf("                  ", t);


		}
		printf("|");


	}
	// --------------------------------End of Time Loop --------------------------------------------------


	return;
}

