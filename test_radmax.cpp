#include <iostream> 
#include <ctime> 
#include <cstdlib>
using namespace std;

double *x, *y, *z;

void create_particles(int NUM){
    x = (double *)malloc(sizeof(int)*NUM);
    y = (double *)malloc(sizeof(int)*NUM);
    z = (double *)malloc(sizeof(int)*NUM);


    srand((unsigned)time(0)); 
    float lowest=1, highest=10; 
    float range=(highest-lowest)+1; 
    for(int index=0; index<NUM; index++){ 
        x[index] = lowest+float(range*rand()/(RAND_MAX + 1.0)); 
        y[index] = lowest+float(range*rand()/(RAND_MAX + 1.0)); 
        z[index] = lowest+float(range*rand()/(RAND_MAX + 1.0)); 
    } 


}

int main() 
{ 
    create_particles(200);
    for(int i=0; i<=20; i++){
        cout<<x[i]<<y[i]<<z[i]<<endl;
    }
}
