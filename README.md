# Parallel MPS Implementation

This repository is a part of a project done under the guidance of Prof. Ahmad Shaekeibinia

## Neighbour Algorithm

The neighbour module can be imported via headers like so - `#include<neighb.cu>`

-------------
The first function available is `neighbour_cuda_1(args).` 
The required arguments are :
- int* - `x`, `y`, `z`
- int - `Xmax`, `Xmin`, `Ymax`, `Ymin`, `Zmax`, `Zmin`, `re`, `DELTA`, `NUM`, `MAX_NEIGHB`
- int** - `neighb`

The function modifies the `neighb` array which is passed by reference in the arguments and populates it with the required values of particle IDs and number of neighbours. For any particle ID `i` the number of neighbours can be found by querying `neighb[i][1]`. And all the neighbour particle IDs can be iterated over by doing the following. 

```C++
for(int j=0; j<neighb[i][1]; j++){
  neighb[i][j+2];
}
```
----
There is another function available which is called `neighbour_cuda_2(args)`.
The required arguments are :
- int* - `x`, `y`, `z`, `particleHash`, `particleID`, `cellStart`, `cellEnd`
- int - `Xmax`, `Xmin`, `Ymax`, `Ymin`, `Zmax`, `Zmin`, `re`, `DELTA`, `NUM`, `MAX_NEIGHB`

This function modifies `particleHash`, `particleID`, `cellStart`, `cellEnd` just like the previous function modifies `neighb` but these 4 arrays combined take less space on the global memory than the `neighb` array and hence constitute a preferable mode of implementation. Together, they make up O(2N + 2M) space where N is the total number of particles and M is the total number of cells. Whereas the `neighb` takes up O(N * MAX_NEIGHB) space.

The neighbours in this can can be looped over for all particles like so:
- For an index `i` find the Particle Id and particle cell number from `particleiD` and `particleHash`. 
- Find the coordinates of the cell in terms of `i`, `j`, and `k`. Here we use `Cnum = (i-1) + (j-1)*ncx + (k-1)*ncx*ncy`. 
- Find the neighbouring cell numbers and iterate over the particles in those cells using `cellStart`, `cellEnd`, and `particleId`. `cellstart` and `cellEnd` are already populated according to the key-sorted `particleHash` with `particleId` as the key-array.



