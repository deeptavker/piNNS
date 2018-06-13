# Parallel MPS Implementation

This repository is a part of a project done under the guidance of Prof. Ahmad Shaekeibinia

## Neighbour Algorithm

The neighbour module can be imported via headers. 

-------------
The first function available is `neighbour_cuda_1(args).` 
The required arguments for now are :
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
