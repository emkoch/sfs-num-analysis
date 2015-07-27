# sfs-num-analysis
This repository contains R scripts for numerical analysis of the Wright-Fisher diffusion 
under pretty general conditions of demography and selection in a single population. 
This is based off of the theory and scheme described by 
[Evans et al. (2007)](http://www.sciencedirect.com/science/article/pii/S0040580906000785)


Individual runs can be performed using 

```
Rscript run_solution_single.r <>
```

The necessary arguments are described in the comments to run_solution_single.r. All 
time parameters are given in coalescent units according to the initial population 
size. 

## Setting the allele frequency grid
The numerical algorithm at present requires that when the step size in the frequency 
grid changes that it does so by a factor of two. A nonuniform grid is greatly desired 
because of the very steep gradient in the spectrum at low frequencies. Another limitation 
of the grid is that it should start at 0 and end at 1. I have not generally solved the 
problem of finding a grid that satisfies these two requirements. Instead, I've found 
two grids that seem to work well enough. The files `x_setup.r` and `x_setup_test.r` 
provide a large and small grid respectively. The large grid works fine for relatively 
we selection, but when S becomes large it's necessary to use the small one. Using the
small grid causes the program to run slower. 