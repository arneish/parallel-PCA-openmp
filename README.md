# col380_lab2_suite

Problem Statement: Implement Principal Component Analysis with Singular Vector Decomposition

## Directories and files
`testcase/`: contains python script `gen_testcase.py` for sample testcase generation  
`lab2_io.h` and `lab2_io.c`: functions to read matrix from file and check the correctness of the result  
`main_omp.c`: function `main()`  
`lab2_omp.h`: header file for the functions to be implemented  
`lab2_omp.c`: implement the functions in this file  
Refer to respective files for furthur details.  
**Note: Do not change the directory structure and prototype of functions.**

## Building
```
g++ -fopenmp lab2_io.c lab2_omp.c main_omp.c -o pca
```

## Running
The program takes two command line arguments:
- arg1: input filename (consist M, N and D)  
- arg2: retention (percentage of information to be retained by PCA) 

Please refer to `main_omp.c` for more details.

## Generating testcases
Script `gen_testcase.py` generates testcases as per the parameters and output the generated testcase in file `testcase_<M>_<N>` in the desired format. You might need to change the values of variables `M` and `N` in the script. Read the comments in the script for more information.
```
python3 gen_testcase.py
```
