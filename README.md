# Overview 
- **Contains: *A highly optimised Parallelised patterns-matching algorithm implementation in C++ using the MPI distributed computing standard*** 
- **Problem:** Given a text string and a set of patterns to be matched, find all occurrences of every pattern in the text
- **Strategy:** Executing a linear-time pattern matching algorithm on multiple patterns concurrently
- MPI C++ file: *lab4_mpi.cpp*
- **Algorithm:** Described in included problem statement doc
- I/O: I/O format can be understood from the included header files and sample testcase files.
- Literature reference used included in repository. 
- Assignment attempted as a part of coursework requirements in *COL380: Introduction to Parallel Programming and Distributed Computing (Sem-II, 2018-19)* (Instructor: Prof. Subodh V. Sharma) at Indian Institute of Technology (IIT), Delhi. 
- The problem statement is included. The following sections describe the assignment submission requirements and how to use the starter codes.

## Problem Statement & Starter Codes: col380_lab4_suite
- Problem Statement: Implement Parallel Periodic Pattern Matching using MPI
- Cloned from: https://github.com/dvynjli/col380_lab4_suite/

## Directories and files
`testcase/`: contains python script `gen_testcase.py` for sample testcase generation  
`lab4_io.h` and `lab4_io.c`: functions to read text & pattern from file and check the correctness of the result  
`main_mpi.c`: function `main()`  
`lab4_mpi.h`: header file for the functions to be implemented  
`lab4_mpi.c`: implement the function in this file  
Refer to respective files for further details.  
**Do not change the directory structure and prototype of functions.**

## Building and Executing
```
mpicc -lm main_mpi.c lab4_mpi.c lab4_io.c -o ppm
```
#### Command Line Arguments
The program takes one command line arguments:
- arg: input filename (consist text and patterns)  

To run the program:
```
mpirun -np 4 ppm <input filename>
```
Example:
```
mpirun -np 4 ppm testcase/testcase_10000_10
```

## Generating testcases
Script `gen_testcase.py` generates testcases as per the parameters and output the generated testcase in file `testcase_<n>_<num_patterns>` in the desired format. You might need to change the values of variables `n` and `num_patterns` in the script. Read the comments in the script for more information.
```
python3 gen_testcase.py
```

## Input-Output Specifications
#### Input dataset specifications
- n : length of text
- text : text in which pattern is to be matched
- num_patterns : #patterns to be matched in the text
- m_set : lengths of patterns in pattern_set
- p_set : periods of patterns in pattern_set
- pattern_set : set of patterns to be matched

The first line of the input file contains `n` followed by `num_patterns`. The second line contains the `text`. Third and fourth lines contain length of patterns `m_set`, and period of patterns `p_set` respectively. Next `num_patterns` lines contain the patterns to be matched. All the values in one line are space separated.  

#### Output Specification
Your program should find all the matches of all the patterns in the given text and store the results in the variables given in the program. We will check the correctness by calling the functions from the program. You should compute following:  
- match_counts : #match of `pattern[i]` in text
- matches : set of all matches of each `pattern[i]` in text

Refer to `lab4_mpi.h` for more details. **Your program should not output anything on `stdout`.**  

## Submission Instructions
- You are supposed to submit only one file named `lab4_mpi.c`. Please make sure all the functions you have used are in this file.
- Do not submit other files
- Your code should build and execute as per the instructions given above. Please make sure that your code doesn't need any Makefile.
- Your program should not output anything in `stdout`.

We will not consider the submissions that don't comply with these guidelines.

