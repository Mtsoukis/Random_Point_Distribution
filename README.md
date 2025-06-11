# Random Point Distribution in Julia

## Overview \& Introduction

This repository outlines some thoughts on trying to find the closest and furthest points in a plain. The basic question is the following: given n points within a unit square, what is the fastest and most efficient way to find the two closest and furthest points?

## Repository Structure

- `IO_Rust_Julia.ipynb`  
  A Jupyter notebook containing all data preprocessing, estimation routines, and result visualizations.
- `IO_Rust_Julia.jl`  
  A Julia Script. It is not very optimized (hopefully soon), as it just takes and simplifies the jupyter notebook. 
- `rust_data_2020.csv`  
  Dataset with fields:
  - `bus_id`: Bus identifier  
  - `period_id`: Time period  
  - `y_it`: Replacement decision (1 = replace, 0 = keep)  
  - `x_it`: Mileage state  
- `Rust_Problems.pdf`  
  Problem set.
- `generate_data.jl`  
  Generates bus data
- `README.md`  
  This file.

## Results

The notebook outputs parameter estimates, computation times, and EV function visualizations for both NFXP and MPEC. I show that using JuMP indeed does solve this problem; Chris used AMPL. 

## Extensions/ To Do

I plan to add the following:

- What if instead of a unit square- I did this within a unit circle or triangle? Would threshold for changing algorithms be higher or lower?
- Would the strategy outlined in main.ipynb still be the most efficient?
- What happens if objective is to find m closest/farthest points? Using Euclidean Distance- perimeter.
- Suppose I want to change my definition of farthest points; I now want to select the m points forming the largest polygon possible. How does the algorithm change; How does the complexity change; how do the results change; how does this change as we increase the number of total points n and the number of points we select m? In what fraction of instances are the same points selected as when using Euclidean Distance?


## License

This project is released under the MIT License.

## Author

**Marios Tsoukis**