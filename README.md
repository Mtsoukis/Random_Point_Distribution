# Random Point Distribution in Julia

## Overview \& Introduction

This repository outlines some thoughts on trying to find the closest and furthest points in a plain. The basic question is the following: given n points within a unit square, what is the fastest and most efficient way to find the two closest and furthest points? Two examples, using the number of points, n=20,1000, and the number of repetitions, r=50, are shown below. 

![alt text](https://github.com/Mtsoukis/Random_Point_Distribution/blob/main/Figures/rpdj20.png)
![alt text](https://github.com/Mtsoukis/Random_Point_Distribution/blob/main/Figures/rpdj1000.png)

The main file is Random_Point_Distribution.jl 

## Results

The brute-force algorithm runs in O(n²) time, while the divide-and-conquer convex-hull algorithm runs in O(n log n) time but with higher overhead. They run in equivalent time at around 500 points; the brute force method is faster for lower number of points. 

![alt text](https://github.com/Mtsoukis/Random_Point_Distribution/blob/main/Figures/compare.png)

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