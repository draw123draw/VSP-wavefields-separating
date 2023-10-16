The catalog consts of:
1.The VSP data in the paper.
2.The first break time picked for the VSP data. The unit is the number of the sampling point.
3.Source code of the algorithm.

The first break time data is only for estimating the delta_t in formula (15) in the paper. The delta_t can also be editted in the 66th row in find_window_size.cpp 
find_window_size.cpp is for calculating the window size to input to slope_scanning_median.cpp
slope_scanning_median.cpp is for calculating the slopes and median filtering the VSP data.

All parameter variables are global variables defined after the declaration of header file.

To run the code, a MinGW compiler is required. Generate the executable file with the Makefile and then run the executable builded by Makefile. 

How to use the Makefile: input "make" or "mingw32-make" in the terminal corresponding to the directory of this project. Remember to have a MinGW compiler in the system path.
