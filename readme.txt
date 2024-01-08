The catalog consts of:
1.The VSP data in the paper.
2.The first break time picked for the VSP data. The unit is the number of the sampling point.
3.Source code of the algorithm.

The first break time data is only for estimating the delta_t in formula (15) in the paper. The delta_t can also be editted in the 86th row in 'find_window_size.cpp' or 66th row in 'find_window_size_fftw.cpp'
'find_window_size.cpp' is for calculating the window size to input to slope_scanning_median.cpp. 'find_window_size_fftw.cpp' is the accelerated version using the 'fftw3' library, however, this can only be used in a system of Windows.
'slope_scanning_median.cpp' is for calculating the slopes and median filtering the VSP data.

All parameter variables are global variables defined after the declaration of header file.

To run the code, a MinGW compiler package is required. 
Note that only systems of Linux (Ubuntu 23.04) and Windows (10,11) are tested to be available for running the code. Ensure there is a 'make' or 'mingw32-make.exe' in the system route. 
Generate the executable file with the Makefile. Then run 'find_window_size(.exe)' to find the suitable window size and 'slope_scanning_median(.exe)' to separate wavefields. 
