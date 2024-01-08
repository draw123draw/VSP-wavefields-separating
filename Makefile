all:find_window_size separation
CC:=g++
find_window_size:
	$(CC) find_window_size.cpp -o find_window_size -lm -O3
# -I"./include" -L"./lib" -lfftw3-3
separation:
	$(CC) slope_scanning_median.cpp -o slope_scanning_median -lpthread -O3
