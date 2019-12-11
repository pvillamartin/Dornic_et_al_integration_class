main: example_main.cpp
	g++ -O3 example_main.cpp `pkg-config --cflags --libs gsl` -o example_main
