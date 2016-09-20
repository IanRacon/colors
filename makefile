app: main.cpp
	g++ -std=c++11 -O3 -ffast-math $^ -o app
clear:
	rm app
