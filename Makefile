all: bball project

bball: bball.cpp ppm.cpp
	g++ bball.cpp ppm.cpp libggfonts.a -Wall -Wextra -o bball -lX11 -lGL -lGLU -lm

project: project.cpp ppm.cpp
	g++ project.cpp ppm.cpp libggfonts.a -Wall -Wextra -o project -lX11 -lGL -lGLU -lm


clean:
	rm -f bball project
	rm -f *.o

