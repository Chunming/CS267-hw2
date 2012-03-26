CC = g++

OSTYPE = $(shell uname -s)
ifeq ($(OSTYPE), Darwin)
  INCLUDES = -I/Library/Frameworks/SDL.framework/Versions/A/Headers/
  LIBS = -lSDLmain -lSDL
  FRAMEWORKS = -framework SDL -framework OpenGL -framework Cocoa
else
  INCLUDES = -I/usr/include/SDL
  LIBS = -lSDLmain -lSDL -lGL -lGLU
  FRAMEWORKS =
endif

visualize: visualize.cpp
	$(CC) $(INCLUDES) $(FRAMEWORKS) $(LIBS) -o $@ $< 

clean:
	rm visualize