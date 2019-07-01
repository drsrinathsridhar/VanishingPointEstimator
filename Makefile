INCDIR = -I. -I/usr/local/include/opencv
DBG    = -g
CC    = gcc
CFLAGS = $(DBG) $(INCDIR) -O3 # Maximum compiler optimization
LINK   = -lm -lopencv_core -lopencv_imgproc -lopencv_calib3d -lopencv_video -lopencv_features2d

all: findVanPoints

findVanPoints: sample_usage.o VanPoints.o
	$(CC) $(CFLAGS) $(LINK) sample_usage.o VanPoints.o -o findVanPoints

sample_usage.o: sample_usage.cpp
	$(CC) $(CFLAGS) -c sample_usage.cpp

VanPoints.o: VanPoints.cpp VanPoints.h
	$(CC) $(CFLAGS) -c VanPoints.cpp

.PHONY: clean
clean:
	rm -f *.o
