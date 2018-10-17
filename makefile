CC = gcc -Wall -pedantic -I$(HOME)/include
LIBS = -lm -L$(HOME)/lib -lsofa_c

BIN = .

SRCS = skypath.c vecmat_utils.c skypath_lib.c
OBJS = skypath.o vecmat_utils.o skypath_lib.o

skypath.o: skypath.h vecmat_utils.h skypath_lib.h
vecmat_utils.o:

skypath: $(OBJS) $(SRCS)
	$(CC) -o $(BIN)/$@ $(OBJS) $(LIBS)

clean:
	rm -f *.o $(BIN)/skypath

.c.o:
	$(CC) -c $<
