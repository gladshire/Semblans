CC = g++
CFLAGS = '-Wl,-rpath,$$ORIGIN/../lib/' -g -pthread
LIBS = -L../lib/ -lboost_system -lboost_filesystem -ldl -lconfini
INCLUDE_PATH = -I../lib -I../include
OBJ_LINK = assemble.o sra.o sra_toolkit.o transcript.o ini_parse.o trinity_wrap.o print_info.o

../bin/assemble: $(OBJ_LINK)
	$(CC) $(CFLAGS) -o ../bin/assemble $(OBJ_LINK) $(LIBS)
assemble.o: assemble.cpp assemble.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c assemble.cpp $(LIBS)
sra_toolkit.o: sra_toolkit.cpp sra_toolkit.h ini_parse.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c sra_toolkit.cpp $(LIBS)
ini_parse.o: ini_parse.cpp ini_parse.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c ini_parse.cpp $(LIBS)
sra.o: sra.cpp sra.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c sra.cpp $(LIBS)
transcript.o: transcript.cpp transcript.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c transcript.cpp $(LIBS)
trinity_wrap.o: trinity_wrap.cpp trinity_wrap.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c trinity_wrap.cpp $(LIBS)
print_info.o: print_info.cpp print_info.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c print_info.cpp $(LIBS)
clean:
	rm assemble.o
	rm trinity_wrap.o
