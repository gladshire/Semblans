CC = g++
CFLAGS = '-Wl,-rpath,$$ORIGIN/../lib/' -g -pthread
LIBS = -L../lib/ -lboost_system -lboost_filesystem -ldl -lconfini
BOOST_PATH = -I../lib/boost_1_80_0
INCLUDE_PATH = -I../lib -I../include
OBJ_LINK = assemble.o sra.o sra_toolkit.o transcript.o ini_parse.o trinity_wrap.o

../bin/assemble: $(OBJ_LINK)
	$(CC) $(CFLAGS) $(BOOST_PATH) -o ../bin/assemble $(OBJ_LINK) $(LIBS)
assemble.o: assemble.cpp assemble.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c assemble.cpp $(LIBS)
sra.o: sra.cpp sra.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c sra.cpp $(LIBS)
sra_toolkit.o: sra_toolkit.cpp sra_toolkit.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c sra_toolkit.cpp $(LIBS)
transcript.o: transcript.cpp transcript.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c transcript.cpp $(LIBS)
ini_parse.o: ini_parse.cpp ini_parse.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c ini_parse.cpp $(LIBS)
trinity_wrap.o: trinity_wrap.cpp trinity_wrap.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c trinity_wrap.cpp $(LIBS)
clean:
	rm assemble.o
	rm sra.o
	rm transcript.o
	rm ini_parse.o
	rm trinity_wrap.o
