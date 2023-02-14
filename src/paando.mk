CC = g++
CFLAGS = '-Wl,-rpath,$$ORIGIN/../lib/' -g -pthread
LIBS = -L../lib -lboost_system -lboost_filesystem -ldl -lconfini
BOOST_PATH = -I../lib/boost_1.80_0
INCLUDE_PATH = -I../lib/ -I../include
OBJ_LINK = paando.o

../bin/Paando: paando.o
	$(CC) $(OBJ_LINK) $(CFLAGS) $(BOOST_PATH) -o ../bin/paando $(LIBS)
paando.o: paando.cpp paando.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c paando.cpp $(LIBS)
clean:
	rm paando.o
