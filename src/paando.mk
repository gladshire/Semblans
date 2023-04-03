CC = g++
CFLAGS = '-Wl,-rpath,$$ORIGIN/../lib/' -pthread
LIBS = -L../lib -lboost_system -lboost_filesystem -ldl -lconfini
INCLUDE_PATH = -I../lib/ -I../include
OBJ_LINK = paando.o print_info.o ini_parse.o sra.o

../bin/Paando: paando.o
	$(CC) $(OBJ_LINK) $(CFLAGS) -o ../bin/paando $(LIBS)
paando.o: paando.cpp paando.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c paando.cpp $(LIBS)
print_info.o: print_info.cpp print_info.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c print_info.cpp $(LIBS)
ini_parse.o: ini_parse.cpp ini_parse.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c ini_parse.cpp $(LIBS)
sra.o: sra.cpp sra.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c sra.cpp $(LIBS)
clean:
	rm paando.o
	rm print_info.o
	rm ini_parse.o
	rm sra.o
