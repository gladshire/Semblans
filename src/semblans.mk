CC = g++ -std=c++11
CFLAGS = '-Wl,-rpath,$$ORIGIN/../lib/' -g -pthread
LIBS = -L../lib -lboost_system -lboost_filesystem -ldl -lconfini
INCLUDE_PATH = -I../lib/ -I../include
OBJ_LINK = semblans.o log.o print_info.o ini_parse.o sra.o transcript.o

../bin/semblans: semblans.o
	$(CC) $(OBJ_LINK) $(CFLAGS) -o ../bin/semblans $(LIBS)
semblans.o: semblans.cpp semblans.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c semblans.cpp $(LIBS)
log.o: log.cpp log.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c log.cpp $(LIBS)
print_info.o: print_info.cpp print_info.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c print_info.cpp $(LIBS)
ini_parse.o: ini_parse.cpp ini_parse.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c ini_parse.cpp $(LIBS)
sra.o: sra.cpp sra.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c sra.cpp $(LIBS)
transcript.o: transcript.cpp transcript.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c transcript.cpp $(LIBS)
clean:
	rm semblans.o
	rm print_info.o
	rm ini_parse.o
	rm sra.o
	rm transcript.o
