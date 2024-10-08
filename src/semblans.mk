CC = g++ -std=c++11
CFLAGS = '-Wl,-rpath,$$ORIGIN/../lib/' -pthread
LIBS = -L../lib -lboost_system -lboost_filesystem -lboost_regex -ldl -lconfini
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
	rm -f semblans.o
	rm -f print_info.o
	rm -f ini_parse.o
	rm -f sra.o
	rm -f transcript.o
