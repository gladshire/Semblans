CC = g++ -std=c++11
CFLAGS = '-Wl,-rpath,$$ORIGIN/../lib/' -pthread -g
#LIBS = -L../lib/ -lboost_system -lboost_filesystem -lboost_iostreams -ldl -lconfini
LIBS = -L../lib -Wl,-Bstatic -lboost_atomic -lboost_system -lboost_filesystem -lboost_iostreams -Wl,-Bdynamic -ldl -lconfini
INCLUDE_PATH = -I../lib -I../include
OBJ_LINK = assemble.o log.o sra.o sra_toolkit.o transcript.o ini_parse.o trinity_wrap.o star_wrap.o print_info.o seq.o llist.o seq_hash.o

../bin/assemble: $(OBJ_LINK)
	$(CC) $(CFLAGS) -o ../bin/assemble $(OBJ_LINK) $(LIBS) -lz
assemble.o: assemble.cpp assemble.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c assemble.cpp $(LIBS)
log.o: log.cpp log.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c log.cpp $(LIBS)
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
star_wrap.o: star_wrap.cpp star_wrap.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c star_wrap.cpp
print_info.o: print_info.cpp print_info.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c print_info.cpp $(LIBS)
seq.o: seq.cpp seq.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c seq.cpp $(LIBS)
llist.o: llist.cpp llist.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c llist.cpp $(LIBS)
seq_hash.o: seq_hash.cpp seq_hash.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c seq_hash.cpp $(LIBS)
clean:
	rm assemble.o
	rm trinity_wrap.o
