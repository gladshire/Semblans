CC = g++ -std=c++11
CFLAGS = '-Wl,-rpath,$$ORIGIN/../lib/' -pthread

#LIBS = -L../lib/ -lboost_atomic -lboost_system -lboost_filesystem -lboost_iostreams -ldl -lconfini
LIBS = -L../lib/ -Wl,-Bstatic -lboost_atomic -lboost_system -lboost_filesystem -lboost_regex -lboost_iostreams -Wl,-Bdynamic -ldl -lconfini
INCLUDE_PATH = -I../lib/ -I../include
OBJ_LINK = preprocess.o sra_toolkit.o ini_parse.o log.o sra.o transcript.o fastqc_wrap.o rcorr_wrap.o rem_unfixable.o trimm_wrap.o kraken_wrap.o rem_overrep.o print_info.o seq.o llist.o seq_hash.o


../bin/preprocess: $(OBJ_LINK)
	$(CC) $(CFLAGS) -o ../bin/preprocess $(OBJ_LINK) $(LIBS) -lz
preprocess.o: preprocess.cpp preprocess.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c preprocess.cpp $(LIBS)
sra_toolkit.o: sra_toolkit.cpp sra_toolkit.h ini_parse.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c sra_toolkit.cpp $(LIBS)
ini_parse.o: ini_parse.cpp ini_parse.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c ini_parse.cpp $(LIBS)
log.o: log.cpp log.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c log.cpp $(LIBS)
sra.o: sra.cpp sra.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c sra.cpp $(LIBS)
transcript.o: transcript.cpp transcript.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c transcript.cpp $(LIBS)
fastqc_wrap.o: fastqc_wrap.cpp fastqc_wrap.h preprocess.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c fastqc_wrap.cpp $(LIBS)
rcorr_wrap.o: rcorr_wrap.cpp fastqc_wrap.h preprocess.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c rcorr_wrap.cpp $(LIBS)
rem_unfixable.o: rem_unfixable.cpp rem_unfixable.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c rem_unfixable.cpp $(LIBS)
trimm_wrap.o: trimm_wrap.cpp trimm_wrap.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c trimm_wrap.cpp $(LIBS)
kraken_wrap.o: kraken_wrap.cpp kraken_wrap.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c kraken_wrap.cpp $(LIBS)
rem_overrep.o: rem_overrep.cpp rem_overrep.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c rem_overrep.cpp $(LIBS)
print_info.o: print_info.cpp print_info.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c print_info.cpp $(LIBS)
seq.o: seq.cpp seq.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c seq.cpp $(LIBS)
llist.o: llist.cpp llist.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c llist.cpp $(LIBS)
seq_hash.o: seq_hash.cpp seq_hash.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c seq_hash.cpp $(LIBS)
clean:
	rm preprocess.o
	rm fastqc_wrap.o
	rm rcorr_wrap.o
	rm rem_unfixable.o
	rm trimm_wrap.o
	rm kraken_wrap.o
	rm rem_overrep.o
	rm print_info.o
