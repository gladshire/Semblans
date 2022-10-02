CC = g++
CFLAGS = -Wl,-rpath,../lib/ -g
LIBS = -L../lib/ -lboost_system -lboost_filesystem -lconfini
BOOST_PATH = -I../lib/boost_1_80_0
INCLUDE_PATH = -I../lib/ -I../include

preprocess: ini_parse.o sra_toolkit.o sra.o fastqc_wrap.o
	$(CC) $(CFLAGS) $(BOOST_PATH) -o preprocess ini_parse.o sra_toolkit.o sra.o fastqc_wrap.o $(LIBS)
sra_toolkit.o: sra_toolkit.cpp sra_toolkit.h ini_parse.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c sra_toolkit.cpp $(LIBS)
ini_parse.o: ini_parse.cpp ini_parse.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c ini_parse.cpp $(LIBS)
sra.o: sra.cpp sra.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c sra.cpp $(LIBS)
fastqc_wrap.o: fastqc_wrap.cpp
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c fastqc_wrap.cpp $(LIBS)
clean:
	rm ini_parse.o
	rm sra_toolkit.o
	rm sra.o
	rm fastqc_wrap.o
