CC = g++
CFLAGS = -Wl,-rpath,../lib/ -g
LIBS = -L../lib/ -lboost_system -lboost_filesystem -lconfini
BOOST_PATH = -I../lib/boost_1_80_0
INCLUDE_PATH = -I../lib/ -I../include
OBJ_LINK = preprocess.o sra.o ini_parse.o sra_toolkit.o fastqc_wrap.o rcorr_wrap.o print_info.o


preprocess: $(OBJ_LINK)
	$(CC) $(CFLAGS) $(BOOST_PATH) -o preprocess $(OBJ_LINK) $(LIBS)
preprocess.o: preprocess.cpp preprocess.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c preprocess.cpp $(LIBS)
sra_toolkit.o: sra_toolkit.cpp sra_toolkit.h ini_parse.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c sra_toolkit.cpp $(LIBS)
ini_parse.o: ini_parse.cpp ini_parse.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c ini_parse.cpp $(LIBS)
sra.o: sra.cpp sra.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c sra.cpp $(LIBS)
fastqc_wrap.o: fastqc_wrap.cpp fastqc_wrap.h preprocess.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c fastqc_wrap.cpp $(LIBS)
rcorr_wrap.o: rcorr_wrap.cpp fastqc_wrap.h preprocess.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c rcorr_wrap.cpp $(LIBS)
print_info.o: print_info.cpp print_info.h
	$(CC) $(CFLAGS) $(BOOST_PATH) $(INCLUDE_PATH) -c print_info.cpp $(LIBS)
clean:
	rm ini_parse.o
	rm sra_toolkit.o
	rm sra.o
	rm fastqc_wrap.o
	rm rcorr_wrap.o
	rm preprocess.o
	rm print_info.o
