CC = g++
CFLAGS = '-Wl,-rpath,$$ORIGIN/../lib/' -g -pthread
LIBS = -L../lib -lboost_system -lboost_filesystem -lboost_iostreams -ldl -lconfini -lcurl
INCLUDE_PATH = -I../lib/ -I../include
OBJ_LINK = postprocess.o sra.o sra_toolkit.o ini_parse.o transcript.o seq.o seq_hash.o ncbi_blast.o diamond.o rem_chimera.o salmon_wrap.o corset_wrap.o filter_corset.o transdecoder_wrap.o print_info.o thread_pool.o ips_client.o ips_job_man.o


../bin/postprocess: $(OBJ_LINK)
	$(CC) $(OBJ_LINK) $(CFLAGS) -o ../bin/postprocess $(LIBS)
postprocess.o: postprocess.cpp postprocess.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c postprocess.cpp $(LIBS)
sra.o: sra.cpp sra.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c sra.cpp $(LIBS)
sra_toolkit.o: sra_toolkit.cpp sra_toolkit.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c sra_toolkit.cpp $(LIBS)
ini_parse.o: ini_parse.cpp ini_parse.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c ini_parse.cpp $(LIBS)
transcript.o: transcript.cpp transcript.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c transcript.cpp $(LIBS)
seq.o: seq.cpp seq.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c seq.cpp $(LIBS)
seq_hash.o: seq_hash.cpp seq_hash.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c seq_hash.cpp $(LIBS)
ncbi_blast.o: ncbi_blast.cpp ncbi_blast.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c ncbi_blast.cpp $(LIBS)
diamond.o: diamond.cpp diamond.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c diamond.cpp $(LIBS)
rem_chimera.o: rem_chimera.cpp rem_chimera.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c rem_chimera.cpp $(LIBS)
salmon_wrap.o: salmon_wrap.cpp salmon_wrap.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c salmon_wrap.cpp $(LIBS)
corset_wrap.o: corset_wrap.cpp corset_wrap.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c corset_wrap.cpp $(LIBS)
filter_corset.o: filter_corset.cpp filter_corset.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c filter_corset.cpp $(LIBS)
transdecoder_wrap.o: transdecoder_wrap.cpp transdecoder_wrap.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c transdecoder_wrap.cpp $(LIBS)
print_info.o: print_info.cpp print_info.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c print_info.cpp $(LIBS)
thread_pool.o: thread_pool.cpp thread_pool.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c thread_pool.cpp $(LIBS)
ips_client.o: ips_client.cpp ips_client.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c ips_client.cpp $(LIBS)
ips_job_man.o: ips_job_man.cpp ips_job_man.h
	$(CC) $(CFLAGS) $(INCLUDE_PATH) -c ips_job_man.cpp $(LIBS)
clean:
	rm $(OBJ_LINK)
