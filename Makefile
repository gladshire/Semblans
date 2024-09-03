MAKE = make -C./src -f
#$(info $$HOME is [${HOME}])

all:
	mkdir -p bin
	$(MAKE) preprocess.mk
	$(MAKE) assemble.mk
	$(MAKE) postprocess.mk
	$(MAKE) semblans.mk

clean:
	$(MAKE) preprocess.mk clean
	$(MAKE) assemble.mk clean
	$(MAKE) postprocess.mk clean
	$(MAKE) semblans.mk clean
