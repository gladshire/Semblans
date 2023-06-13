MAKE = make -C./src -f
#$(info $$HOME is [${HOME}])
all:
	$(MAKE) preprocess.mk
	#$(MAKE) preprocess.mk clean
	$(MAKE) assemble.mk
	#$(MAKE) assemble.mk clean
	$(MAKE) postprocess.mk
	#$(MAKE) postprocess.mk clean
	$(MAKE) semblans.mk
	#$(MAKE) semblans.mk clean
