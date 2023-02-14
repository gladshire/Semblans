MAKE = make -C./src -f
$(info $$HOME is [${HOME}])
all:
	$(MAKE) preprocess.mk
	$(MAKE) assemble.mk
	$(MAKE) postprocess.mk
	$(MAKE) paando.mk
