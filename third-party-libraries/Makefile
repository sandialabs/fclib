



CC = gcc
CFLAGS = -g -O2



NETCDF = netcdf-3.6.2
EXODUS = exodusii-4.71
LIBXML2 = libxml2-2.6.31

BUNDLE_DIR = $(PWD)/bundle

THIS_DIR = dependency_libs
TAR_NAME = $(THIS_DIR).tar

ALL_FILES = 	$(THIS_DIR)/src/$(NETCDF).tar.gz \
		$(THIS_DIR)/src/$(EXODUS).tar.gz \
		$(THIS_DIR)/src/$(LIBXML2).tar.gz \
		$(THIS_DIR)/patches/patch_netcdf.patch \
		$(THIS_DIR)/patches/patch_$(EXODUS).patch \
		$(THIS_DIR)/Makefile \
		$(THIS_DIR)/README.txt




all:	netcdf exodus
everything: netcdf exodus libxml2

#==[ NetCDF Stuff ]============================================================
$(NETCDF):
	rm -R -f $(NETCDF);
	gunzip -c src/$(NETCDF).tar.gz | tar -xf -;
	patch -p0 < patches/patch_netcdf.patch

$(BUNDLE_DIR)/lib/libnetcdf.a: $(NETCDF)
	cd $(NETCDF);				\
	CC=$(CC) ./configure --prefix=$(BUNDLE_DIR) --disable-f77 --disable-f90; \
	make install;

# Note: NetCDF's --enable-c-only option prevents some stuff with globalvars
#  from being built. Use the disable commands specifically for fortran


netcdf: $(BUNDLE_DIR)/lib/libnetcdf.a


#==[ Exodus Stuff ]============================================================
$(EXODUS):
	rm -R -f $(EXODUS);
	gunzip -c src/$(EXODUS).tar.gz | tar -xf - ;
	patch -p0 < patches/patch_$(EXODUS).patch;


$(BUNDLE_DIR)/lib/libexoIIv2c.a: $(EXODUS)
	cd $(EXODUS); 						\
	make -f Makefile.standalone BUNDLE_DIR=$(BUNDLE_DIR) CC=$(CC) CCOPTIONS="$(CFLAGS)" libexoIIv2c.a;	\
	cp libexoIIv2c.a $(BUNDLE_DIR)/lib;			\
	cp cbind/include/*.h $(BUNDLE_DIR)/include;


exodus: $(BUNDLE_DIR)/lib/libexoIIv2c.a



release:
	cd ..;	\
	rm -f $(TAR_NAME).gz; \
	tar -cf $(TAR_NAME) $(ALL_FILES); \
	gzip $(TAR_NAME); 

#==[ libxml2 ]=================================================================
# This lib is usually installed on a machine, but if you nee to recompile for
# icc, here's something to start with.
$(LIBXML2):
	rm -R -f $(LIBXML2);
	gunzip -c src/$(LIBXML2).tar.gz | tar -xf -;

$(BUNDLE_DIR)/lib/libxml2.a: $(LIBXML2)
	cd $(LIBXML2); \
	CC=$(CC) ./configure --prefix=$(BUNDLE_DIR); \
	make install;

libxml2: $(BUNDLE_DIR)/lib/libxml2.a



#==[ Cleanup ]=================================================================

clean_netcdf:
	rm -f $(BUNDLE_DIR)/lib/libnetcdf.a;
	if [ -d $(NETCDF) ]; then 			\
		cd $(NETCDF);				\
		make clean; 				\
	fi;

clean_exodus:
	rm -f $(BUNDLE_DIR)/lib/libexoIIv2c.a $(EXODUS)/libexoIIv2c.a;
	if [ -d $(EXODUS) ]; then 			\
		cd $(EXODUS); 				\
		make -f Makefile.standalone clean; 	\
	fi;

clean_libxml2:
	rm -f $(BUNDLE_DIR)/lib/libxml2.a;
	if [ -d $(LIBXML2) ]; then 			\
		cd $(LIBXML2); 				\
		make clean; 				\
	fi;


clean: clean_netcdf clean_exodus clean_libxml2




distro: clean_netcdf clean_exodus clean_libxml2
	rm -R -f $(BUNDLE_DIR) $(NETCDF) $(EXODUS) $(LIBXML2);
