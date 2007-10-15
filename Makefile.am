# Copyright (C) 2000 Andreas Heger <heger@ebi.ac.uk>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


EXTRA_DIST = configure
SUBDIRS = m4 doc src tests bench

cvs-dist:
	pkg=`echo "@PACKAGE@" | tr a-z A-Z`; \
	ver=`echo "@VERSION@" | sed 's/\./_/g'`; \
	tag="$$pkg-$$ver"; \
	echo tag=$$tag; \
	if cvs -n log -h README| egrep -e $$tag > /dev/null; then \
		echo "VERSION not new; not releasing" 1>&2; \
		exit 1; \
	else :; \
	fi; \
	cvs tag -c $$tag
	$(MAKE) dist
