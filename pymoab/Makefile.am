
AUTOMAKE_OPTIONS=foreign

EXTRA_DIST = $(srcdir)/cmake/cython_test.pyx $(srcdir)/pymoab/* $(srcdir)/tests/*
EXTRA_DIST += $(srcdir)/interface_notes.txt $(srcdir)/pyproject.toml $(srcdir)/requirements.txt

test_files = test_core.py test_hom.py test_scd.py test_range.py test_tag.py test_topo_util.py test_skinner.py

local_install_path = local_inst

all-local: local-install

local-install:
	(CFLAGS='-w -fPIC' PYTHONPATH=$(PWD)/$(local_install_path)/:$(PYTHONPATH) \
	 pip install -e $(PWD) --prefix $(PWD))

install-exec-local: all-local
	$(MKDIR_P) $(DESTDIR)$(pythondir)
	(PYTHONPATH=$(DESTDIR)$(pythondir):$(PYTHONPATH) pip install $(PWD) \
		--prefix=$(DESTDIR)$(prefix) )

check: all-local
	cd tests; \
	for i in $(test_files); \
	do \
		(@SLIBENVPATH@=$(abs_top_builddir)/src/.libs/:$(@SLIBENVPATH@) PYTHONPATH=$(abs_builddir)/../pymoab/:$(PYTHONPATH) $(PYTHON) $$i) || exit 1; \
	done

clean-local:
	-rm -rf build pymoab/*.cpp pymoab/__pycache__ pymoab/*.so pymoab/*.pyc tests/*.pyc pymoab.egg-info dist

uninstall:
	rm -rf $(DESTDIR)$(pythondir)/*

CLEANFILES = tests/quads_test.vtk tests/outfile.h5m tests/write_tag_test.h5m
