
AUTOMAKE_OPTIONS=foreign

EXTRA_DIST = $(srcdir)/cmake/cython_test.pyx $(srcdir)/pymoab/* $(srcdir)/tests/*
EXTRA_DIST += $(srcdir)/interface_notes.txt $(srcdir)/pyproject.toml $(srcdir)/requirements.txt

test_files = test_core.py test_hom.py test_scd.py test_range.py test_tag.py test_topo_util.py test_skinner.py

local_install_path = $(subst $(prefix)/,,$(pythondir))

all-local: local-install

local-install:
	$(MKDIR_P) $(abs_top_builddir)/pymoab/$(local_install_path)
	(CFLAGS='-w -fPIC' PYTHONPATH=$(abs_top_builddir)/pymoab/$(local_install_path)/:$(PYTHONPATH) \
		$(PYTHON) -m pip install -e $(PWD) --prefix $(abs_top_builddir)/pymoab)

install-exec-local: all-local
	$(MKDIR_P) $(DESTDIR)$(pythondir)
	(PYTHONPATH=$(DESTDIR)$(pythondir):$(PYTHONPATH) $(PYTHON) -m pip install $(PWD) \
		--prefix=$(DESTDIR)$(prefix) )

check: all-local
	cd tests; \
	for i in $(test_files); \
	do \
		(@SLIBENVPATH@=$(abs_top_builddir)/src/.libs/:$(@SLIBENVPATH@) PYTHONPATH=$(abs_top_builddir)/pymoab:$(PYTHONPATH) $(PYTHON) $$i) || exit 1; \
	done

clean-local:
	-rm -rf build lib pymoab/*.cpp pymoab/__pycache__ pymoab/*.so pymoab/*.pyc tests/*.pyc pymoab.egg-info

uninstall:
	rm -rf $(DESTDIR)$(pythondir)/*

CLEANFILES = tests/quads_test.vtk tests/outfile.h5m tests/write_tag_test.h5m $(local_install_path)/*
