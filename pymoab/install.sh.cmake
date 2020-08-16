cd @CMAKE_BINARY_DIR@/pymoab

PYMOAB_INSTALL_PREFIX=@PYMOAB_INSTALL_PREFIX@

export PYTHONPATH="$PYMOAB_INSTALL_PREFIX:$PYTHONPATH"

if [[ ! -d $PYMOAB_INSTALL_PREFIX ]]; then
  mkdir -p @PYMOAB_INSTALL_PREFIX@
fi


if [ -x "$(command -v python)" ]; then
  python_exe=python
elif [ -x "$(command -v python3)" ]; then
  python_exe=python3
fi

${python_exe} setup.py install --prefix=@CMAKE_INSTALL_PREFIX@ --record @PYMOAB_INSTALL_PREFIX@/install_files.txt
