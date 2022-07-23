#!/bin/bash
# Usage:
#
# Use clang-format to enforce uniform coding style throughout the entire repository or per-file basis
#   config/dev/run_clang_format.sh <source_file_name>
# Or use the following command to reformat the entire repository
#   git ls-tree --name-only HEAD -r |  grep "cpp\|hpp" | xargs config/dev/run_clang_format.sh
#   git ls-tree --name-only HEAD -r |  grep "\.c$" | xargs config/dev/run_clang_format.sh
CLANGFORMAT_EXE=`which clang-format`
SEDWITHOPTIONS="sed -i.bak"

if [ "$#" -eq 0 ]; then
    echo "Illegal number of parameters"
    echo "Usage: run_clang_format.sh <list_of_source_files>"
    exit 1
fi

allfiles=$*

function process_source()
{
	srcfile=$1

	echo "Processing $srcfile ..."
	eval `${CLANGFORMAT_EXE} -i -style=file ${srcfile}`
	eval `${SEDWITHOPTIONS} -n -e 'H;${x;s/;\n *MB_CHK/;MB_CHK/g;p;}' ${srcfile}`
	eval `${SEDWITHOPTIONS} -n -e 'H;${x;s/;\n *RR;/;RR;/g;p;}' ${srcfile}`
	eval `${SEDWITHOPTIONS} -n -e 'H;${x;s/;\n *CHECK_ERR/;CHECK_ERR/g;p;}' ${srcfile}`
	eval `${SEDWITHOPTIONS} -n -e 'H;${x;s/;\n *CHK_MPI_ERR/;CHK_MPI_ERR/g;p;}' ${srcfile}`
	eval `${SEDWITHOPTIONS} -n -e 'H;${x;s/;\n *CHKERR/;CHKERR/g;p;}' ${srcfile}`
	eval `${SEDWITHOPTIONS} -n -e 'H;${x;s/;\n *CHK_ERR/;CHK_ERR/g;p;}' ${srcfile}`
	eval `${SEDWITHOPTIONS} -n -e 'H;${x;s/;\n *ERRORR/;ERRORR/g;p;}' ${srcfile}`
	eval `${SEDWITHOPTIONS} -n -e 'H;${x;s/;\n *MSQ_CHKERR/;MSQ_CHKERR/g;p;}' ${srcfile}`
	eval `${SEDWITHOPTIONS} -n -e 'H;${x;s/;\n *MSQ_ERRRTN/;MSQ_ERRRTN/g;p;}' ${srcfile}`
	eval `${SEDWITHOPTIONS} '/./,$!d' ${srcfile}`
	eval `rm ${srcfile}.bak`

}

# Process all given user input source files
for srcfile in ${allfiles}
do
	process_source ${srcfile}
done

exit 0
