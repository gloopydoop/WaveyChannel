#!/bin/bash


function error_quit {
    echo -e "$@"
    echo
    echo -e 'Usage:'
    echo -e './compile_script --clean'
    echo -e '   To clean build direcrtory. Makenek will ask for cleaning 3rd-party libraries.'
    echo
    echo -e './compile_script --all'
    echo -e '   To compile the code.'
    exit 1
}

SETUP="channel"
# commands
cp_cmd=cp
mv_cmd=mv
rm_cmd=rm
# paths
SOURCE_ROOT="/scratch/nobis/Nek5000"
#modified files
NEKFILES=''
INCFILES=''
# makenek variables
F77="mpif77"
export F77
CC="mpicc"
export CC
PPLIST=""
export PPLIST
USR="frame.o mntrlog_block.o mntrlog.o mntrtmr_block.o mntrtmr.o rprm_block.o rprm.o io_tools_block.o io_tools.o chkpoint.o chkpt_mstp.o conht_tools.o IBM_tools.o parallfil.o"
USR+=" map2D.o stat.o stat_IO.o"
export USR

# arguments
args=("$@")
argsnr=$#

# check arguments
# parameters number check
if [ $[argsnr] -ne 1 ]; then
    error_quit 'Wrong arguments number!'
fi

for il in "$@"
do
case $il in
        --clean)
                ${SOURCE_ROOT}/bin/makenek clean
                shift
                ;;
        --all)
                ${SOURCE_ROOT}/bin/makenek ${CASE}
                shift
                ;;
        *) 
                error_quit 'Wrong argument.'
                ;;
esac
done
