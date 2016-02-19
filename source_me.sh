export LIBRARY_PATH=`pwd`/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
export LD_INCLUDE_PATH=`pwd`/lib:$LD_INCLUDE_PATH
export C_INCLUDE_PATH=`pwd`/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=`pwd`/include:$CPLUS_INCLUDE_PATH
export INCLUDE_PATH=`pwd`/include:$INCLUDE_PATH
export PATH=`pwd`/bin:$PATH
if [[ "x$CC" == "x" ]]
then
    export CC=$(which gcc)
fi
    
if [[ "x$CXX" == "x" ]]
then
    export CXX=$(which g++)
fi
