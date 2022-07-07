#!/usr/bin/env bash

host=$(hostname | sed -e 's/\..\+//')
flavor="RelWithDebInfo"
intel_flags="-fp-model precise "

case $host in
    dhcp1-24)
	host="jarvis"
	;;
    gurney|tagir|licallo)
	host="licallo"
	;;
    compute-*)
	if grep -q atlantico /etc/hosts
	then
	    host=atlantico
	fi
	;;
    atlantico)
	echo "Atlantico frontal node not supported, please connect to a compute node." >&2
	exit 1
	;;
esac

usage() {
    echo "usage $(basename $0) [option]+"
    echo "Configure, through cmake, the build structure FWT3D."
    echo "Do NOT use it in source diretory (there are still legacy Makefiles in there."
    echo "with options:"
    echo "  -l: specify the hostname (default: based on hostname, here '$host')"
    echo "  -f <flavor>: one of 'Debug', 'Release', 'RelWithDebInfo' (default)."
    echo "  -g: same as '-f Debug'."
    echo "  -r: same as '-f Release'."
    echo "  -c: clean cmake cache before generating."
    echo "  -e: specify alternate cmake engine (default cmake)."
    echo "  -b <boost dir.>: Directory of boost)."
    echo "  -h: this help."
    
}

cmake_engine="cmake"

while getopts l:f:grce:b:h opt
do
    case $opt in 
	i)
	    profiling=$OPTARG
	    ;;
	l)
	    host=$OPTARG
	    ;;
	f)
	    flavor=$OPTARG
	    ;;
	g)
	    flavor="Debug"
	    ;;
	r)
	    flavor="Release"
	    ;;
	b)
	    boost_root=$OPTARG
	    if [ ! -d "$boost_root" ]
	    then
		echo "$boost_root does not exist." >&2
		exit 1
	    fi
	    ;;
	c)
	    if [ -f CMakeCache.txt ]
	    then
		echo removing cmake cache
		cp -v CMakeCache.txt CMakeCache.bak
		rm -v CMakeCache.txt
	    fi
	    ;;
	e)
	    cmake_override=$OPTARG
	    cmake_engine=$cmake_override
	    ;;
	h)
	    usage
	    exit 0
	    ;;
	\?)
	    echo "unknown option '$opt'" >&2
	    exit 1
    esac
done

src=$(readlink -f $(dirname $0))
if ! which pdflatex 2> /dev/null
then
    echo "Not pdflatex, safer to skip documentation processing"
    echo "Remove $src/SKIP_DOC if you want to retry."
    touch $src/SKIP_DOC
fi

if [ -x $(pwd)/$(basename $0) ]
then
    echo "We told you NOT to run this in distrib directory, go in a separate build directory and start again."
    exit
fi

cmd="$cmake_engine $src"

case $host in
    gurney|tagir|licallo)
	host="licallo"
	# Sanity check 
	if [[ "$CXX" != *icpc ]]
	then
	    echo "Could not find Intel C++ compiler in your environement" >&2
	    echo "Did you 'source /softs/env-intel15.0.2-impi5.0.3.sh' ?" >&2
	    echo "Don't forget to re run $0 with -c option ($0 -h for more information)." >&2
	    exit 1;
	fi
	
	impi_flags=$(mpicc -show | sed -e 's/icc //')
	echo "mpi flags: \"$impi_flags\""
	impi_compile_flags=$(echo $impi_flags | sed -e 's/\-L.*//')
	echo "mpi compile flags: \"$impi_compile_flags\""
	impi_link_flags=$impi_flags
	echo "mpi link flags: \"$impi_link_flags\""
	if [ -z "$boost_root" ]
	then
	    boost_root=/softs/boost-1.59.0-intel15.0.2-impi5.0.3
	fi
	echo "Using Boost located at \"$boost_root\"".
	cmd="$cmake_engine $prof_options -DCMAKE_CXX_FLAGS:STRING=\"$intel_flags\" -DMPI_COMPILE_FLAGS=\"$impi_compile_flags\" -DMPI_LINK_FLAGS=\"$impi_link_flags\" -DMPIEXEC=mpiexec.hydra -DMPIEXEC_NUMPROC_FLAG="-np" -DCMAKE_BUILD_TYPE=$flavor -DBOOST_ROOT=$boost_root $src"
	;;
    atlantico)
	export CXX=icpc
	export CC=icc
	if which mpicxx | grep -q -v /share/apps/openmpi-1.8.4-intel15
	then
	    echo "unsupported environment on $host." >&2
	    echo "Please do:" >&2
	    echo "    \$ module purge"
	    echo "    \$ module load intel-15 openmpi-1.8.4_ifb-intel-15"
	    echo "And retry" >&2
	    exit 1
	fi
	if [ -z "$cmake_override" ]
	then
	    #cmake_engine=/home/alain/install/cmake-3.1.3/bin/cmake
	    cmake_engine=cmake28
	fi
	if [ -z "$boost_root" ]
	then
	    boost_root=/share/apps/boost_intel15-opmi184
	    if [ "$flavor" == "Debug" ]
	    then
		boost_root="${boost_root}-dbg"
	    fi
	fi
	echo "Using Boost located at \"$boost_root\"".
	impi_compile_flags=$(mpicc -showme:compile)
	impi_link_flags=$(mpicc -showme:link)
	cmd="$cmake_engine $prof_options -DCMAKE_CXX_FLAGS:STRING=\"$intel_flags\" -DMPI_COMPILE_FLAGS='$impi_compile_flags' -DMPIEXEC=mpiexec -DMPIEXEC_NUMPROC_FLAG="-np"  -DMPI_LINK_FLAGS=\"$impi_link_flags\" -DCMAKE_BUILD_TYPE=$flavor -DBOOST_ROOT=$boost_root $src"
	;;
    thera)
	host="thera"
	mpicc=$(which mpicc)
	scpswenv=/user/alainm/seiscope
	if [[ "$mpicc" != "$scpswenv/"* ]]
	then
	    echo "Found mpicc='$mpicc' not in $scpswenv You do not have a supported compiler/mpi environment"
	    echo "Please source /user/alainm/seiscope/openmpi-1.8.4-gcc492/set-env.sh and try again"
	    echo "or /user/alainm/seiscope/openmpi-1.8.4-mt-gcc492/set-env.sh and try again"
	    exit 0
	fi
	if [ -z "$cmake_override" ]
	then
	    # no user override, adjust default to 2.8 or later
	    cmake_engine=cmake28
	fi
	mpi_compile_flags=$(mpif90 -showme:compile)
	mpi_link_flags=$(mpif90 -showme:link)
	if [ -z "$boost_root" ]
	then
	    boost_root=$scpswenv/boost
	fi
	cmd="$cmake_engine -DMPI_COMPILE_FLAGS=\"$mpi_compile_flags\" -DMPI_LINK_FLAGS=\"$mpi_link_flags\" -DCMAKE_BUILD_TYPE=$flavor $src"
	;;
    jarvis)
	if [ -z "$boost_root" ]
	then
	    . /opt/boost-1.58.0/env.sh
	    boost_root=/opt/boost-1.58.0
	fi
	echo "Using Boost located at \"$boost_root\"".
	cmd="$cmake_engine -DBOOST_ROOT=$boost_root -DCMAKE_BUILD_TYPE=$flavor $src .."
	;;
    *)
	echo "unsupported host '$host'" >&2
	exit 1
	;;
esac

echo "Will configure on host \"$host\" with:"
echo "    $cmd"

eval $cmd
