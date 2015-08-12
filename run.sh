#!/bin/bash
incoll=""
mc=$2
case x$1 in
    xlocal)
	runmode="local"
	case x$2 in
	    xpythia)
		incoll="${PWD}/pythia/input_files.dat"
		;;
	    xphojet)
		incoll="${PWD}/phojet/input_files.dat"
		;;
	    xdipsy)
		incoll="${PWD}/dipsy/input_files.dat"
		;;

	esac
	;;
    xlite)
	runmode="lite"
	case x$2 in
	    xpythia)
		incoll="${PWD}/pythia4/input_files.dat"
		;;
	    xphojet)
		incoll="${PWD}/phojet/input_files.dat"
		;;
	    xdipsy)
		incoll="${PWD}/dipsy/input_files.dat"
		;;
	esac
	;;
    xpod)
	runmode="pod"
	case x$2 in
	    xpythia)
		incoll="pythia"
		;;
	    xphojet)
		incoll="phojet"
		;;
	    xdipsy)
		incoll="pythia"
		;;
	esac
	;;
    x--help)
	echo "Usage: $0 [local|lite|proof]" ; exit 0
	;; 
    *)
	echo "No or Unknown type ($1) specified" ; exit 1
	;;
esac
shift


#rm -rf $runmode
#set -x
# runmode, nmax, debug, outfilename


outfile="hmtf_mc_mult_${mc}"
root -l -x "runProof.C(\"${runmode}\", -1, 0, \"${incoll}\", \"${outfile}\")"

#igprof -pp -z -o profiling.pp.gz root -l -x "\"runProof.C(\"${runmode}\", -1, 0, \"${incoll}\", \"${outfile}\")\""
#runTrain --class=HMTFMCMultEstTrain --name=$name --type=ESD --url="$url" $@ 

#
# EOF
#
