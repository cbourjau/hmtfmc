#!/bin/bash
incoll=""
case x$1 in
    xlocal)
	name="local"
	case x$2 in
	    xpythia)
		incoll="${PWD}/pythia/input_files.dat"
		;;
	    xphojet)
		incoll="${PWD}/phojet/input_files.dat"
		;;
	esac
	;;
    xlite)
	name="lite"
	case x$2 in
	    xpythia)
		incoll="${PWD}/pythia/input_files.dat"
		;;
	    xphojet)
		incoll="${PWD}/phojet/input_files.dat"
		;;
	esac
	;;
    xpod)
	name="pod"
	case x$2 in
	    xpythia)
		incoll="pythia"
		;;
	    xphojet)
		incoll="phojet"
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

#rm -rf $name
#set -x
root -l -x "runProof.C(\"${name}\", -1, 0, \"${incoll}\")"
#runTrain --class=HMTFMCMultEstTrain --name=$name --type=ESD --url="$url" $@ 

#
# EOF
#
