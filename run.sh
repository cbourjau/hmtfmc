#!/bin/sh
wks="&workers=4"
rv="&pattern=galice.root&recursive"
case x$1 in
    xlocal)
	name="local"
	wks=
	url="local://${PWD}/input" ;;    
    xlite)
	name="lite"
	url="lite://${PWD}/input" ;;
    xproof|xcaf)
	name="caf"
	# rv="&root=v5-34-08"
	wks= #"&workers=10"
	rv="&aliphysics=vAN-20150323&root=v5-34-08&aliroot=v5-06-19&plugin"
	url="proof://alice-caf.cern.ch//default/ilakomov/TestDownScale1" ;;
    xalien|xgrid)
	name="alien"
	rv="&aliphysics=last,regular&pattern=*/galice.root&files=../inputs"
	dir="/alice/cern.ch/user/p/pwgpp_mc/2015/17_Week/TestMultiplicity/Test2/Test_1M_events_iter1"
	url="alien://${dir}"
	;;
    x--help)
	echo "Usage: $0 [local|lite|proof]" ; exit 0 ;; 
    *)
	echo "No or Unknown type ($1) specified" ; exit 1 ;;
esac
shift

url="${url}?mc${rv}${wks}#TE"

rm -rf $name
set -x
runTrain --class=HMTFMCMultEstTrain --name=$name --type=ESD --url="$url" $@ 

#
# EOF
#
