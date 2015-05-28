#!/bin/sh
wks="&workers=10"
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
    x--help)
	echo "Usage: $0 [local|lite|proof]" ; exit 0 ;; 
    *)
	echo "No or Unknown type ($1) specified" ; exit 1 ;;
esac
shift

url="${url}?mc${rv}${wks}#TE"

rm -rf $name
set -x
runTrain --class=HMTFMCTrain --name=$name --type=ESD --url="$url" $@ 

#
# EOF
#
