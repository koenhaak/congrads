#!/bin/bash

# Wrapper script for conmap.py and trendsurf.py; run without arguments for usage.
#
# Developed at DCCN (Donders Centre for Cognitive Neuroimaging), Donders Institute
# for Brain, Cognition and Behaviour. Radboud University, Nijmegen, The Netherlands
#
# Authors: KV Haak, AF Marquand, CF Beckmann.
#
# If you use CONGRADS in your research, please quote the following journal reference:
#
# Haak KV, Marquand AF, Beckmann CF (2018) Connectopic mapping with resting-state fMRI. 
# NeuroImage 170:83-94. 

usage () { echo -e "\nUsage:\n$0 -i <input> -r <roi> -m <mask> -o <out> [options]"
		   echo -e "\nCompulsory arguments:"
		   echo -e "\t-i\t\tinput file name [4D image]"
		   echo -e "\t-r\t\tregion-of-interest [3D binary image]"
		   echo -e "\t-m\t\tmask [3D binary image]"
		   echo -e "\t-o\t\tpath to output directory"
		   echo -e "\nOptional arguments:"
		   echo -e "\t-n <int>\tnumber of connectopic maps [default=1]"
		   echo -e "\t-s\t\tsave eta-squared matrix"
		   echo -e "\t-z\t\tnormalise output maps to range [0-1]"
		   echo -e "\t-p\t\tproject connectopic maps into mask"
		   echo -e "\t-f <int>\tfit spatial model of order <int>"	 
		   echo -e "\t-F <int>\tfit spatial model of order <int>"
		   echo -e "\t\t\tto pre-estimated connectopic map(s)"
		   echo -e "\t\t\tcompulsory arguments: -i -r -o\n"
         }

OPTIND=1

nev=1
fit=0
Fit=0

while getopts "i:r:m:o:n:f:F:pzs" opt; do
	case "$opt" in
	i)	fun+=("$OPTARG")
		;;
	r)	roi=$OPTARG
		;;
	m)	mas=$OPTARG
		;;
	o)	out=$OPTARG 
		;;
	n)	nev=$OPTARG
		;;
	s)	eta="--save_eta2"
		;;
	p)  pro="--project"
		;;
	z)	nor="--norm"
		;;
	f)  fit=$OPTARG
		;;
	F)  Fit=$OPTARG
	esac
done

shift $((OPTIND-1))

[ "$1" =  "--" ] && shift

type fslpython >/dev/null 2>&1 || { echo >&2 "Congrads requires fslpython but it's not installed. Aborting."; exit 1; }

if [ ! $Fit -gt 0 ]; then

	if [ ! "$fun" ] || [ ! "$roi" ] || [ ! "$mas" ] || [ ! "$out" ]; then
		usage
		exit 1
	fi

	infiles=()
	for fn in ${fun[@]}; do
		fn="$(readlink -f $fn)"
		if [ ! -f $fn ]; then
			echo "Can't find $fn. Aborting."
			exit 1
		else
			infiles+=($fn)
		fi
	done
	
	roi="$(readlink -f $roi)"
	if [ ! -f $roi ]; then
		echo "Can't find $roi. Aborting."
		exit 1
	fi

	mas="$(readlink -f $mas)"
	if [ ! -f $mas ]; then
		echo "Can't find $mas. Aborting."
		exit 1
	fi

	mkdir -p $out	
	out="$(readlink -f $out)"

	fslpython conmap.py -i ${infiles[@]} -r $roi -m $mas -o $out --nmaps $nev $eta $nor $pro

	if [ $fit -gt 0 ]; then
		cmap=${out}/$(basename "${roi}" .nii.gz).cmaps.nii.gz
		fslpython trendsurf.py -i $cmap -r $roi -b $fit -o $out
	fi

else

	if [ ! "$fun" ] || [ ! "$roi" ] || [ ! "$out" ]; then
		usage
		exit 1
	fi

	fun="$(readlink -f $fun)"
	if [ ! -f $fun ]; then
		echo "Can't find $fun. Aborting."
		exit 1
	fi

	roi="$(readlink -f $roi)"
	if [ ! -f $roi ]; then
		echo "Can't find $roi. Aborting."
		exit 1
	else
		echo "Region-of-interest: $roi"
	fi

	out="$(readlink -f $out)"
	echo "Writing results to: $out"
	mkdir -p $out

	fslpython trendsurf.py -i $fun -r $roi -o $out -b $Fit 

fi

# EOF


