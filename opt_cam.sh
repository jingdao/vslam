#!/bin/bash

SIFT=/home/jd/Downloads/siftDemoV4/sift
BAG=$(rospack find lidar_mapping)/bag2img.py
process_image=false
use_g2o=false

if [ "$#" -ne "2" ] || ! [ -f $1 ] || ! [ -d $2 ]
then
	echo "./opt_cam.sh data.bag image_dir/"
	exit
fi

if ! [ -f $2/hector_pose.txt ]
then
	echo "Require $2/hector_pose.txt"
fi

if $process_image
then
	rm $2/*.ppm
	rm $2/*.pgm
	rm $2/*.key
	$BAG $1 $2 /camera/image_raw -s 2 -t 10 -r 0.5 -p $2/hector_pose.txt
	i=0
	while true
	do
		if [ -f $2/$i.ppm ]
		then
			convert $2/$i.ppm $2/$i.pgm
			$SIFT < $2/$i.pgm > $2/$i.key
			((i++))
		else
			break
		fi
	done
	./improc $2/key.match $2/*.key
fi
./check_epipole $2/pose_stamped.txt $2/key.match $2/valid.match
if $use_g2o
then
	./match_g2o $2/pose_stamped.txt $2/valid.match $2/map_point.txt
else
	./match_solver $2/pose_stamped.txt $2/valid.match $2/map_point.txt
fi
./viz_cam $2/pose_stamped.txt $2/map_point.txt
