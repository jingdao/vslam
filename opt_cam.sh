#!/bin/bash

SIFT=/home/jd/Downloads/siftDemoV4/sift
BAG=$(rospack find lidar_mapping)/bag2img.py

if [ "$#" -ne 2 ] || ! [ -f $1 ] || ! [ -d $2 ]
then
	echo "./opt_cam.sh data.bag image_dir/"
	exit
fi

if ! [ -f $2/hector_pose.txt ]
then
	echo "Require $2/hector_pose.txt"
fi

rm $2/*.ppm
rm $2/*.pgm
rm $2/*.key
$BAG $1 $2 /camera/image_raw -s 5 -t 5 -r 3 -p $2/hector_pose.txt
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
./match_g2o $2/pose_stamped.txt $2/key.match $2/map_point.txt
./viz_cam $2/pose_stamped.txt $2/map_point.txt
