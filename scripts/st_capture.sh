#!/bin/bash

# Test for device
if [ -z $1 ]; then
    device=/dev/video0
else
    device=$1
fi

# Go to shared memory
mkdir -p /dev/shm
cd /dev/shm

# Remove old files
rm img*.pgm

# Start capture
ffmpeg -f video4linux2 -i $device -s 720x576 -r 25 img%06d.pgm

