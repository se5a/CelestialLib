#!/bin/bash

echo "Starting"

#compile
if gcc celestial.c -o celestial; then
	echo "Success!";
	#clear screen and run
	clear
	./celestial
else 
	echo "Failure"; 
fi
