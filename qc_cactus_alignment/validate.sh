#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
	echo "Usage: ./validate.sh <input hal file>"
	exit -1
fi


hal=$1

# Check if the HAL database is valid
halValidate $hal


