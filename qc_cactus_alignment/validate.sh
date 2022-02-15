#!/bin/bash



# Author: @cb46



if [ -z $1 ]; then
	echo "Usage: ./validate.sh <cactus.alignment>"
	exit -1
fi


HAL=$1

# Check if the HAL database is valid; we can check this by simply running halValidate
halValidate $HAL
