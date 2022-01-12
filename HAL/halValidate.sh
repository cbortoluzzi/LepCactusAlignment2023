#!/bin/bash

# Author : @cb46


if [ -z $1 ]; then
  echo "Usage: ./halValidate.sh <cactus.alignment>"
  exit -1
fi


echo "halValidate $HAL"
halValidate $HAL

