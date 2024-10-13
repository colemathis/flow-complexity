#!/bin/bash

# Define the variable NUM
NUM=10

# Create logs directory if it doesn't exist
mkdir -p ./data/logs

# Loop from 1 to NUM and print each number
for (( i=1; i<=NUM; i++ ))
do
    echo julia launch.jl ./data/array.csv ${i} '>>' ./data/logs/launch.jl_${i}.log
done
