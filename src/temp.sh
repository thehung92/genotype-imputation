#!/bin/bash
# test parallel
sleep 5 &
sleep 10 &

# wait for all background to finish
wait `jobs -p` && echo "execute after background process finished" 
echo "execute now"
