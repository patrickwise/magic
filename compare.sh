#!/bin/bash

# get the differences
paste $1 $2 | awk '{print $3-$10, $4-$11, $5-$12, $6-$13}'

#
