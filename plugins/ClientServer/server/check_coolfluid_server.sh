#!/bin/bash

count=$(ps aux | grep "app_server $1" | awk '$11!="grep"{print $11}' | wc -l)

echo -n $count

[ $count -eq 0 ] && exit 0
exit 1
