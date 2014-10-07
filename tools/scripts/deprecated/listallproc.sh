#!/bin/bash
# @author Tiago Quintino

SCRIPT_NAME=$(basename $0)

on_error () {
	echo $SCRIPT_NAME - list all the processes with the given name
	echo usage: $SCRIPT_NAME process_name
	exit 1
}

if [ -z $1 ]; then
	on_error
fi

ps u -C "$1" | grep $USER
