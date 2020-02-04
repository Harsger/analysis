#!/usr/bin/env bash
#

source /etc/profile.d/modules.sh

module load root

command=""

for input; do
    command="$command ${input}"
done

echo " @ $(hostname) : $command"

date +%Y-%m-%dT%H:%M:%S
$command
date +%Y-%m-%dT%H:%M:%S

