#!/bin/bash

# This is an entry script that will continuously copy everything using tar from the temporary directory while ips is running

"$@" &
ips_pid=$!

# get which config is used from ips command
config=$(python3 -c "import argparse, sys; p=argparse.ArgumentParser(); p.add_argument('-i','-j','--config','--simulation',type=str,action='store'); print(p.parse_known_args(sys.argv)[0].config)" "$@")

# get TMPXFS, RUN_ID and RANK from config file
source <(grep = "$config" | sed 's/ //g')

INTERVAL=10

HOSTNAME=$(hostname)

run_tar () {
    file_name="$RUN_ID"/$(date +"%Y-%m-%dT%H.%M.%S")_"$HOSTNAME".tar.gz
    pushd .
    cd "$TMPXFS/$RANK"
    echo tar caf "$OLDPWD"/"${file_name}" SUMMARY run?????.{log,out,config} ips_?????.log ipslog.?????
    time tar caf "$OLDPWD"/"${file_name}" SUMMARY run?????.{log,out,config} ips_?????.log ipslog.?????
    popd
    # Remove all *.tar.gz except most recent 2
    ls -t "$RUN_ID"/*_"$HOSTNAME".tar.gz | tail +3 | xargs rm -f
}

if [ "$TMPXFS" ];
then
    sleep $INTERVAL
    run_tar
    while ps -p ${ips_pid} > /dev/null ;
    do
        sleep $INTERVAL
        run_tar
    done
    # Create one final tar file since ips.py is now done
    run_tar
else
    wait ${ips_pid}
fi

exit
