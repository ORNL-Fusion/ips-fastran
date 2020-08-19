#!/bin/bash

# This is an entry script that will continuously copy everything using tar from the temporary directory while ips is running

echo $(date --iso-8601=seconds) - ips_continuous_copy started

"$@" &
ips_pid=$!

# get which config is used from ips command
config=$(python3 -c "import argparse, sys; p=argparse.ArgumentParser(); p.add_argument('-i','-j','--config','--simulation',type=str,action='store'); print(p.parse_known_args(sys.argv)[0].config)" "$@")

# get TMPXFS, RUN_ID, RANK and SUMMARY from config file
source <(grep = "$config" | sed 's/ //g')

CHECK_INTERVAL=${CHECK_INTERVAL:-10}
TAR_INTERVAL=${TAR_INTERVAL:-300}
FILES_TO_ARCHIVE="${FILES_TO_ARCHIVE:-SUMMARY run?????.log run?????.out run?????.config ips_?????.log ipslog.?????}"

run_tar () {
    echo $(date --iso-8601=seconds) - tar started
    file_name="$RUN_ID"_$(date +"%Y-%m-%dT%H.%M.%S").tar.gz
    pushd .
    cd "$TMPXFS/$RANK"
    echo tar caf "$SUMMARY"/"${file_name}" $FILES_TO_ARCHIVE
    time tar caf "$SUMMARY"/"${file_name}" $FILES_TO_ARCHIVE
    popd
    # Remove all *.tar.gz except most recent 2
    ls -t "$SUMMARY"/"$RUN_ID"_*.tar.gz | tail +3 | xargs rm -f
    echo $(date --iso-8601=seconds) - tar finished
}

last_tar_time=$(date +%s)

if [ "$TMPXFS" ];
then
    while ps -p ${ips_pid} > /dev/null ;
    do
        if [ $TAR_INTERVAL -ne "0" -a "$(($(date +%s)-$last_tar_time))" -gt $TAR_INTERVAL ]
        then
            run_tar
            last_tar_time=$(date +%s)
        fi
        sleep $CHECK_INTERVAL
    done
    echo $(date --iso-8601=seconds) - ips run finished
    # Create one final tar file since ips.py is now done
    run_tar
else
    wait ${ips_pid}
fi

echo $(date --iso-8601=seconds) - ips_continuous_copy finished

exit
