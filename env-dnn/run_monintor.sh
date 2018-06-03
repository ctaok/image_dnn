#!/bin/bash

if [ $# != 2 ]; then
	echo -e "Usage:\n\t$0 <result_dir> <auc_file>"
    exit;
fi

result_dir=$1
auc_log=$2
model_dir="tmp_model"
test_log="log.test"
config_test="test.conf"

mkdir -p $result_dir
mkdir -p $model_dir

while true
do
    for model in `ls -rt lr_model_epoch_* 2> /dev/null`
    do
        sleep 1s
		echo $model
        #result_file=`echo $model | awk -F "[_-.]" '{printf "result_epoch_%d_%d.txt\n", $4, $7}'`
        result_file=`echo $model | awk -F "[_]" '{printf "result_epoch_%d.txt\n", $4}'`
        if [ -z "`grep modelFile $config_test`" ]; then
            perl -p -i -e "s/^(jobType.*)/modelFile = $model/" $config_test
        else
            perl -p -i -e "s/^(modelFile.*)/modelFile = $model/" $config_test
        fi
        #perl -p -i -e "s/^(labelFile.*)/labelFile = $result_dir\/$result_file/" $config_test
        perl -p -i -e "s#^labelFile.*#labelFile = $result_dir/$result_file#" $config_test
        ./dnn_lr $config_test >> $test_log 2>&1
        ./cal_auc $result_dir/$result_file >> $auc_log
        mv $model $model_dir/
		rm lr_model_epoch_0.model*
        tmp_pid=`tail -10 $test_log | grep "total sentence" | awk '{print $3}'`
    done
    sleep 30s
done
