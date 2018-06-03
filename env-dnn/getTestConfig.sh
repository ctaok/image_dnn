#!/bin/bash

if [ $# != 2 ]; then
    echo -e "Usage:\n\t$0 <trainconfig> <testconfig>"
    exit;
fi

TRAIN_FILE=$1
TEST_FILE=$2

if [ -f $TEST_FILE ]; then
	rm $TEST_FILE
fi
cp $TRAIN_FILE $TEST_FILE
if [ -z "`grep modelFile $TEST_FILE`" ]; then
    perl -p -i -e "s/^(jobType.*)/\1\nmodelFile = xxx/" $TEST_FILE
else
    perl -p -i -e "s/^(modelFile.*)/#\1\nmodelFile = xxx/" $TEST_FILE
fi
perl -p -i -e "s/^epoch.*/epoch = 1/" $TEST_FILE
perl -p -i -e "s/^(keepHistory.*)/\1\nlabelFile = xxx/" $TEST_FILE
perl -p -i -e "s/testSentBp/predictSentBp/" $TEST_FILE
perl -p -i -e "s/^dataFileList.*/dataFileList = list\/test.lst/" $TEST_FILE

