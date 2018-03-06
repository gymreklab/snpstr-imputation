#!/bin/bash

set -e

AWS_ACCESS_KEY=$1
AWS_SECRET_KEY=$2
CHROM=$3
PART_START=$4
PART_END=$5
BATCH_SIZE=$6
KEYNAME=$7

usage()
{
    BASE=$(basename -- "$0")
    echo "Launch amazon instance to run hipstr on a set of loci
Usage:
    $BASE <aws_access_key> <aws_secret_key> <chrosome> <part_start> <part_end> <batch_size> <keyname>
"
    exit 1
}
#test -z ${AWS_ACCESS_KEY} && usage
#test -z ${AWS_SECRET_KEY} && usage
#test -z ${CHROM} && usage
#test -z ${PART_START} && usage
#test -z ${PART_END} && usage
#test -z ${BATCH_SIZE} && usage
#test -z ${KEYNAME} && usage

# Instance details
SPOT_PRICE=0.07
INSTANCE_TYPE=c4.xlarge
IMAGE_ID=ami-80861296

STARTUP_SCRIPT=$(cat run_from_aws.sh | \
    sed "s/\=\$1/\=${AWS_ACCESS_KEY}/" | sed "s~\=\$2~\=${AWS_SECRET_KEY}~" | \
    sed "s~\=\$3~\=${CHROM}~" | \
    sed "s~\=\$4~\=\"${PART_START}\"~" | \
    sed "s~\=\$5~\=${PART_END}~" | \
    sed "s~\=\$6~\=${BATCH_SIZE}~" | \
    sed "s/\=\$7/\=${KEYNAME}/")
STARTUP_SCRIPT_ENCODE="$(echo "${STARTUP_SCRIPT}" | gbase64 -w 0)"

LAUNCH_SPEC="{\"EbsOptimized\":true, \"ImageId\":\"${IMAGE_ID}\",\"Placement\":{\"AvailabilityZone\": \"us-east-1b\"},\"SecurityGroupIds\":[\"sg-5e914222\"], \"KeyName\":\"${KEYNAME}\",\"InstanceType\":\"${INSTANCE_TYPE}\", \"UserData\":\"${STARTUP_SCRIPT_ENCODE}\", \"BlockDeviceMappings\": [ {\"DeviceName\": \"/dev/sdf\",\"Ebs\": {\"VolumeSize\": 1000,\"DeleteOnTermination\": true,\"VolumeType\": \"gp2\"}}]}"

aws ec2 request-spot-instances \
    --spot-price ${SPOT_PRICE} \
    --instance-count 1 \
    --type one-time \
    --launch-specification "${LAUNCH_SPEC}"