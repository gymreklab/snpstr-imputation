#!/bin/bash

set -e

AWS_ACCESS_KEY=$1
AWS_SECRET_KEY=$2
CHROM=$3
KEYNAME=$4

usage()
{
    BASE=$(basename -- "$0")
    echo "Launch amazon instance to run hipstr on a set of loci
Usage:
    $BASE <aws_access_key> <aws_secret_key> <chrosome> <keyname>
"
    exit 1
}
#test -z ${AWS_ACCESS_KEY} && usage
#test -z ${AWS_SECRET_KEY} && usage
#test -z ${CHROM} && usage
#test -z ${KEYNAME} && usage

# Instance details
SPOT_PRICE=0.15
INSTANCE_TYPE=m4.2xlarge
IMAGE_ID=ami-80861296

STARTUP_SCRIPT=$(cat run_from_aws.sh | \
    sed "s/\=\$1/\=${AWS_ACCESS_KEY}/" | sed "s~\=\$2~\=${AWS_SECRET_KEY}~" | \
    sed "s~\=\$3~\=${CHROM}~" | \
    sed "s/\=\$4/\=${KEYNAME}/")
STARTUP_SCRIPT_ENCODE="$(echo "${STARTUP_SCRIPT}" | gbase64 -w 0)"

LAUNCH_SPEC="{\"EbsOptimized\":true, \"ImageId\":\"${IMAGE_ID}\",\"Placement\":{\"AvailabilityZone\": \"us-east-1b\"},\"SecurityGroupIds\":[\"sg-5e914222\"], \"KeyName\":\"${KEYNAME}\",\"InstanceType\":\"${INSTANCE_TYPE}\", \"UserData\":\"${STARTUP_SCRIPT_ENCODE}\", \"BlockDeviceMappings\": [ {\"DeviceName\": \"/dev/sdf\",\"Ebs\": {\"VolumeSize\": 2,\"DeleteOnTermination\": true,\"VolumeType\": \"gp2\"}}]}"

aws ec2 request-spot-instances \
    --spot-price ${SPOT_PRICE} \
    --instance-count 1 \
    --type one-time \
    --launch-specification "${LAUNCH_SPEC}"
