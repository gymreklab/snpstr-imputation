#!/bin/bash

AWS_ACCESS_KEY=$1
AWS_SECRET_KEY=$2
CHROM=$3


HOMEDIR=/root/
AWS_DIR=${HOMEDIR}/.aws
AWS_CONFIG_FILE=${AWS_DIR}/config
AWS_CRED_FILE=${AWS_DIR}/credentials


OUTBUCKET=s3://gatk-results

usage()
{
    BASE=$(basename -- "$0")
    echo "Run HipSTR
Usage:
    $BASE <aws access key> <aws secret key> <chrom> <part_start> <part_end> <batch_size>
       - aws access key and aws secret keys are for AWS configuration
       - chrom you are calling from
Does the following:
1. Set up AWS configuration
2. Download necessary files
3. Create jobs
4. Run those jobs
5. Upload results to S3 bucket
6. Terminate
"
    terminate
    exit 1
}

terminate() {

    INSTANCE_ID=$(curl http://169.254.169.254/latest/meta-data/instance-id)
    # Get log
    aws s3 cp --output table /var/log/cloud-init-output.log ${OUTBUCKET}/l1o_log/${CHROM}.${INSTANCE_ID}.log
    # Terminate instance
    echo "Terminating instance ${INSTANCE_ID}"
    aws ec2 terminate-instances --output table --instance-ids ${INSTANCE_ID}
    exit 1 # shouldn't happen
}

#test -z ${AWS_ACCESS_KEY} && usage
#test -z ${AWS_SECRET_KEY} && usage
#test -z ${CHROM} && usage


die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@" >&2
    terminate
    exit 1
}

# Install things
sudo apt-get update || die "Could not update"
sudo apt-get -y install awscli || die "Could not install aws"
sudo apt-get -y install git || die "Could not install git"
sudo apt-get -y install make gcc libz-dev libncurses5-dev libbz2-dev liblzma-dev libcurl3-dev libssl-dev autoconf g++ python2.7 parallel python-pip openjdk-8-jre || die "Could not install devtools"
sudo pip install cyvcf2

cd ${HOMEDIR}
git clone https://github.com/samtools/htslib
cd htslib
git log --pretty=format:'%h' -n 1
autoheader
autoconf
./configure --enable-libcurl
make
sudo make install

cd ${HOMEDIR}
wget https://github.com/samtools/bcftools/releases/download/1.7/bcftools-1.7.tar.bz2
tar -xvjf bcftools-1.7.tar.bz2
cd bcftools-1.7
make
sudo make install

# Set up AWS credentials
echo "Setting up AWS credentials in ${AWS_DIR}"
mkdir -p ${AWS_DIR} || die "Could not create AWS dir"
echo "[default]" > ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CONFIG_FILE}"
echo "output = table" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CONFIG_FILE}"
echo "region = us-east-1" >> ${AWS_CONFIG_FILE}  || die "Could not write to ${AWS_CONFIG_FILE}"
echo "aws_secret_access_key = ${AWS_SECRET_KEY}" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_access_key_id = ${AWS_ACCESS_KEY}" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "[default]" > ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_secret_access_key = ${AWS_SECRET_KEY}" >> ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_access_key_id = ${AWS_ACCESS_KEY}" >> ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"

export AWS_ACCESS_KEY_ID=${AWS_ACCESS_KEY}
export AWS_SECRET_ACCESS_KEY=${AWS_SECRET_KEY}

# Set ulimit
echo "fs.file-max = 13107" | sudo tee -a /etc/sysctl.conf
sudo sysctl -p
echo "* soft     nproc          13107" | sudo tee -a /etc/security/limits.conf
echo "* hard     nproc          13107" | sudo tee -a /etc/security/limits.conf
echo "* soft     nofile         13107" | sudo tee -a /etc/security/limits.conf
echo "* hard     nofile         13107" | sudo tee -a /etc/security/limits.conf
echo "root soft     nproc          13107" | sudo tee -a /etc/security/limits.conf
echo "root hard     nproc          13107" | sudo tee -a /etc/security/limits.conf
echo "root soft     nofile         13107" | sudo tee -a /etc/security/limits.conf
echo "root hard     nofile         13107" | sudo tee -a /etc/security/limits.conf
echo "session required pam_limits.so" | sudo tee -a /etc/pam.d/common-session
ulimit -n 13107
ulimit -Hn 13107
echo "ulimit:"
ulimit -n


# setup ebs
sudo mkfs -t ext4 /dev/xvdf
sudo mkdir /storage
sudo mount /dev/xvdf /storage/
sudo chmod 777 /storage/

# Download files
cd /storage/ || die "Could not go to storage dir"

# clone str-imputation repository
git clone https://github.com/shubhamsaini/str-imputation.git || die "Could not clone git repo"
cd str-imputation/L1O_Analysis/

# download shapeit haplotype
aws s3 cp s3://gatk-results/shapeit_haplotypes/shapeit.chr${CHROM}.with.ref.reorder.vcf.gz . || die "Could not download shapeit_haplotypes"
aws s3 cp s3://gatk-results/shapeit_haplotypes/shapeit.chr${CHROM}.with.ref.reorder.vcf.gz.csi .
bcftools index shapeit.chr${CHROM}.with.ref.reorder.vcf.gz

# download phased hipstr panel
aws s3 cp s3://gatk-results/hipstr_phased/hipstr.chr${CHROM}.phased.vcf.gz . || die "Could not download shapeit_haplotypes"
aws s3 cp s3://gatk-results/hipstr_phased/hipstr.chr${CHROM}.phased.vcf.gz.csi .
bcftools index hipstr.chr${CHROM}.phased.vcf.gz

# run L1O L1O_Analysis
./runme.sh hipstr.chr${CHROM}.phased.vcf.gz shapeit.chr${CHROM}.with.ref.reorder.vcf.gz snp.str.chr${CHROM}.vcf.gz ${CHROM}

aws s3 cp l1o.results.chr${CHROM}.csv ${OUTBUCKET}/l1o_log/
aws s3 cp l1o.imputed.chr${CHROM}.str.vcf.gz ${OUTBUCKET}/l1o_log/
aws s3 cp l1o.imputed.chr${CHROM}.str.vcf.gz ${OUTBUCKET}/l1o_log/

terminate
