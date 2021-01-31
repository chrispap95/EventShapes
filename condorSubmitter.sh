for i in HT1000to1500 HT1500to2000 HT2000toInf;
do
cat > condor_${i}.jdl << "EOF"
universe = vanilla
Executable = condor-exec.csh
EOF
echo "Arguments = ${i} bkg CMSSW_11_2_0" >> condor_${i}.jdl
cat >> condor_${i}.jdl << "EOF"
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
request_cpus = 4
Transfer_Input_Files = condor-exec.csh, CMSSW_11_2_0.tgz
Output = logs/python_$(Cluster)_$(Process).stdout
Error = logs/python_$(Cluster)_$(Process).stderr
Log = logs/python_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
EOF
condor_submit condor_${i}.jdl
done

for i in mMed-1000_mDark-2_temp-2_decay-darkPhoHad;
do
cat > condor_${i}.jdl << "EOF"
universe = vanilla
Executable = condor-exec.csh
EOF
echo "Arguments = ${i} sig CMSSW_11_2_0" >> condor_${i}.jdl
cat >> condor_${i}.jdl << "EOF"
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
request_cpus = 4
Transfer_Input_Files = condor-exec.csh, CMSSW_11_2_0.tgz
Output = logs/python_$(Cluster)_$(Process).stdout
Error = logs/python_$(Cluster)_$(Process).stderr
Log = logs/python_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
EOF
condor_submit condor_${i}.jdl
done
