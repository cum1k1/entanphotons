 #!/usr/bin/bash

#should launch as bash nprocs_calibration.sh ONLY!!!!!!!!!

export MACRO_PATH=/entanglement_new_files/macro
export SOURCE_PATH=/entanglement_new_files/entangled_scattered_low_gagg_zs/
export temp_DIR=calibrated_files/
rm -rf  $SOURCE_PATH/$temp_DIR
echo $SOURCE_PATH
mkdir -p  $SOURCE_PATH/$temp_DIR
counter=0
NUM_PROCS0=$(ps -A | grep root |wc -l)
echo NUM_PROCS0 $NUM_PROCS0
NUM_PROCS=4+$NUM_PROCS0/2
for filename in ${SOURCE_PATH}07a8de9a*.root;
do
    fname=$(basename $filename .root)
    echo $fname
    if [ -d \"$SOURCE_PATH/$temp_DIR/$fname\" ]
    then 
        echo "already calibrated"
    else
    ((counter=counter+1))
    mkdir -p  $SOURCE_PATH/$temp_DIR/$fname
    cp -r $MACRO_PATH  $SOURCE_PATH/$temp_DIR/$fname/macro
    cd $SOURCE_PATH/$temp_DIR/$fname/macro
    touch log.out
    nohup root -l -b -q "Waveforms_calibration_nprocs.cpp(\"${SOURCE_PATH}\",\"$filename\",\"$temp_DIR\",\"_$fname\")" &
    
    nprocs=$(ps -A | grep root |wc -l)
    echo $nprocs        
        while [[ $(ps -A | grep root |wc -l) -ge 2*$NUM_PROCS ]]
        do
            sleep 5
        done
    fi
done
wait
#hadd -f ${SOURCE_PATH}calibrated_time.root ${SOURCE_PATH}${temp_DIR}calibrated_time_*.root
cd $MACRO_PATH

#rm -rf files/
#touch data_info.txt
#echo "$RUNNAME">>data_info.txt

#rm `find . -name "Converted_fitted*_FIT_QA.pdf" ! -name "Converted_fitted_FIT_QA.pdf"`
#rm Converted_fitted*.root
#rm -rf  $SOURCE_PATH/$temp_DIR


#    nprocs=$(ps -A | grep root |wc -l)
