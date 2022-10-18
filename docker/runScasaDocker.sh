#!/bin/sh

syntax="bash runScasaDocker.sh -param [parameter input file]"

echo -e "
                                                                                                                                                                           
---------------------------------------------------------------------------
                     You are running Scasa v1.0.1 using docker ....
---------------------------------------------------------------------------
"


while [[ $# -gt 1 ]]
do
key="$1"

case $key in
     
     -param|--parameter)
     param=$2
     shift
     ;;

     *)

esac

done

if [[ -n "$param" ]]; then
   echo ""
   echo "Loading parameters from file..."
   if [[ -f "$param" ]]; then
       source $param
   else
       echo "Your parameter file does not exist. Please check and try again!"
       exit
   fi
   
fi


printf "nthreads=${nthreads}\nwhitelist=${whitelist}\nindex=${index}\nindex_dir=${index_dir}" > Scasa_temp

path=`pwd`
params_path="${path}/Scasa_temp"


if [ -v index_dir ]; then
    cmd="docker run -it"
    cmd+=" -v ${ref}:/source/referenceDB/refMrna.fa.gz"
    cmd+=" -v ${INPUT}:/source/input"
    cmd+=" -v ${OUTPUT}:/source/output"
    cmd+=" -v ${index_dir}:/source/index_dir"
    cmd+=" -v ${params_path}:/source/params.sh"  
    cmd+=" nghiavtr/scasa:v1.0.1"
else
    cmd="docker run -it"
    cmd+=" -v ${ref}:/source/referenceDB/refMrna.fa.gz"
    cmd+=" -v ${INPUT}:/source/input"
    cmd+=" -v ${OUTPUT}:/source/output"
    cmd+=" -v ${params_path}:/source/params.sh"  
    cmd+=" nghiavtr/scasa:v1.0.1"
fi


eval $cmd
rm Scasa_temp
