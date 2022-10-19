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
shift
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


mycmd="printf '"
if [[ -v nthreads ]]; then
      mycmd+=" nthreads=${nthreads}\n"
fi
if [[ -v nwhitelist ]]; then
      mycmd+=" whitelist=${whitelist}\n"
fi
if [[ -v index ]]; then
      mycmd+=" index=${index}\n"
fi
if [[ -v index_dir ]]; then
      mycmd+=" index_dir=${index_dir}\n"
fi
if [[ -v tech ]]; then
      mycmd+=" tech=${tech}\n"
fi
if [[ -v cellthreshold ]]; then
      mycmd+=" cellthreshold=${cellthreshold}\n"
fi
if [[ -v project ]]; then
      mycmd+=" project=${project}\n"
fi
if [[ -v samplesheet ]]; then
      mycmd+=" samplesheet=${samplesheet}\n"
fi
if [[ -v mapper ]]; then
      mycmd+=" mapper=${mapper}\n"
fi
if [[ -v xmatrix ]]; then
      mycmd+=" xmatrix=${xmatrix}\n"
fi
if [[ -v postalign_dir ]]; then
      mycmd+=" postalign_dir=${postalign_dir}\n"
fi
if [[ -v createxmatrix ]]; then
      mycmd+=" createxmatrix=${createxmatrix}\n"
fi

mycmd+="' > Scasa_temp"
eval $mycmd



#printf "nthreads=${nthreads}\nwhitelist=${whitelist}\nindex=${index}\nindex_dir=${index_dir}" > Scasa_temp

path=`pwd`
params_path="${path}/Scasa_temp"


cmd="docker run -it"
if [ "$index" = "NO" ]; then
#    echo "USE INDEX DIR"
    cmd+=" -v ${index_dir}:/source/index_dir"
else
#    echo "USE INDEX FASTA"
    cmd+=" -v ${ref}:/source/referenceDB/refMrna.fa.gz"    
fi
cmd+=" -v ${INPUT}:/source/input"
cmd+=" -v ${OUTPUT}:/source/output"
cmd+=" -v ${params_path}:/source/params.sh"  
cmd+=" nghiavtr/scasa:v1.0.1"


eval $cmd
rm Scasa_temp
