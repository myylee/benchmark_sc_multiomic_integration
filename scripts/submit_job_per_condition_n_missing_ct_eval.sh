#!/bin/bash

#BSUB -J dweisepytest #LSF Job Name
#BSUB -q mingyao_normal
#BSUB -o pytestdweise.%J.txt #Name of the job output file
### -- Default: use 8 cores --
#BSUB -n 8
#BSUB -R "span[hosts=1]"
### -- specify that we need 32GB of memory per core/slot --
#BSUB -R "rusage[mem=32GB]"
### -- specify that we want the job to get killed if it exceeds 32 GB per core/slot --
#BSUB -M 32GB
### -- send notification at completion --
#BSUB -N



############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "general script to run a python/R script (with two arguments) under a certain conda env"
   echo
   echo "Syntax: scriptTemplate [-i|w|c|s|p|r|e|f|t|l]"
   echo "options:"
   echo "i     input address."
   echo "w     output address"
   echo "c     conda environment name"
   echo "s     script path"
   echo "p     script is written in Python"
   echo "r     script is written in R"
   echo "e     eval script path"
   echo "m     method key"
   echo "f     path to result matrix (specific to each method)"
   echo "t     path to barcode-to-cell type table"
   echo "l     number of clusters"
   echo "a     path to R script that performs gene-peak association using the predicted gene/peak expression"
   echo "b     path to gene-peak association pair list, or false if no association to be evaluated"
   echo "g     path to file storing the rare cell types being evaluated in this script (must be present, no default)"
   echo
}

############################################################
############################################################
# Main program                                             #
############################################################
############################################################

py="false"
r="false"

#  If a character is followed by a colon (e.g. f:), that option is expected to have an argument. thus here, p and r are not expected to have an argument.
while getopts i:w:c:s:pre:f:m:t:l:a:b:g: flag
do
    case "${flag}" in
        i) in_dir=${OPTARG};;
        w) out_dir=${OPTARG};;
        c) conda_env=${OPTARG};;
        s) script_path=${OPTARG};;
        p) py="true";;
        r) r="true";;
        e) eval_path=${OPTARG};;
        f) file_path=${OPTARG};;
        m) method_key=${OPTARG};;
        t) ct_ref=${OPTARG};;
        l) nclust=${OPTARG};;
        a) eval_path_gp=${OPTARG};;
        b) gp_truth=${OPTARG};;
        g) rare_ct_path=${OPTARG};;
    esac
done

echo "Working Directory: $PWD";
echo "Input Directory: $in_dir";
echo "Output Directory: $out_dir";
echo "Conda Environment: $conda_env";
echo "Script: $script_path";
echo "Running in Python: $py";
echo "Running in R: $r";
echo "Evaluation script: $eval_path";
echo "Result matrix path: $file_path";
echo "Method key: $method_key";
echo "Cell type reference path: $ct_ref";
echo "nclust: $nclust";
echo "gene-pair eval script: $eval_path_gp";
echo "gene-pair truth list: $gp_truth";
echo "path to rare cell type list: $rare_ct_path";

#module load R 
source ~/anaconda3/etc/profile.d/conda.sh
conda activate $conda_env

echo "Conda env activated";

#create output directory if it doesn't exist 
mkdir -p $out_dir


# white spaces in between everything! 
if [[ $py == "true" ]] && [[ $r == "false" ]]
then
    echo "running a Python script."
    python $script_path $in_dir $out_dir
elif [[ $py == "false" ]] && [[ $r == "true" ]]
then
    echo "running a R script."
    Rscript --vanilla $script_path $in_dir $out_dir $nclust
else 
    echo "please specify the language for the script to be run in, can either be python or r but not both or none."
fi 

echo "Running evaluation";
conda deactivate
conda activate scib2

if [[ "$file_path" == *"seurat4"* ]] || [[ "$file_path" == *"liger"* ]]
then
    echo "Not clustering" 
    python $eval_path $out_dir $file_path $ct_ref $rare_ct_path
else 
    echo "Will cluster using latent embedding" 
    python $eval_path $out_dir $file_path $ct_ref $rare_ct_path $nclust
fi 

if [[ "$eval_path_gp" != "false" ]] || [[ "$gp_truth" != "false" ]]
then
    echo "Evaluating predicted profiles"
    Rscript --vanilla $eval_path_gp $in_dir $out_dir $method_key $gp_truth
else 
    echo "Not evaluating predicted profiles, END!"
fi  
    

