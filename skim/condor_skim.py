#!/usr/bin/env python3
"""
This program takes care of skimming of DTNtuples:
- it checks for the presence of already skimmed files
- it allows parallel skimming using HTCondor
- it provides a status summary of the skim process
"""

import argparse
from datetime import datetime
from os import path, makedirs 
from sys import exit
import subprocess
import glob

#----------------
# Variables
#----------------

DEBUG = False

FILES_PER_BLOCK = 5

EOS_BASE_FOLDER = "/eos/user/c/cbattila/snd_analysis/"

INPUT_FOLDERS = { "TB" : "/eos/experiment/sndlhc/convertedData/commissioning/testbeam_June2023_H8/",
                  "TI18" : "/eos/experiment/sndlhc/convertedData/physics/2023/"}

#-----------------
# Helper functions
#-----------------

def non_skimmed_files(in_folder, out_folder):
    
    results = []
    
    paths_in = glob.glob(path.join(in_folder, "*.root"))
    files_in = [path.basename(path_in) for path_in in paths_in]

    paths_out = glob.glob(path.join(out_folder, "*.root"))
    files_out = [path.basename(path_out) for path_out in paths_out]

    for file_in, path_in in zip(files_in, paths_in):
        if file_in not in files_out:
            results.append(path_in)

    return results

def condor_create_sh(job_folder, files, out_folder, i_block):

    sh_file_name = path.join(job_folder, f"run_skim_{i_block}.sh")
    sh_file = open(sh_file_name,"w")

    sh_file.write("#! /usr/bin/bash\n")
    sh_file.write("cd /afs/cern.ch/user/c/cbattila/private/\n")
    sh_file.write("source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/setUp.sh\n")
    sh_file.write("eval `alienv load sndsw/latest`\n")
    sh_file.write("cd /afs/cern.ch/user/c/cbattila/private/sndlhc_bo_tbanalysis/skim/\n")

    for file in files:
        command = ["./run_skim.sh", 
                   f"{file}",
                   f"{out_folder}",
                   "true",
                   "\n"]

        sh_file.write(" ".join(command))

    sh_file.close()

    return sh_file_name

def condor_create_jdl(job_folder, sh_file_name, i_block):

    jdl_file_name = path.join(job_folder, f"run_skim_{i_block}.jdl")
    jdl_file = open(jdl_file_name, "w")

    jdl_file.write(f"executable = {sh_file_name}\n")
    jdl_file.write("universe = vanilla\n")
    jdl_file.write(f"output = {job_folder}/condor_{i_block}_$(Cluster)_$(Process).out\n")
    jdl_file.write(f"error  = {job_folder}/condor_{i_block}_$(Cluster)_$(Process).err\n")
    jdl_file.write(f"log    = {job_folder}/condor_{i_block}_$(Cluster)_$(Process).log\n")
    jdl_file.write("queue 1\n")

    jdl_file.close()

    return jdl_file_name

def condor_skim_command(jdl_file_name):

    command = ["condor_submit", jdl_file_name]

    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    output, error = process.communicate()

    if DEBUG:
        print("[condor OUTPUT]:")
        print(output.decode("utf8"))

    if error:
        print("[condor ERROR]:")
        print(error.decode("utf8"))

    return


def condor_skim(job_folder, run, type):

    in_folder = path.join(INPUT_FOLDERS[type],f"run_{run}")
    out_folder = path.join(EOS_BASE_FOLDER, f"{type}/run_{run}")

    files_tbp = non_skimmed_files(in_folder,out_folder)
    file_blocks = [files_tbp[i_file:i_file + FILES_PER_BLOCK] for i_file in range(0, len(files_tbp), FILES_PER_BLOCK)]
    
    for file_block, i_block in zip(file_blocks, range(1, len(file_blocks)+1)):
        
        if DEBUG:
            print(f"[condor_skim] creating .sh script for file block # {i_block}.")
        sh_file_name = condor_create_sh(job_folder, file_block, out_folder, i_block)
        
        if DEBUG:
            print(f"[condor_skim] creating condor .jdl file for .sh script # {i_block}.")
        jdl_file_name = condor_create_jdl(job_folder, sh_file_name, i_block)
        
        if DEBUG:
            print(f"[condor_skim] running condor_submit for .jdl file # {i_block}.")
        condor_skim_command(jdl_file_name)

        print(f"Submitted job [{i_block:3} / {len(file_blocks)}]")
    return

def status(run, type):
    
    in_folder = path.join(INPUT_FOLDERS[type],f"run_{run}")
    out_folder = path.join(EOS_BASE_FOLDER, f"{type}/run_{run}")

    if not path.isdir(out_folder):
        print(f"[status] folder : {out_folder} does not exist.")
        exit(100)

    n_files_to_process = len(non_skimmed_files(in_folder, out_folder))

    if n_files_to_process:    
        print(f"[status] {n_files_to_process} files are not yet skimmed.")
        print("[status] if you have no running jobs, please re-run this macro.")
    else:
        print("[status] all files have been skimmed, you can now proceed running hadd.")

    return

if __name__ == '__main__':

    #----------------------
    # Setup argument parser
    #----------------------

    PARSER = argparse.ArgumentParser(description=__doc__)
    
    PARSER.add_argument("command",
                        help="Either: 'condor_skim' or 'status'")
    
    PARSER.add_argument("run_number",
                        help="run number to be analysed")
    
    PARSER.add_argument("type",
                        help="run number to be analysed")

    ARGS = PARSER.parse_args()

    #----------------
    # Variables
    #----------------

    COMMANDS = ["condor_skim", "status"]

    if ARGS.command not in COMMANDS:
        print(f"[generic] command : {ARGS.command} must be either in: {COMMANDS}.")
        exit(100)

    if ARGS.type not in INPUT_FOLDERS.keys():
        print(f"[generic] command : {ARGS.type} must be either in: {INPUT_FOLDERS.keys()}.")
        exit(100)


    job_time = datetime.now().strftime("%d%m%Y_%H%M%S")
    job_folder = path.join("./jobs", f"run_{ARGS.run_number}_{ARGS.type}_{job_time}")
    skim_folder = path.join(EOS_BASE_FOLDER, f"{ARGS.type}/run_{ARGS.run_number}")

    OUT_FOLDERS = []
    
    if "condor_skim" in ARGS.command:
        OUT_FOLDERS.append(job_folder)
        OUT_FOLDERS.append(skim_folder)

    for folder in OUT_FOLDERS:
        folder_path = path.join(folder)
        if not path.exists(folder_path):
            makedirs(folder_path)

    if ARGS.command == "condor_skim":
        condor_skim(job_folder, ARGS.run_number, ARGS.type)
    elif ARGS.command == "status":
        status(ARGS.run_number, ARGS.type)
        
