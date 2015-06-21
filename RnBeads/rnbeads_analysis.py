#!/usr/bin/env python
__author__ = 'Bjorn Wouters'
__email__ = "bjorn-wouters@hotmail.com"

"""
Description: Automatic adding and/or analysis of a specific genome in the R RnBeads package.
Version: 1.0.0.
Known bugs:
- stdout fix for Galaxy (not showing in Galaxy except on the end).
Dependencies:
- RnBeads package (edited for this script).
- Template package.
- prepare_analysis script.
- Assembly template folder.
Work to be done:
- Adding of the assembly to the source package isn't stable (can't delete assemblies folder or you'll get an error).
- Adding of more options for the analysis instead of using the default.
- Automatic sample file creation.
"""

import subprocess
import sys
import os
from shutil import rmtree, move, copyfile
import tarfile
import argparse
import prepare_analysis
import check_file_formats
import time


def main():
    tmp_files = list()
    # Parse all the arguments given be the argparse module.
    args = parse_args()
    try:
        # Means the argparse subprocess "only_analysis" is chosen.
        if "fasta" not in args:
            # Run analysis with (optional) added annotation.
            analysis_script = run_analysis(args)
            tmp_files.append(analysis_script)
            sys.exit(0)
        # Checks the files on the right syntax.
        check_file_formats.check_fasta(args.fasta)
        check_file_formats.check_bed(args.bed)
        check_file_formats.check_sample_file(args.sample_file)
        # Prepare analysis files
        if args.bed:
            prepare_bed_analysis(args)
        # Coverts given .Fasta file to a 2bit file and writes it to the /tmp/ directory.
        twobit_file = fasta_to_2bit(args)
        tmp_files.append(twobit_file)
        # Makes description file out of the species and genus information given by the argparse arguments and
        # writes it to the /tmp/ directory.
        description = make_rnbeads_description(twobit_file, args)
        tmp_files.append(description.name)
        # Creates the new package of the given fasta via the BSgenome.forge method in R.
        #remove all existing BSgenome.* in /tmp
        cmd = "rm -rf %s"%os.path.join(args.temp_directory,'BSgenome*')
        run_subprocess(cmd,"removing all existing BSgenome packages")
        forge_script = forge_genome_file(description, args)
        tmp_files.append(forge_script)
        # Installs the genome file to R on your system.
        genome_folder, genome_package = install_genome_file(args)
        tmp_files.extend([genome_folder, genome_package])
        # Append the RnBeads sourcecode with the new assembly.
        append_source_code(args, genome_folder)
        # Append assembly file of the RnBeads package.
        append_script = append_assembly(args)
        tmp_files.append(append_script)
        # Makes the description file for the given assembly.
        make_assembly_description(args, genome_folder)
        # Install the new RnBeads package so appended source code will be used in the analysis.
        install_rnbeads(args)
        # Calculates all the CpG sites for each chromosomes and appends it to the assembly package.
        site_data, r_script, region_file = get_cpg_sites(args, genome_folder)
        tmp_files.extend([site_data, r_script, region_file])
        # Installs the given assembly with the annotated site data
        install_assembly(args)
        # Means the "add_and_analysis" is chosen.
        if "cores" in args:
            # Run analysis with (optional) added annotation.
            analysis_script = run_analysis(args)
            tmp_files.append(analysis_script)
    finally:
        # Clear all files written to the tmp directory.
        clear_tmp(tmp_files)


def run_subprocess(cmd, log_message):
    """
    Run subprocess under standardized settings
    force the cmds to be a string
    """
    sys.stdout.write("now starting:\t%s\n\n" % log_message)
    sys.stdout.write('running:\t%s\n\n' % cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
    exit_code = p.wait()
    stdout = p.stdout.read().replace('\r', '\n')
    stderr = p.stderr.read().replace('\r', '\n')
    if stdout:
        sys.stdout.write('stdout:\n%s\n' % stdout)
    if stderr:
        sys.stdout.write('stderr:\n%s\n' % stderr)
    sys.stdout.write('finished:\t%s\n\n' % log_message)
    if exit_code:
        return Exception("Call of %s failed with \n %s" % (cmd, stderr))
    else:
        return 0


def compress_folder(args, cur_time):
    """
    For Galaxy: After finishing the analysis, the folder will be zipped in tar.gz format.
    Replaces the Galaxy .dat file with the compressed analysis folder.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Gets the folder destination of the current script.
    # script_dir = script_dir.replace(' ','\ ')
    analysis_folder = os.path.join(script_dir, "assemblies", args.assembly_code, cur_time)
    if not os.path.isdir(analysis_folder):
        return
    tar_file = tarfile.open(os.path.join(args.temp_directory, cur_time+"_analysis.tar.gz"), "w:gz")
    tar_file.add(analysis_folder, arcname=cur_time+"_analysis")
    tar_file.close()
    move(tar_file.name, args.output)


def prepare_bed_analysis(args):
    """
    :argument: args, all arguments from argparse.
    Preparing analysis which means:
    - Making .bed files with each EPP annotated for every sample.
    - Making chromosomes file with each valid chromosome.
    - Making the analysis folder.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Gets the folder destination of the current script.
    # script_dir = script_dir.replace(' ','\ ')
    sys.stdout.write("Starting: Preparing analysis files.\n")
    output_dir = os.path.join(script_dir, "assemblies", args.assembly_code, "bs.bed/")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Create folder if not already exists.
    # Copy's the given sample file to the analysis folder.
    copyfile(args.sample_file, os.path.join(script_dir, "assemblies", args.assembly_code, "sample.csv"))

    sample_file = open(os.path.join(script_dir, "assemblies", args.assembly_code, "sample.csv"))
    given_samples = list(sample.split(",")[0] for sample in sample_file)[1:]
    sample_file.close()

    input_file, output_dict, samples = prepare_analysis.ParseFiles(args.bed, output_dir)  # Make all .bed files
    # Fill all .bed files with formatted info.
    invalid_samples = prepare_analysis.IgvToRnBeads(input_file, output_dict, samples, output_dir, given_samples,
                                                    args.minimal_reads)

    # If there are invalid samples (that have less than 5% of the reads of the maximum sample`s reads), Then
    # there will be a new samples file created with these filtered out.
    if invalid_samples:
        original_sample_file = open(sample_file.name)
        lines = original_sample_file.readlines()
        original_sample_file.close()
        new_sample_file = open(sample_file.name, "w")
        for line in lines:
            if line.split(",")[0] not in invalid_samples:
                new_sample_file.write(line)
        new_sample_file.close()

        sys.stdout.write("There are: " + str(len(invalid_samples)) + """ samples filtered out of the sample file because
        they have to few sites.\nSamples that are filtered out: """ + " and ".join(invalid_samples) + "\n")

    sys.stdout.write("Finished: Preparing analysis files.\n")
    # sys.exit(0)


def run_analysis(args):
    """
    :argument: args, all arguments from argparse.
    :return: Analyse directory.
    Analysis function for the R RnBeads package. Prepares the template .R script run.analysis.
    After the analysis, the output folder will be zipped in a given directory.
    """
    cur_time = time.strftime("%d_%m_%Y_%H:%M")
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Gets the folder destination of the current script.
    # script_dir = script_dir.replace(' ','\ ')
    package_name = "".join(["BSgenome.", args.species_name, args.genus_name, ".NIOO.v", args.version])
    if args.annotation:  # Creates R readable code for the file locations of the analysis.
        annotation_files = list()
        annotation_names = list()
        #TODO: adding multiple annotation files in galaxy does not yet work.
        for annotation_file in args.annotation:
            annotation_files.append('"'+annotation_file+'"')
            annotation_names.append('"'+os.path.basename(os.path.splitext(annotation_file)[0])+'"')
        annotation_files = ",".join(annotation_files)
        annotation_names = ",".join(annotation_names)
    else:
        annotation_files = "NULL"
        annotation_names = "NULL"
    # Load template variable files and fill with correct parameters from script_dict
    script_dict = {
        "annotation_files": annotation_files,
        "annotation_names": annotation_names,
        "assembly": args.assembly_code,
        "package": os.path.basename(package_name),
        "species": args.species_name,
        "lib_path": args.lib_path,
        "directory": os.path.join(script_dir, "assemblies", args.assembly_code+"/"),
        "cores": args.cores,
        "time": cur_time}
    # Puts the dict in the file.
    template_script = open(os.path.join(script_dir,"templates","run.analysis.R")).read() % script_dict
    # Write filled in analysis R script to tmp dir.
    r_script = open(os.path.join(args.temp_directory, "run.analysis.R"), "w")
    r_script.write(template_script)
    r_script.close()
    # Load r script in R and run consecutively without user input
    command = "R < "+r_script.name+" --no-save"
    log_message = "Running RnBeads analysis, this could take a while."
    run_subprocess(command, log_message)  # Running the main analysis.

    if args.output:
        # Compresses the RnBeads output folder and changes it with the .dat file Galaxy
        compress_folder(args, cur_time)

    return r_script.name


def install_assembly(args):
    """
    Install the template assembly package with the given data and DESCRIPTION file.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Gets the folder destination of the current script.
    script_dir = script_dir.replace(' ','\ ')
    assembly_package = script_dir+"/assembly"
    # If the library path consists out of a empty list (R syntax) than the packages will be written to the standard
    # package installation folder (mostly /usr/local/bin/R). Otherwise it will be written to the specific folder.
    if args.lib_path == "c()":
        specific_path = ""
    else:
        specific_path = " -l " + args.lib_path + " ; R_LIBS=" + args.lib_path + "; export R_LIBS"
    command = "R CMD INSTALL " + assembly_package + specific_path
    log_message = "Installing the new assembly package"
    run_subprocess(command, log_message)  # Installs the assembly with code:command


def get_cpg_sites(args, package_name):
    """
    :argument: args, all arguments from argparse.
    :param package_name: Name of the package.
    Creates the sites and the region output for the assembly.
    The sites.RData will contain the information of each CpG position and the regions.RData will contain
    data of the optional annotation data.
    After the installation of the assembly, the files will be deleted.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Gets the folder destination of the current script.
    # script_dir = script_dir.replace(' ','\ ')
    output_file_name = args.assembly_code+".CpG.RData"
    output_dir = os.path.join(script_dir,'assembly','data',output_file_name)
    #check if data exists as a subfolder of assembly!
    if not os.path.exists(os.path.join(script_dir,'assembly','data')):
        os.mkdir(os.path.join(script_dir,'assembly','data'))
    region_output = os.path.join(script_dir.replace(' ','\ '),'assembly','data',args.assembly_code+'.regions.RData')
    script_dict = {
        "assembly": args.assembly_code,
        "package": os.path.basename(package_name),
        "species": args.species_name,
        "sites_output": output_dir.replace(' ','\ '),
        "regions_output": region_output,
        "lib_path": args.lib_path,
        "chromosomes": os.path.join(script_dir, "assemblies", args.assembly_code+"/bs.bed/chromosomes.txt")}
    template_script_file = open(os.path.join(script_dir, "templates/get.sites.R"))
    # Puts the dict in the string.
    template_script = template_script_file.read() % script_dict
    template_script_file.close()
    # Makes a temporary R script and executes it after it is appended with the template script.
    r_script = open(os.path.join(args.temp_directory,"get.sites.R"), "w")
    r_script.write(template_script)
    r_script.close()
    command = "R < "+r_script.name+" --no-save"
    log_message = "Adding CpG sites to assembly annotation."
    run_subprocess(command, log_message)  # Runs the annotation script.
    return [output_dir, r_script.name, region_output]  # For the cleanup of the /tmp/ folder.


def make_assembly_description(args, package_name):
    """
    :argument: args, all arguments from argparse.
    :param package_name: Name of the package.
    Makes the description for the new assembly. The description is already made in a template though
    the assembly and the package name needs to be specified.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Gets the folder destination of the current script.
    # script_dir = script_dir.replace(' ','\ ')
    description_dict = {
        "assembly": args.assembly_code,
        "package": os.path.basename(package_name)}
    # Opens the template and puts the missing values in the string and writes it to a temporary file.
    template_description_file = open(os.path.join(script_dir, "templates/assembly_description.DCF"))
    template_description = template_description_file.read() % description_dict
    template_description_file.close()
    new_description = open(args.temp_directory+"DESCRIPTION", "w")
    new_description.write(template_description)
    new_description.close()
    move(new_description.name, script_dir+"/assembly/DESCRIPTION")  # Puts the new description file in the folder.


def install_rnbeads(args):
    """
    Installs the appended RnBeads package to the /home/R folder of the user.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    script_dir = script_dir.replace(' ','\ ')
    rnbeads_package = script_dir+"/RnBeads"
    # If the library path consists out of a empty list (R syntax) than the packages will be written to the standard
    # package installation folder (mostly /usr/local/bin/R). Otherwise it will be written to the specific folder.
    if args.lib_path == "c()":
        specific_path = ""
    else:
        specific_path = " -l " + args.lib_path + " ; R_LIBS=" + args.lib_path + "; export R_LIBS"
    command = "R CMD INSTALL " + rnbeads_package + specific_path
    log_message = "reinstalling the appended RnBeads package"
    run_subprocess(command, log_message)  # Runs the process.


def append_assembly(args):
    """
    :argument: args, all arguments from argparse.
    Appends the RnBeads package with the assembly annotation and appends the sourcecode to its destination folder.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Gets the folder destination of the current script.
    # script_dir = script_dir.replace(' ','\ ')
    annotation_file = "".join(['"', script_dir, "/RnBeads/data/annotations.RData", '"'])
    assembly_dict = {
        "file": annotation_file,
        "assembly": args.assembly_code,
        "tmp": args.temp_directory
    }
    # Opens the template and puts the missing values in the string and writes it to a temporary file.
    template_script_file = open(os.path.join(script_dir, "templates/add.assembly.R"))
    template_script = template_script_file.read() % assembly_dict
    template_script_file.close()
     # Writes the new string to a temporary .R file and at the end of the script it will be deleted.
    add_assembly_script = open(os.path.join(args.temp_directory, "add.assembly.R"), "w")
    add_assembly_script.write(template_script)
    add_assembly_script.close()
    command = "R < "+add_assembly_script.name+" --no-save"
    log_message = "Adding assembly to RnBeads package."
    run_subprocess(command, log_message)  # Executes the command and appends the assembly data.
    move(os.path.join(args.temp_directory,"annotations.RData"),\
    os.path.join(script_dir,"RnBeads","data","annotations.RData"))
    return add_assembly_script.name


def append_source_code(args, folder_name):
    """
    Appends new genome Rfile to the existing R source code so that Rnbeads can be run on the 'new' assembly.
    """
    #TODO: investigate safer / better method to check for doubles / existing genome assembly names. Do not add
    #TODO: assemblies that already exist!
    chrom_sizes = prepare_analysis.chrom_sizes(args.fasta, args.temp_directory, args.assembly_code)
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Gets the folder destination of the current script.
    # script_dir = script_dir.replace(' ','\ ')
    move(chrom_sizes.name, script_dir+"/RnBeads/inst/extdata/chromSizes/"+os.path.basename(chrom_sizes.name))
    chromosome_var = "".join(["\n", args.assembly_code, '.chr <- read.table("',
                              os.path.join(script_dir, "assemblies", args.assembly_code+"/"),
                              'bs.bed/chromosomes.txt")\n', '##%(chromosomes)s'])
    assembly_var = "".join(["\n,", "'"+args.assembly_code+"'", " = ", args.assembly_code, ".chr[[1]]\n",
                           '##%(assembly_table)s'])
    package_var = "".join(["\nelse if (assembly == ", "'", args.assembly_code, "') {\n",
                           "suppressPackageStartupMessages(require(", os.path.basename(folder_name), "))\n",
                           "genome.data <- ", args.species_name, "\n}\n", "##%(assembly_package)s"])
    annotation_dict = {
        "chromosomes": chromosome_var,
        "assembly_table": assembly_var,
        "assembly_package": package_var}
    # Reads the template script and puts the variable of the annotation dict in the script as a string.
    annotation = open(script_dir+"/RnBeads/R/assemblies.R").read() % annotation_dict
    # Makes a temporary file where the new changed R script will be written to and executed. Will be deleted at
    # the end of the script.
    new_annotation = open(script_dir+"/RnBeads/R/assemblies.R", "w")
    new_annotation.write(annotation)
    new_annotation.close()


def install_genome_file(args):
    """
    Installs and builds the new BSgenome (R library) for the given .fasta file.
    """
    tmp = args.temp_directory

    # Opens the DESCRIPTION file for getting the package name of the genome.
    description = open(os.path.join(args.temp_directory,"description.DCF"))
    genome_folder = description.readline().strip("Package: \n")
    description.close()

    version_name = args.version+".0"  # Syntax needs to be version + .0
    command = ''.join(["cd ", tmp, ";R CMD build ", genome_folder, "/"])
    log_message = "Building genome file"
    run_subprocess(command, log_message)  # Builds the genomic annotation file.
    genome_package = ''.join([genome_folder, "_", version_name, ".tar.gz"])
    # If the library path consists out of a empty list (R syntax) than the packages will be written to the standard
    # package installation folder (mostly /usr/local/bin/R). Otherwise it will be written to the specific folder.
    #c() is an empty list in R.
    if args.lib_path == "c()":
        specific_path = ""
    else:
        specific_path = " -l " + args.lib_path + " ; R_LIBS=" + args.lib_path + "; export R_LIBS"
    command = "".join(["cd ", tmp, ";R CMD INSTALL ", genome_package, specific_path])
    log_message = "Installing genome file"
    run_subprocess(command, log_message)  # After building the folder to a tar, the package can be installed.
    return [genome_folder, genome_package]


def forge_genome_file(description, args):
    """
    Creates the new package of the given fasta via the BSgenome.forge method in R.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Gets the folder destination of the current script.
    #script_dir = script_dir.replace('Thomas ','Thomas\ ')
    file_dict = {"DCF": description.name,
                 "tmp": args.temp_directory}
    file_template = open(os.path.join(script_dir, "templates/forge.genome.R"))
    forge_template = file_template.read()
    file_template.close()
    adjusted_forge_template = forge_template % file_dict  # Puts the dict in the template string.
    # Makes an empty R script and fills it with the appended template script.
    # After it is appended, the script will be executed and afterwards it will be deleted.
    forge_r_script = open(os.path.join(args.temp_directory,"forge.genome.R"), "w")
    forge_r_script.write(adjusted_forge_template)
    forge_r_script.close()
    # Makes from the description file, a folder which can be build and installed via bash.
    command = "R < "+forge_r_script.name+" --no-save"
    log_message = "Forging genome file"
    run_subprocess(command, log_message)  # Runs the appended template script.
    return forge_r_script.name


def make_rnbeads_description(twobit_file, args):
    """
    Makes the description file for the RnBeads package.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Gets the folder destination of the current script.
    #script_dir = script_dir.replace(' ','\ ')
    file_dict = {
        "species": args.species_name,
        "genus": args.genus_name,
        "twobit_dir": args.temp_directory,
        "twobit_name": os.path.basename(twobit_file),
        "templ_dir": args.temp_directory,
        "version": args.version}
    with open(os.path.join(script_dir, "templates/DESCRIPTION_TEMPLATE.DCF")) as description_template:
        description = description_template.read()
        adjusted_template = description % file_dict  # Puts the dict in the string gotten from the template.
        # Writes the adjusted template to a new file that is written to the given temp folder.
        description = open(os.path.join(args.temp_directory,"description.DCF"), "w")
        description.write(adjusted_template)
        description.close()
    return description


def fasta_to_2bit(args):
    """
    Coverts the given fasta to a .2bit file via the faToTwoBit executable.
    """
    #TODO: build fato2bit from source on installation to make sure it runs.
    #Source for fat2bit for mac osx is here: http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/
    #TODO: or list dependency in help and check for executable upon running analysis.
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Gets the folder destination of the current script.
    script_dir = script_dir.replace(' ','\ ')
    sys.stdout.write("""Adding the genome of """+args.species_name+"""" to the RnBeads package\n
                Starting with converting the .fasta to a .2bit file.\n""")
    fasta = args.fasta
    twobit_name = os.path.basename(os.path.splitext(fasta)[0])+'.2bit'
    twobit_file = os.path.join(args.temp_directory, twobit_name)
    # Writes the fasta to a twobit file in the tmp folder and is deleted after the analysis.
    log_message = "Converts the given fasta to a .2bit file"
    command = " ".join([os.path.join(script_dir, "templates/faToTwoBit"), fasta, twobit_file])
    run_subprocess(command, log_message)
    return twobit_file


def parse_args():
    """
    Processes the given argparse arguments.
    """
    parser = argparse.ArgumentParser(description="""Add new genome to the RnBeads R module.
     Or, if already done, run only the analysis.""")
    subparsers = parser.add_subparsers(help='Choose your analysis type:')
    #TODO: add option to list already run analysis for existing genomes. Skip existing genomes.
    # If chosen: only the run_analysis function will be executed.
    analysis_only = subparsers.add_parser('analysis_only', help="""If assembly and genome already added:
    analysis only is possible.""")
    analysis_only.add_argument('-c', '--cores', help='Number of cores you want to use for the analysis', default="2")
    analysis_only.add_argument('-a', '--annotation', default=None, nargs="*",
                               help="""Annotation file(s) to be added for the annotation (optional).
                               Extra files can be added by separating them with a whitespace.
                               File syntax: e.g. chr1\nchr2\n.
                               The name of the annotation will be the same as the filename.""")
    analysis_only.add_argument('-tmp', '--temp_directory', help='temp directory', default="/tmp/")
    analysis_only.add_argument('-s', '--species_name', help='Name of the species belonging to the genome')
    analysis_only.add_argument('-g', '--genus_name', help='Name of the genus from where the organism stems from')
    analysis_only.add_argument('-v', '--version', help='Version of the genome to be forged', default="1")
    analysis_only.add_argument('-ac', '--assembly_code', help="""Assembly code used for your organism in the RnBeads
                                                              package""")
    analysis_only.add_argument('-lp', '--lib_path', help='Library installation folder for R packages.', default="c()")
    analysis_only.add_argument('-o', '--output', help='Output file (needed for galaxy)', default=None)

    # If chosen all the main functions will be executed.
    add_and_analysis = subparsers.add_parser('add_and_analysis', help="""Add genome and assembly AND analyse it
                                                                afterwards.""")
    add_and_analysis.add_argument('-c', '--cores', help='Number of cores you want to use for the analysis', default="2")
    add_and_analysis.add_argument('-f', '--fasta', help='Fasta input file of the new genome', default=None)
    add_and_analysis.add_argument('-b', '--bed', help='Bed input file to be analysed. (optional if already done)',
                                  default=None)
    add_and_analysis.add_argument('-a', '--annotation', default=None, nargs="*",
                                  help="""Annotation file(s) to be added for the annotation (optional).
                                  Extra files can be added by separating them with a whitespace.
                                  File syntax: e.g. chr1\nchr2\n.
                                  The name of the annotation will be the same as the filename.""")
    add_and_analysis.add_argument('-sf', '--sample_file', help='Sample file location.')
    add_and_analysis.add_argument('-tmp', '--temp_directory', help='temp directory', default="/tmp/")
    add_and_analysis.add_argument('-s', '--species_name', help='Name of the species belonging to the genome')
    add_and_analysis.add_argument('-g', '--genus_name', help='Name of the genus from where the organism stems from')
    add_and_analysis.add_argument('-v', '--version', help='Version of the genome to be forged', default="1")
    add_and_analysis.add_argument('-ac', '--assembly_code', help="""Assembly code used for your organism in the RnBeads
                                                            package""")
    add_and_analysis.add_argument('-o', '--output', help='Output file (needed for galaxy)', default=None)
    add_and_analysis.add_argument('-lp', '--lib_path', help='Library installation folder for R packages.',
                                  default="c()")
    add_and_analysis.add_argument('-mr', '--minimal_reads', help='Number of minimal reads per sample on one CpG site',
                                  default=5)

    # If chosen: every main function fill be executed except for the run_analysis function.
    add_only = subparsers.add_parser('add_only', help="Only add the genome and assembly to the RnBeads package.")
    add_only.add_argument('-f', '--fasta', help='Fasta input file of the new genome', default=None)
    add_only.add_argument('-b', '--bed', help='Bed input file to be analysed. (optional if already done)', default=None)
    add_only.add_argument('-sf', '--sample_file', help='Sample file location.')
    add_only.add_argument('-tmp', '--temp_directory', help='temp directory', default="/tmp/")
    add_only.add_argument('-s', '--species_name', help='Name of the species belonging to the genome')
    add_only.add_argument('-g', '--genus_name', help='Name of the genus from where the organism stems from')
    add_only.add_argument('-v', '--version', help='Version of the genome to be forged', default="1")
    add_only.add_argument('-ac', '--assembly_code', help='Assembly code used for your organism in the RnBeads package')
    add_only.add_argument('-lp', '--lib_path', help='Library installation folder for R packages.', default="c()")
    add_only.add_argument('-mr', '--minimal_reads', help='Number of minimal reads per sample on one CpG site',
                          default=5)

    # Parses the arguments to the args variable.
    args = parser.parse_args()
    return args


def clear_tmp(trash):
    """
    Clears all the items that are appended to the trash variable.
    """
    for item in trash:
        if os.path.isdir(item):
            rmtree(item)
        else:
            try:
                os.remove(item)
            except OSError:
                pass


# Calls the main function if the script is executed.
if __name__ == '__main__':
    main()