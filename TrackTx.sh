#!/bin/bash
cd -- "$(dirname "$0")" > /dev/null

#----------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------

### Instructions

    # This is the master script for aligning PROseq-data, mapping functional genomic regions, comparing gene expression, and studying sat III repeats.

    #----------------------------------------------------------------------------------------------------------------------------------------------

    # The input files for this script are:

    # All files can be downloaded from to generate bigWig files, you need to know which organism you're working with and which 
    # samples you're analysing (either their SRR numbers from the Gene Expression Omnibus or have the files localy)

    #----------------------------------------------------------------------------------------------------------------------------------------------

    # Required packages and how to install them are listed here:

    # Homebrew and conda are easiest to use if you're running MacOS and are not experienced with installing and adding packages to your PATH.

    #--------------------------------

    # These are required packages:

    # You need to install conda on your own. You can find instructions for conda here: https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html
    # then install the following packages if they are not installed already (if you don't the script will do it for you)

    # wget              |   conda install -c anaconda wget
    # bowtie2           |   conda install -c bioconda bowtie2
    # SRA-tools         |   conda install -c bioconda SRA-tools  
    # fastx_toolkit     |   conda install -c bioconda fastx_toolkit 
    # samtools          |   conda install -c bioconda samtools
    # bedtools          |   conda install -c bioconda bedtools
    # fetchchromsizes   |   conda install -c bioconda ucsc-fetchchromsizes
    # bedGraphToBigWig  |   conda install -c bioconda ucsc-bedGraphToBigWig

#----------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------

clear; sleep 1
echo "Welcome to TrackTx - the transcription tracking data analysis pipeline!"
sleep 1

### Check that the required packages are installed and working before starting the script

echo; echo "Checking if all required packages are installed..."
sleep 1

# List of commands to check
commands=("bedGraphToBigWig" "bowtie2" "fasterq-dump" "fastx_clipper" "samtools" "bedtools" "fetchChromSizes")

# Initialize an array to store non-installed commands
missing_commands=()

# Loop through the list of commands
for cmd in "${commands[@]}"; do
    if ! command -v "$cmd" &> /dev/null; then
        if [ "$cmd" == "bedGraphToBigWig" ]; then
            actual_cmd="ucsc-bedGraphToBigWig"
        elif [ "$cmd" == "fasterq-dump" ]; then
            actual_cmd="sra-tools"
        elif [ "$cmd" == "fastx_clipper" ]; then
            actual_cmd="fastx_toolkit"
        elif [ "$cmd" == "fetchChromSizes" ]; then
            actual_cmd="ucsc-fetchchromsizes"
        else
            actual_cmd="$cmd"
        fi
        missing_commands+=("$actual_cmd")
    fi
done

# Display non-installed commands and ask the user if they want to install
if [ "${#missing_commands[@]}" -gt 0 ]; then
    echo "The following packages are not available:"
    for missing_cmd in "${missing_commands[@]}"; do
        echo "- $missing_cmd"
    done

    read -p "Do you want to install the missing packages? (y/n): " choice
    if [ "$choice" == "y" ]; then
        
        # Function to install packages
        install_packages() {
            local packages=("$@")
            for package in "${packages[@]}"; do
                if conda list "$package" -n "$ENV_NAME" | grep -q "$package"; then
                    echo "Package '$package' is already installed. Skipping installation."
                else
                    echo "Installing package '$package'..."
                    conda install -y -c bioconda "$package"
                fi
            done
        }
            
        ENV_NAME="TrackTxw"

        # Check if Conda environment exists
        if conda info --envs | grep -q "$ENV_NAME"; then
            echo "Activating conda environment '$ENV_NAME'."
            # Activate Conda environment
            eval "$(conda shell.bash hook)"
            conda activate "$ENV_NAME"
            install_packages "${missing_commands[@]}"
            conda_environment_active=true
        else
            echo "Creating Conda environment '$ENV_NAME'..."
            conda create -y -n "$ENV_NAME"

            # Activate Conda environment
            eval "$(conda shell.bash hook)"
            conda activate "$ENV_NAME"

            conda_environment_active=true
        fi

        # Install packages
        install_packages "${missing_commands[@]}"
    else
        echo "No packages will be installed."
    fi
else
    echo "All required packages are available."
fi

sleep 1; clear; sleep 1

#----------------------------------------------------------------------------------------------------------------------------------------------

options=("Homo Sapiens (T2T-CHM13v2/hs1)" "Homo Sapiens (hg38)" "Mus Musculus (GRCmm39)" "Drosophila Melanogaster (dm6)" "Canis Lupus Familiaris (canFam6)" "Other" "Quit")
echo "Select the organism (and genome version) you'd like to work with: "
select opt in "${options[@]}"; do
    case "$REPLY" in
        1) organism="hs1"; break ;;
        2) organism="hg38"; break ;;
        3) organism="mm39"; break ;;
        4) organism="dm6"; break ;;
        5) organism="canFam6"; break ;;
        6) read -p "Which organism would you like to work with (answer as a genome version, eg. hs1 or canFam6)? " organism; break ;;
        7) echo;echo "Goodbye!"; exit 0; ;;
        *) echo "Invalid option. Try another one." >&2
    esac
done

echo

mkdir -p $organism 

#-----------------------------------------------------------------------

### The following parts are for the preparation of data 

#-----------------------------------------------------------------------

genome () {

    mkdir -p ${organism}/genome
    cd ${organism}/genome

    echo "Downloding the reference genome for $organism and running bowtie2-build - this might take some time!"
    # download the reference genome from NCBI, rename the file "reference_genome.fa.gz" and unzip the .gz FASTA file:
    if wget --no-check-certificate -cO genome_$organism.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/${organism}/bigZips/${organism}.fa.gz; then
        gunzip genome_$organism.fa.gz
    else
        read -p "The genome could not be found using the golden path - please enter the URL to the desired genome with the .fa.gz filetype (or press enter to use local genome): " URL
        if [ ! -z "${URL}" ]; then 
            wget --no-check-certificate -cO genome_$organism.fa.gz $URL
            gunzip genome_$organism.fa.gz 
        else 
            echo "assuming the genome is avialible in the folder named 'genome' with the name formated as 'genome_${organism}.fa'"
            sleep 5
        fi

    fi

    bowtie2-build genome_$organism.fa bowtie_$organism

    cd ../../
}

#-----------------------------------------------------------------------

sample () {
    mkdir -p ${organism}/samples

    options=("Input SRR numbers" "Load local files" "Quit")
    echo "Do you want to download data from Gene Expression Omnibus using SRR numbers or load local fasta/q files?"
    echo "Choose an option: "
    select opt in "${options[@]}"; do
        case "$REPLY" in
        1) fastq_source="Input_SRR";;
        2) fastq_source="Load_files";;
        $((${#options[@]}))) echo "Goodbye!"; exit 0; ;;
        *) echo "Invalid option. Try another one.";continue;;
        esac
        break
    done

    #-----------------------------------------------------------------------

    mkdir -p ${organism}/data

    if [ $fastq_source == "Input_SRR" ]; then     #Input your own SRR numbers from GEO
        cd ${organism}/data
    
        fasta_numbers=() # fasta_numbers is an array
        sample_names=() # sample_names is an array

        read -p "Do you have replicates that you want to concatenate? (y/n)? " answer
        case ${answer:0:1} in
            y|Y )
                read -p "How many replicates of each samples do you want to analyse?: " cvar1
                read -p "How many samples do you want to analyse?: " cvar2
                name_loop=1; name_loop2=0
                while [ $name_loop -le $cvar2 ] ; do
                    read -p "What would you like to call sample ${name_loop}? " sample_names[$name_loop2]
                    name_loop=$((name_loop+1))
                    name_loop2=$((name_loop2+1))
                done
                loop0=1 ; loop1=1 ; loop2=1
                while [ $loop0 -le $(($cvar1 * $cvar2)) ] ; do
                    while [ $loop2 -le $cvar2 ] ; do
                        while [ $loop1 -le $cvar1 ] ; do
                            read -p "Enter the SRR digits (without the letters 'SRR') of sample ${loop2} replicate ${loop1}: " fasta_numbers[$loop0]
                            loop1=$((loop1+1))
                            loop0=$((loop0+1))
                        done
                        loop1=1
                        loop2=$((loop2+1))
                    done
                done ;;
            * )
                read -p "How many samples (without replicates) do you want to analyse?: " cvar3
                name_loop=1; name_loop2=0
                while [ $name_loop -le $cvar3 ] ; do
                    read -p "What would you like to call sample ${name_loop}? " sample_names[$name_loop2]
                    name_loop=$((name_loop+1))
                    name_loop2=$((name_loop2+1))
                done
                loop3=1
                while [ $loop3 -le $cvar3 ] ; do
                    read -p "Enter the digits of the SRR #${loop3}: " fasta_numbers[$loop3]
                    loop3=$((loop3+1))
                done;;
        esac

        fasta_numbers=("${fasta_numbers[@]/#/SRR}")

        if [ $cvar1 -gt 1 ] ; then
            for nr in ${fasta_numbers[@]}; do
                fasterq-dump ${nr} -p
            done

            loop1=0; loop2=1; loop3=0;
            while [ $loop2 -le $cvar2 ] ; do
                cat ${fasta_numbers[$loop1]}.fastq ${fasta_numbers[$((loop1+1))]}.fastq > ~-/${organism}/samples/${sample_names[$loop3]}.fastq
                loop1=$((loop1+2))
                loop2=$((loop2+1))
                loop3=$((loop3+1))
            done
            echo "Samples have been downloaded and replicates have been concatinated"
        else
            loop1=1; loop2=0;
            while [ $loop1 -le ${#fasta_numbers[@]} ] ; do
                for nr in ${fasta_numbers[@]}; do
                    fasterq-dump ${nr}
                    mv ${nr}.fastq ~-/${organism}/samples/${sample_names[$loop2]}.fastq
                    loop1=$((loop1+1))
                    loop2=$((loop2+1))
                done
            done
            echo "Samples have been downloaded"
        fi

    elif [ $fastq_source == "Load_files" ]; then    #Use your own fastq-files
        cd ${organism}/data

        echo "Move your files into the folder named 'data', it's found in the folder '${organism}'"
        echo "Name the files/samples so that sample 1 contains 'sample_1', sample 2 contains 'sample_2' etc."
        echo
        read -p "Do you need to split the files based on a barcode? (y/n): " answer2
        read -p "Do you have replicates that you want to concatenate? (y/n)? " answer4
        read -p "Have you moved all the files into the specified folder? (y/n): " answer

        #--------------------------------------------------

        case ${answer4:0:1} in
                y|Y )
                    read -p "How many samples do you want to analyse?: " cvar2
                    name_loop=1; name_loop2=0
                    while [ $name_loop -le $cvar2 ] ; do
                        read -p "What would you like to call sample ${name_loop}? " sample_names[$name_loop2]
                        name_loop=$((name_loop+1))
                        name_loop2=$((name_loop2+1))
                    done;;
                * ) echo;;
        esac

        #--------------------------------------------------

        case ${answer2:0:1} in
            y|Y )
                read -p "Have you prepaired a file named 'barcodes' with the barcodes and sample names and put the txt-file in the folder 'data'? (y/n): " answer3
                case ${answer3:0:1} in
                    y|Y ) cat *.fastq | fastx_barcode_splitter.pl --bcfile barcodes.txt --bol --mismatches 2 --partial 2 --prefix ../samples/ --suffix ".fastq" --quiet;;
                    * ) echo "Prepare the neccesary files and then come back!"; exit;; 
                esac;;
            * ) echo
            #* ) cp *.fastq ~-/${organism}/samples/;; 
        esac

        #cd ../samples

        #--------------------------------------------------

        case ${answer:0:1} in
            y|Y )
                case ${answer4:0:1} in
                y|Y )
                    loop3=0
                    for ((i=1;i<=cvar2;i++)); do
                        cat *sample_${i}*.fastq > ~-/${organism}/samples/${sample_names[$loop3]}.fastq
                        loop3=$((loop3+1))
                    done
                ;;
                * )
                    loop3=0;
                    for ((i=1;i<=cvar2;i++)); do
                        cp sample_${i}.fastq ~-/${organism}/samples/${sample_names[$loop3]}.fastq
                        loop3=$((loop3+1))
                    done;;
                esac;;
            * )
                echo "move the files and come back later"
            ;;
        esac

    fi

    cd ../../
}

#-----------------------------------------------------------------------

check_quality () {

    mkdir -p ${organism}/data/quality_reports 

    sample_list=()  # Declare an empty array

    for file in ${organism}/samples/*.fastq; do
        filename=$(basename "$file")
        filename=${filename%.*}
        sample_list+=("$filename")  # Append the filename to the array
    done
    

    cd ${organism}/samples

    for x in "${sample_list[@]}"; do 

        # Output quality report
        output_report="quality_report_${x}.txt"

        # Run fastq_quality_stats
        fastx_quality_stats -i "${x}.fastq" -o "$output_report"

        echo "Quality report generated: $output_report"
        mv $output_report ../data/quality_reports 
    done

    cd ../../../
}

#-----------------------------------------------------------------------

convert () {
    read -p "What is your PCR adapter sequence? (just press ENTER if there is no adapter, type molgen for 'our' standard adapter): " adapter_sequence
    if [ -z "$adapter_sequence" ]; then
        echo 
    elif [ $adapter_sequence == "molgen" ]; then 
        adapter_sequence="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC" # for K562 and MEF
    fi
    
    read -p "Does your sample contain UMIs (y/n)?: " UMI_status
    case ${UMI_status:0:1} in
    y|Y ) UMI_status=true;;
    * ) UMI_status=false;; esac

    read -p "Does your sample contain a barcode (y/n)?: " barcode_status
    case ${barcode_status:0:1} in
    y|Y ) barcode_status=true
          read -p "How long is the barcode (how many bp)?: " barcode_length;;
    * ) barcode_status=false;; esac
    
    options=("Homo Sapiens (T2T-CHM13v2/hs1)" "Homo Sapiens (hg38)" "Mus Musculus (GRCm39)" "Drosophila Melanogaster (DM6)" "Canis Lupus Familiaris (canFam6)" "No spike in")
    PS3="Select the genome you'd like to use as 'spike in': "
    IFS=$'\n'
    select opt in "${options[@]}" "Quit"; do 
        case "$REPLY" in
        1) spikeIn_org="hs1";; 
        2) spikeIn_org="hg38";; 
        3) spikeIn_org="mm39";;
        4) spikeIn_org="dm6";; 
        5) spikeIn_org="canFam6";; 
        6) spikeIn_org=false;;
        $((${#options[@]}))) echo "Goodbye!"; exit 0; ;;
        *) echo "Invalid option. Try another one.";continue;;
        esac
        break
    done
    echo

    #-----------------------------------------------------------------------

    # Make a list of your samples, without the .fastq extension

    sample_list=()  # Declare an empty array

    for file in ${organism}/samples/*; do
        filename=$(basename "$file")
        filename=${filename%.*}
        sample_list+=("$filename")  # Append the filename to the array
    done

    cd ${organism}/genome

    if [ ! -f chrm_sizes_${organism}.txt ]; then
        if fetchChromSizes ${organism} > chrm_sizes_$organism.temp; then
            grep -v "_" chrm_sizes_$organism.temp > chrm_sizes2_$organism.temp
            #grep -v "chrM" chrm_sizes2_$organism.temp > chrm_sizes_$organism.txt
            sort -V -o chrm_sizes_${organism}.txt chrm_sizes_${organism}.txt
            rm *.temp
        else
            echo; echo "Could not find a list of chromosome sizes - create a list in txt.format with the chromosomes and the size of each chromosome, name it chrm_sizes_${organism} and put it in the 'genome' folder before you start the conversion to bigWig."
            exit
        fi
    else
        echo
    fi
    cd ../../

    #-----------------------------------------------

    project_root=$(pwd)

    mkdir -p ${organism}/bigWig
    cd ${organism}/bigWig

    #-----------------------------------------------

    #Add the path to the genome version that was previously built with bowtie2
    bowtie="${project_root}/${organism}/genome/bowtie_${organism}"
    chromSizes="${project_root}/${organism}/genome/chrm_sizes_${organism}.txt"
    spike_in_seq="${project_root}/${spikeIn_org}/genome/bowtie_${spikeIn_org}" 

    for x in "${sample_list[@]}"; do 

    # ---------------[PART 1]-----------------------------

    ### Convert the file to a fasta file
        if [[ ${project_root}/${organism}/samples/${x}.fastq == *.fastq ]]; then
            echo; echo "Converting sample ${x} to fasta file"
            sed -n '1~4s/^@/>/p;2~4p' ${project_root}/${organism}/samples/${x}.fastq > ${x}.fa; fi

    ### Remove the adapter sequence:
        if [ -n "$adapter_sequence" ]; then
            echo "removing adapter sequence from ${x}:"
            fastx_clipper -v -n -a ${adapter_sequence} -l 12 -i ${x}.fa -o ${x}_clipped.fa 
        else
            rsync -ah --progress ${x}.fa ${x}_clipped.fa; fi

    ### Remove the UMI:
        if [ "$UMI_status" = true ]; then
            echo "removing PCR duplicates and collapsing reads from ${x}:" 
            fastx_collapser -v -i ${x}_clipped.fa -o ${x}_clipped_UMI.fa 
            mv ${x}_clipped_UMI.fa ${x}_clipped.fa
        fi
    
    ### Remove the 7bp barcode:
        if [ "$barcode_status" = true ]; then    
            echo trimming the ${barcode_length} bp molecular barcodes on left from ${x} :
            fastx_trimmer -f ${barcode_length} -v -i ${x}_clipped.fa -o ${x}_clipped_barcode.fa  
            mv ${x}_clipped_barcode.fa ${x}_clipped.fa
        fi
            
    ### Create a reverese complement of the sequences:
        echo "making reverse complement of the clipped sequences of ${x}"
        fastx_reverse_complement -i ${x}_clipped.fa -o ${x}_RC.fa 
        echo

        echo "aligning reads from ${x} to the $organism genome with bowtie2"
        bowtie2 --end-to-end -p 7 -x ${bowtie} --un ${x}_hgRemovedTemp.fa -f -U ${x}_RC.fa | awk '$2 != 4 {print}' | sed '/XS:/d' | samtools view -S -b '-' > ${x}.bam
        bowtie2 --end-to-end -p 7 -x ${bowtie} --un ${x}_hgRemovedTemp.fa -f -U ${x}_RC.fa | awk '$2 != 4 {print}' | samtools view -S -b '-' > ${x}_allMap.bam
        ## aligns reads to hg using bowtie2. --end-to-end required the whole read to align. I.e. no soft clipping allowed, maintains the nucleotide-resolution of PRO-seq. -p 7 uses 7 cores, -x ${bowtie} provides the bowtie2 indexed genome, and -f -U ${x}_clipped_RC.fa provides the input file.
        ## the output sam file is piped to awk which removes the read that has no reported alignments (flag 4 in 2nd column of sam file), and the sed '/XS:/d' command gets rid of multimapped reads. --un saves unaligned reads to a separate fasta file that we will use to map the spike-in reads.
        ## samtools view command converts sam file to bam file, -S is for input is SAM, -b is for output BAM.
        echo

        fastx_clipper -l 40 -i ${x}_hgRemovedTemp.fa -o ${x}_hgRemovedTemp_min40nt.fa   ### Removes short reads before spike-in mapping. The short reads are the most like ones to cross-map to distinct genomes.

        if [ "$spikeIn_org" != false ]; then
            if [ -f "${project_root}/${spikeIn_org}/genome/bowtie_${spikeIn_org}.1.bt2" ]; then
                echo "Bowtie2 index for ${spikeIn_org} already exists, moving on..."
            else
                cd ../../

                mkdir -p ${spikeIn_org}/genome
                cd ${spikeIn_org}/genome
                echo "Downloding the reference genome for ${spikeIn_org} and running bowtie2-build - this might take some time!"
                # download the reference genome from NCBI, rename the file "reference_genome.fa.gz" and unzip the .gz FASTA file:
                wget -q --no-check-certificate -cO - https://hgdownload.soe.ucsc.edu/goldenPath/${spikeIn_org}/bigZips/${spikeIn_org}.fa.gz > genome_${spikeIn_org}.fa.gz
                
                fetchChromSizes ${spikeIn_org} > chrm_sizes_${spikeIn_org}.temp
                grep -v "_" chrm_sizes_${spikeIn_org}.temp > chrm_sizes2_${spikeIn_org}.temp
                grep -v "chrM" chrm_sizes2_${spikeIn_org}.temp > chrm_sizes_${spikeIn_org}.txt
                sort -V -o chrm_sizes_${spikeIn_org}.txt chrm_sizes_${spikeIn_org}.txt
                
                gunzip genome_${spikeIn_org}.fa.gz
                bowtie2-build genome_${spikeIn_org}.fa bowtie_${spikeIn_org}

                rm *.temp
                cd ${project_root}/${organism}/bigWig
            fi

            echo "aligning reads that did not match to the ${organism} genome to ${spikeIn_org} genome"
            bowtie2 --end-to-end -p 7 -x ${project_root}/${spikeIn_org}/genome/bowtie_${spikeIn_org} -f -U ${x}_hgRemovedTemp_min40nt.fa | awk '$2 != 4 {print}' | sed '/XS:/d' | samtools view -S -b '-' > ${x}_${spikeIn_org}_spikeIn.bam
        fi

        rm *Temp*

    # ---------------[PART 2]-----------------------------

    ### the bam file is converted to bed and piped to be sorted by chromosome and then start position.
        echo "converting bam to bed file and sorting it"
        bamToBed -i ${x}.bam > ${x}.bed 
        sort -k1,1V -k2,2n ${x}.bed -o ${x}.bed
        bamToBed -i ${x}_allMap.bam > ${x}_allMap.bed
        sort -k1,1V -k2,2n ${x}_allMap.bed -o ${x}_allMap.bed
    
        if [ "$spikeIn_org" != false ]; then 
            bamToBed -i ${x}_${spikeIn_org}_spikeIn.bam > ${x}_${spikeIn_org}_spikeIn.bed
            sort -k1,1V -k2,2n ${x}_${spikeIn_org}_spikeIn.bed -o ${x}_${spikeIn_org}_spikeIn.bed; fi

        cSample=$(grep -c ^ ${x}_allMap.bed)

        echo "total count from .bed file of ${x} is:"
        echo $cSample

        if [ "$spikeIn_org" != false ]; then
            nSpikeIn=$(grep -c ^ ${x}_${spikeIn_org}_spikeIn.bed)
            ## n is calculated by summing the reads, which is equivalent to the number of lines in the dm6 bed file. Since we only retained uniquely mapping reads (see the bowtie2 line), only reads that uniquely mapped to dm6 are counted here.
            echo; echo "Number of reads mapping to spikeIn dm6 genome is: $nSpikeIn"; fi
            # prints the number of reads that uniquely map to spike-in genome

        grep -v "_" ${x}.bed > ${x}_temp.bed 
        grep -v "chrM" ${x}_temp.bed  > ${x}.bed
        grep -v "_" ${x}_allMap.bed > ${x}_allMap_temp.bed 
        grep -v "chrM" ${x}_allMap_temp.bed  > ${x}_allMap.bed

        rm ${x}_temp.bed ${x}_allMap_temp.bed

    # ---------------[PART 3]-----------------------------
   
        echo; echo "generating non-normalized bedgraph files of ${x}. Retains only the three_prime_nt of each read"
        awk '$6 == "+"' ${x}.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizes} > ${x}_unnorm_pl.bedgraph
        awk '$6 == "-"' ${x}.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizes} > ${x}_unnorm_m.bedgraph
        awk '{$4=$4*-1; print}' ${x}_unnorm_m.bedgraph > ${x}_unnorm_mn.bedgraph    # turn the number 1 to -1
        
        awk -F'\t' '!/_/ && $1 != "chrM" {print}' ${x}_allMap.bed > ${x}_allMap_temp.bed
        mv ${x}_allMap_temp.bed ${x}_allMap.bed

        echo; echo "generating non-normalized bedgraph files of allMapping reads of ${x}. Retains only the three_prime_nt of each read"
        awk '$6 == "+"' ${x}_allMap.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizes} > ${x}_allMap_unnorm_pl.temp
        awk '$6 == "-"' ${x}_allMap.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizes} > ${x}_allMap_unnorm_m.temp
        awk '{$4=$4*-1; print}' ${x}_allMap_unnorm_m.temp > ${x}_allMap_unnorm_mn.temp    # turn the number 1 to -1

        awk -F'\t' '!/_/{print}' ${x}_allMap_unnorm_pl.temp > ${x}_allMap_unnorm_pl.bedgraph
        awk -F'\t' '!/_/{print}' ${x}_allMap_unnorm_mn.temp > ${x}_allMap_unnorm_mn.bedgraph

        echo; echo "making bigwig from non-normalized bedgraph files of ${x}"
        bedGraphToBigWig ${x}_unnorm_pl.bedgraph ${chromSizes} ${x}_unnorm_pl.bigWig
        bedGraphToBigWig ${x}_unnorm_mn.bedgraph ${chromSizes} ${x}_unnorm_mn.bigWig
        ## makes bigwig file from the bedgraph files

        echo; echo "making bigwig from non-normalized bedgraph files of allMapping reads of ${x}"
        bedGraphToBigWig ${x}_allMap_unnorm_pl.bedgraph ${chromSizes} ${x}_allMap_unnorm_pl.bigWig
        bedGraphToBigWig ${x}_allMap_unnorm_mn.bedgraph ${chromSizes} ${x}_allMap_unnorm_mn.bigWig
        ## makes bigwig file from the bedgraph files

        rm *.temp

    done

    # ----------------------------------------------------

    cd ../../
}

#-----------------------------------------------------------------------

### The following parts are for the analysis of data ###

#-----------------------------------------------------------------------

functionalgenomicregions () {

    # Check if the gene list exists
    if test -e "${organism}/genome/genes.txt"; then
        echo
    else
        echo; echo "Downloading a gene list from GitHub..."

        if [ $organism == "hs1" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_hs1.txt > ${organism}/genome/genes.txt
        elif [ $organism == "hg38" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_hg38.txt > ${organism}/genome/genes.txt
        elif [ $organism == "mm39" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_mm39.txt > ${organism}/genome/genes.txt
        elif [ $organism == "dm6" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_dm6.txt > ${organism}/genome/genes.txt
        elif [ $organism == "canFam6" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_canFam6.txt > ${organism}/genome/genes.txt
        elif [ $organism == "TAIR" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_TAIR.txt > ${organism}/genome/genes.txt
        else
            read -p "Please enter a URL to a gene list with the .gff.gz filetype: " URL_genes
            wget -q --no-check-certificate -cO genes_${organism}.gff.gz $URL_genes

            #-----------------------------------------

            gunzip genes_${organism}.gff.gz

            #-----------------------------------------

            # extract all rows containing genes:
            awk '$3 == "gene" {print $1, $4, $5, $7, $9}' genes_${organism}.gff > extracted_genes.txt

            # extract all gene names:
            sed -Ee 's/(.*Name=)([^;]*)(;gbkey.*)/\2/' -Ee 's/(.*Name=)([^;]*)(;description.*)/\2/' -Ee 's/(.*ID=gene-Dmel_)([^;]*)(;Dbxref.*)/\2/' extracted_genes.txt > modified_file.txt

            # remove all data except chr, txstart, txend and strand:
            cut -d' ' -f 1-4 extracted_genes.txt > temp.txt; paste -d' ' temp.txt modified_file.txt > genes0.temp

            #-----------------------------------------

            # extract old chromosome names as a list
            awk '{print $1}' genes0.temp > chromosomes.temp
            chromosomes_old=($(cat chromosomes.temp | tr '[:space:]' '[\n*]' | uniq))

            #-----------------------------------------

            # create list of new chromosome lists dependant on type of organism
            echo "Input the name of each chromosome for the genome ${organism}"
            chromosomes_new=()
            while IFS= read -r -p "Next item (end with an empty line): " line; do
                [[ $line ]] || break  # break if line is empty
                chromosomes_new+=("$line")
            done

            #-----------------------------------------

            # determine number of chromosomes in list
            array_length=${#chromosomes_new[@]}

            # rename chromosomes in gene-file from old to new chromosomes
            for i in $(seq 0 $((array_length-1))); do
            sed "s/${chromosomes_old[${i}]}/${chromosomes_new[${i}]}/g" genes0.temp > temporary.temp && mv temporary.temp genes0.temp
            done

            # add titles to the different columns in the file and change all seperators to tabs
            echo "chr txStart txEnd strand geneName" | cat - genes0.temp > genes0_temp.temp && mv genes0_temp.temp genes0.temp
            sed 's/ /\t/g' genes0.temp > temp.temp && mv temp.temp genes_txt.temp

            pattern=$(IFS="|"; echo "${chromosomes_new[*]}")
            grep -E "$pattern" genes_txt.temp > genes.txt
            mv genes.txt ${organism}/genome/

            # remove temporary files
            rm *.temp genes_${organism}.gff temp.txt extracted_genes.txt modified_file.txt 
        fi

        #-----------------------------------------

        echo "List of genes created.."
        echo
    fi

    #-----------------------------------------------------------------------

    #Make a list of your samples, without the .fastq extension

    cd ${organism}/samples
    declare -a sample_list
    for file in *.fastq; do
        sample_list=("${sample_list[@]}" "$file")
    done
    sample_list=("${sample_list[@]%.*}")
    cd ../../

    #---------------------------------------------------

    for sample in ${sample_list[@]}; do 

        echo "Creating functional genomics .bed-file for sample $sample"

        mkdir -p ${organism}/analysis/functionalGenomics_${sample}          # create a folder for output files

        #---------------------------------------------------

        wget -q -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/scripts/R1.R > R1.R
        wget -q -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/scripts/R2.R > R2.R

        Rscript R1.R $organism $sample --save   # R_script that adds coordinates for different features such as TSS, CPS, genebody etc.

        awk '{$4=$4*-1; print}' ${organism}/bigWig/${sample}_allMap_unnorm_mn.bedgraph > ${organism}/bigWig/${sample}_allMap_unnorm_mn_temp.bedgraph    # turn the number 1 to -1
        cat ${organism}/bigWig/${sample}_allMap_unnorm_mn_temp.bedgraph | tr ' ' '\t' > ${organism}/bigWig/${sample}_allMap_unnorm_mn_temp2.bedgraph

        bedtools slop -i ${organism}/bigWig/${sample}_allMap_unnorm_pl.bedgraph -g ${organism}/genome/chrm_sizes_${organism}.txt -b 10 | bedtools merge -d 1 -c 4 -o max > ${organism}/analysis/functionalGenomics_${sample}/pl_bedgraph.bed
        bedtools slop -i ${organism}/bigWig/${sample}_allMap_unnorm_mn_temp2.bedgraph -g ${organism}/genome/chrm_sizes_${organism}.txt -b 10 | bedtools merge -d 1 -c 4 -o max > ${organism}/analysis/functionalGenomics_${sample}/mn_bedgraph.bed
        bedtools intersect -a ${organism}/analysis/functionalGenomics_${sample}/pl_bedgraph.bed -b ${organism}/analysis/functionalGenomics_${sample}/mn_bedgraph.bed > ${organism}/analysis/functionalGenomics_${sample}/divergent_transcription.bed

        rm ${organism}/bigWig/${sample}_allMap_unnorm_mn_temp.bedgraph ${organism}/bigWig/${sample}_allMap_unnorm_mn_temp2.bedgraph

        #-----------------------------------------------

        bedtools intersect -v -a ${organism}/analysis/functionalGenomics_${sample}/divergent_transcription.bed -b ${organism}/analysis/functionalGenomics_${sample}/refGenes_TSSpm500.txt > ${organism}/analysis/functionalGenomics_${sample}/enhancers.bed
        bedtools intersect -wa -a ${organism}/analysis/functionalGenomics_${sample}/refGenes_TSSpm500.txt -b ${organism}/analysis/functionalGenomics_${sample}/divergent_transcription.bed > ${organism}/analysis/functionalGenomics_${sample}/activeGenes.bed
        Rscript R2.R $organism $sample --save

        rm R1.R R2.R

        #-----------------------------------------------

        echo; echo "retaining 3' most coordinate of the bed file."

        cd ${organism}/analysis/functionalGenomics_${sample}

        project_root="$(dirname $(dirname $(dirname $(realpath $0))) )"
        awk '$6 == "+"' ${project_root}/bigWig/${sample}.bed > tempPL.bed
        awk '$6 == "-"' ${project_root}/bigWig/${sample}.bed > tempMN.bed
        awk '{$2 = $3; print}' tempPL.bed > tempPL_3p.bed
        awk '{$3 = $2; print}' tempMN.bed > tempMN_3p.bed

        cat tempPL_3p.bed tempMN_3p.bed | tr ' '  '\t' > PROseq_3pnt.bed

        rm *temp*

        #-----------------------------------------------

        echo; echo "creating .bed files"

        ##counting engaged Pol II at promoter-proximal regions
        bedtools intersect -s -u -wa -a PROseq_3pnt.bed -b ppPolII.txt > PROseq_ppPolII.bed
        bedtools intersect -v -a     PROseq_3pnt.bed -b ppPolII.txt > ppRemoved.bed

        ##counting engaged Pol II at the sites of divergent transcription
        bedtools intersect -S -u -wa -a ppRemoved.bed -b divTx.txt > PROseq_ppDiv.bed
        bedtools intersect -v -a     ppRemoved.bed -b divTx.txt > ppdivRemoved.bed

        ##counting engaged Pol II at enhancers
        bedtools intersect -u -wa -a ppdivRemoved.bed -b enhancers.bed > PROseq_enhancers.bed
        bedtools intersect -v -a     ppdivRemoved.bed -b enhancers.bed > ppdivEnhRemoved.bed

        ##counting engaged Pol II at CPS
        bedtools intersect -s -u -wa -a ppdivEnhRemoved.bed -b CPS.txt > PROseq_CPS.bed
        bedtools intersect -v -a     ppdivEnhRemoved.bed -b CPS.txt > ppdivEnhCPSRemoved.bed

        ##counting engaged Pol II at GB
        bedtools intersect -s -u -wa -a ppdivEnhCPSRemoved.bed -b geneBody.txt > PROseq_GB.bed
        bedtools intersect -v -a     ppdivEnhCPSRemoved.bed -b geneBody.txt > ppdivEnhCPSgbRemoved.bed

        ##counting engaged Pol II at short genes
        bedtools intersect -s -u -wa -a ppdivEnhCPSgbRemoved.bed -b shortGenes.txt > PROseq_SG.bed
        bedtools intersect -v -a     ppdivEnhCPSgbRemoved.bed -b shortGenes.txt > ppdivEnhCPSgbSGRemoved.bed

        ##counting engaged Pol II at termination windows
        bedtools intersect -s -u -wa -a ppdivEnhCPSgbSGRemoved.bed -b TW.txt > PROseq_TW.bed
        bedtools intersect -v -a     ppdivEnhCPSgbRemoved.bed -b TW.txt > PROseq_noGene_noEnh.bed

        rm *Removed.bed

        echo; echo ".bed files created"

        #-----------------------------------------------

        echo; echo "here comes count data:"

        ## Total count of active sites of transcription (uniquely mapping reads) in the dataset:
        echo total reads:
        totC=$(awk 'END { print NR }' PROseq_3pnt.bed)
        echo ${totC}
        ## Engaged Pol II at transcribed enhancers:
        echo enhancers:
        enhC=$(awk 'END { print NR }' PROseq_enhancers.bed)
        echo ${enhC}

        ## Engaged Pol II at divergent transcripts:
        echo divergentTx:
        divC=$(awk 'END { print NR }' PROseq_ppDiv.bed)
        echo ${divC}

        ## Engaged Pol II at promoters:
        echo promoter:
        ppC=$(awk 'END { print NR }' PROseq_ppPolII.bed)
        echo ${ppC}

        ## Engaged Pol II at gene bodies:
        echo gene body:
        GBC=$(awk 'END { print NR }' PROseq_GB.bed)
        echo ${GBC}

        ## Engaged Pol II at the cleavage and polyadenylation regions:
        echo CPS:
        CPSC=$(awk 'END { print NR }' PROseq_CPS.bed)
        echo ${CPSC}

        ## Engaged Pol II at short genes:
        echo short genes:
        SGC=$(awk 'END { print NR }' PROseq_SG.bed)
        echo ${SGC}

        ## Engaged Pol II at the termination windows:
        echo termination window:
        twC=$(awk 'END { print NR }' PROseq_TW.bed)
        echo ${twC}

        ## Engaged Pol II at sites that did not localize to any of the identified functional regions:
        echo noGenesNoEnhancers:
        noGdnoE=$(awk 'END { print NR }' PROseq_noGene_noEnh.bed)
        echo ${noGdnoE}

        #-----------------------------------------------

        echo; echo "creating a bed-file named functionalGenomicRegions for further analysis"

        bedtools intersect -c -wa -a ppPolII.txt -b PROseq_ppPolII.bed > ppPolCounts.tmp
        bedtools intersect -c -wa -a divTx.txt -b PROseq_ppDiv.bed > ppDivCounts.tmp
        bedtools intersect -c -wa -a enhancers.bed -b PROseq_enhancers.bed > enhancerCounts.tmp
        bedtools intersect -c -wa -a geneBody.txt -b PROseq_GB.bed > geneBodyCounts.tmp
        bedtools intersect -c -wa -a CPS.txt -b PROseq_CPS.bed > CPSCounts.tmp
        bedtools intersect -c -wa -a shortGenes.txt -b PROseq_SG.bed > SGCounts.tmp
        bedtools intersect -c -wa -a TW.txt -b PROseq_TW.bed > TerminationWinCounts.tmp

        cut -f 1-4,6- -d " " enhancerCounts.tmp > enhancerCounts.tmp

        #-----------------------------------------------

        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="243,132,0"; print }' ppPolCounts.tmp > ppPolCounts.bed                         # orange
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="178,59,212"; print }' ppDivCounts.tmp > ppDivCounts.bed                        # purple
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="115,212,122"; print }' enhancerCounts.tmp > enhancerCounts.bed                 # green
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="0,0,0"; print }' geneBodyCounts.tmp > geneBodyCounts.bed                       # black
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="103,200,249"; print }' CPSCounts.tmp > CPSCounts.bed                           # blue
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="253,218,13"; print }' SGCounts.tmp > SGCounts.bed                              # yellow
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="255,54,98"; print }' TerminationWinCounts.tmp > TerminationWinCounts.bed       # red

        #-----------------------------------------------

        cat ppPolCounts.bed ppDivCounts.bed enhancerCounts.bed geneBodyCounts.bed CPSCounts.bed SGCounts.bed TerminationWinCounts.bed > catRegions.temp

        awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $7 "\t" $8}' catRegions.temp > catRegions2.temp
        awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $4 "\t" $5 "\t" $2 "\t" $3 "\t" $7}' catRegions2.temp > catRegions3.temp
        awk '!seen[$1,$2,$3]++' catRegions3.temp | sortBed > catRegions4.temp
        #awk -F '\t' -v OFS='\t' '{ $(NF+1) ="."; print }' catRegions3.temp > catRegions4.temp

        touch headerLine.txt
        echo track name="functional_genomic_regions_$sample" itemRgb="On" >> headerLine.txt

        cat headerLine.txt catRegions4.temp > functionalGenomicRegions_${sample}.bed
        awk -v c=${totC} '{ $4 /= c; print }' OFS='\t' functionalGenomicRegions_${sample}.bed > functionalGenomicRegions_${sample}_norm.bed
        awk -i inplace -v OFS='\t' '{ $4 *= 1000000; print }' functionalGenomicRegions_${sample}_norm.bed

        mv functionalGenomicRegions_${sample}.bed ${project_root}/analysis/functionalGenomicRegions_${sample}.bed
        mv functionalGenomicRegions_${sample}_norm.bed ${project_root}/analysis/functionalGenomicRegions_${sample}_norm.bed

        #-----------------------------------------------

        rm *.temp
        rm *.tmp
        rm headerLine.txt
        cd ../../../

        #-----------------------------------------------

    done

    echo "Analysis of functional genomic regions done"
    echo

    #---------------------------------------------------

    cd ../../
}

#-----------------------------------------------------------------------

CompGeneExpression () {

    mkdir -p ${organism}/analysis/differential_expression
    cd ${organism}/samples
    declare -a sample_list
    for file in *.fastq; do
       sample_list=("${sample_list[@]}" "$file")
    done
    sample_list=("${sample_list[@]%.*}")
    cd ../../

    #---------------------------------------------------

    read -p "You need at least three samples to use DEseq2 to determine which genes are up/down regulated. Do you want to continue? (y/n) " continue_DEseq2
    case ${continue_DEseq2:0:1} in
    y|Y ) continue_DEseq2=true;;
    * ) echo "Thank you for using TrackTx"; break 2;; esac
    
    if [ "$continue_DEseq2" = true ]; then

        wget -q -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/scripts/R4.R > R4.R
        wget -q -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/scripts/R5.R > R5.R

        for sample in ${sample_list[@]}; do 
            Rscript R4.R $organism $sample  --save
        done

        # Initialize an empty array for user input
        condition_list=()

        # Loop to read user input until the new list has the same length as sample_list
        for item in "${sample_list[@]}"; do
            read -p "Enter condition for sample $item (eg. "control" or "treatment"): " user_input
            condition_list+=("$user_input")
        done

        # Print the user-inputted list
        echo "User-inputted list: ${condition_list[@]}"

        number_of_samples=${#sample_list[@]}

        Rscript R5.R $organism $number_of_samples ${sample_list[@]} ${condition_list[@]}

        rm R4.R R5.R
    fi
}

#-----------------------------------------------------------------------

map_repeats () {
    cd ${organism}/samples
    declare -a sample_list
    for file in *.fastq; do
        sample_list=("${sample_list[@]}" "$file")
    done
    sample_list=("${sample_list[@]%.*}")
    cd ../../

    mkdir -p ${organism}/analysis/sat_III
    cd ${organism}/analysis/sat_III

    project_root="$(dirname $(dirname $(dirname $(dirname $(realpath $0)))) )"
    bowtie="${project_root}/${organism}/genome/bowtie_${organism}"
    chromSizes="${project_root}/${organism}/genome/chrm_sizes_${organism}.txt"

    options=("A termination motif ('CAACCCGAGT' / 'CAACACGAGT')" "A twice repeating motif (CATTxCATT)" "A three times repeating motif (CATTxCATTxCATT)" "Choose the motif on your own" "Quit")
    
    echo "What type of repeat do you want to analyse?"
    PS3="Pick an option: "
    IFS=$'\n'
    select opt in "${options[@]}"; do 
        case "$REPLY" in
        1) motif="CAACCCGAGT|CAACACGAGT" ;;
        2) motif="CATT.CATT";;
        3) motif="CATT.CATT.CATT";;
        4) motif="own";;
        $((${#options[@]}))) echo "Goodbye!"; exit 0; ;;
        *) echo "Invalid option. Try another one.";continue;;
        esac
        break
    done

    #-----------------------------------------------------------------------

    if [ $motif == "own" ]; then
        read -p "Type your motif here as a regex pattern (eg: 'TATA' or 'CCC.TTT'): " motif; fi

    for x in ${sample_list[@]}; do 

        grep -E -A 2 -B 1 "$motif" ${project_root}/${organism}/samples/${x}.fastq > output_${x}.fastq
        grep -v -- "--" output_${x}.fastq > ${x}.fastq
        echo; echo "filtering done"

        sample_count=$(grep -c "^@" ${x}.fastq)
        echo; echo "Number of matches in ${x}: $sample_count"

        echo; echo "aligning sat III repeats to the hs1 genome with bowtie2"
        bowtie2 --end-to-end -p 7 -x ${bowtie} -q -U ${x}.fastq | awk '$2 != 4 {print}' | samtools view -S -b '-' > ${x}_sat_III_allMap.bam

        echo; echo "converting bam files to bed files"
        bamToBed -i ${x}_sat_III_allMap.bam > ${x}_sat_III_allMap.bed && sort -k1,1V -k2,2n ${x}_sat_III_allMap.bed -o ${x}_sat_III_allMap.bed

        cSatIII=$(grep -c ^ ${x}_sat_III_allMap.bed)
        echo "total number of sat III repeats mapped is: $cSatIII"

        echo; echo "generating non-normalized bedgraph files of allMapping reads. Retains only the three_prime_nt of each read:"
        awk '$6 == "+"' ${x}_sat_III_allMap.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizes} > ${x}_sat_III_allMap_unnorm_pl.bedgraph
        awk '$6 == "-"' ${x}_sat_III_allMap.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizes} > ${x}_sat_III_allMap_unnorm_m.bedgraph
        awk '{$4=$4*-1; print}' ${x}_sat_III_allMap_unnorm_m.bedgraph > ${x}_sat_III_allMap_unnorm_mn.bedgraph

        echo; echo "making bigWig files from non-normalized bedgraph files of allMapping reads"
        bedGraphToBigWig ${x}_sat_III_allMap_unnorm_pl.bedgraph ${chromSizes} ${x}_sat_III_unnorm_pl.bigWig
        bedGraphToBigWig ${x}_sat_III_allMap_unnorm_mn.bedgraph ${chromSizes} ${x}_sat_III_unnorm_mn.bigWig
    
    done
    
    echo
}

#-----------------------------------------------------------------------

download_and_convert () {
    mkdir -p ${organism}/genome
    mkdir -p ${organism}/samples
    mkdir -p ${organism}/data
    mkdir -p ${organism}/data/quality_reports
    mkdir -p ${organism}/analysis

### options

    #-------------------------------------------

    options=("Input SRR numbers" "Load local files" "Quit")
    echo "Do you want to download data from Gene Expression Omnibus using SRR numbers or load local fasta/q files?"
    echo "Choose an option: "
    select opt in "${options[@]}"; do
        case "$REPLY" in
        1) fastq_source="Input_SRR";;
        2) fastq_source="Load_files";;
        $((${#options[@]}))) echo "Goodbye!"; exit 0; ;;
        *) echo "Invalid option. Try another one.";continue;;
        esac
        break
    done

    #-------------------------------------------

    if [ $fastq_source == "Input_SRR" ]; then     #Input your own SRR numbers from GEO
    
        fasta_numbers=() # fasta_numbers is an array
        sample_names=() # sample_names is an array

        read -p "Do you have replicates that you want to concatenate? (y/n)? " answer
        case ${answer:0:1} in
            y|Y )
                read -p "How many replicates of each samples do you want to analyse?: " cvar1
                read -p "How many samples do you want to analyse?: " cvar2
                name_loop=1; name_loop2=0
                while [ $name_loop -le $cvar2 ] ; do
                    read -p "What would you like to call sample ${name_loop}? " sample_names[$name_loop2]
                    name_loop=$((name_loop+1))
                    name_loop2=$((name_loop2+1))
                done
                loop0=1 ; loop1=1 ; loop2=1
                while [ $loop0 -le $(($cvar1 * $cvar2)) ] ; do
                    while [ $loop2 -le $cvar2 ] ; do
                        while [ $loop1 -le $cvar1 ] ; do
                            read -p "Enter the SRR digits (without the letters 'SRR') of sample ${loop2} replicate ${loop1}: " fasta_numbers[$loop0]
                            loop1=$((loop1+1))
                            loop0=$((loop0+1))
                        done
                        loop1=1
                        loop2=$((loop2+1))
                    done
                done ;;
            * )
                read -p "How many samples (without replicates) do you want to analyse?: " cvar3
                name_loop=1; name_loop2=0
                while [ $name_loop -le $cvar3 ] ; do
                    read -p "What would you like to call sample ${name_loop}? " sample_names[$name_loop2]
                    name_loop=$((name_loop+1))
                    name_loop2=$((name_loop2+1))
                done
                loop3=1
                while [ $loop3 -le $cvar3 ] ; do
                    read -p "Enter the digits of the SRR #${loop3}: " fasta_numbers[$loop3]
                    loop3=$((loop3+1))
                done;;
        esac

        fasta_numbers=("${fasta_numbers[@]/#/SRR}")

    elif [ $fastq_source == "Load_files" ]; then    #Use your own fastq-files
        echo "Move your files into the folder named 'data', it's found in the folder '${organism}'"
        echo "Name the files/samples so that sample 1 contains 'sample1', sample 2 contains 'sample2' etc."
        echo
        read -p "Do you need to split the files based on a barcode? (y/n): " answer2
        read -p "Do you have replicates that you want to concatenate? (y/n)? " answer4

        case ${answer4:0:1} in
                y|Y )
                    read -p "How many samples do you want to analyse?: " cvar2
                    name_loop=1; name_loop2=0
                    while [ $name_loop -le $cvar2 ] ; do
                        read -p "What would you like to call sample ${name_loop}? " sample_names[$name_loop2]
                        name_loop=$((name_loop+1))
                        name_loop2=$((name_loop2+1))
                    done;;
                * ) echo;;
        esac

        read -p "Have you moved all the files into the specified folder? (y/n): " answer
    fi

    #-------------------------------------------

    read -p "Do you want to perform a quality control on the downloaded data? (y/n): " ans_QC

    #-------------------------------------------


    read -p "What is your adapter sequence? (just press ENTER if there is no adapter, type molgen for standard adapter): " adapter_sequence
    if [ -z "$adapter_sequence" ]; then
        echo 
    elif [ $adapter_sequence == "molgen" ]; then 
        adapter_sequence="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC" # for K562 and MEF
    fi
    
    read -p "Does your sample contain UMIs (y/n)?: " UMI_status
    case ${UMI_status:0:1} in
    y|Y ) UMI_status=true;;
    * ) UMI_status=false;; esac

    read -p "Does your sample contain a barcode (y/n)?: " barcode_status
    case ${barcode_status:0:1} in
    y|Y ) barcode_status=true
          read -p "How long is the barcode (how many bp)?: " barcode_length;;
    * ) barcode_status=false;; esac
    
    options=("Homo Sapiens (T2T-CHM13v2/hs1)" "Homo Sapiens (hg38)" "Mus Musculus (GRCm39)" "Drosophila Melanogaster (DM6)" "Canis Lupus Familiaris (canFam6)" "No spike in")
    PS3="Select the genome you'd like to use as 'spike in': "
    IFS=$'\n'
    select opt in "${options[@]}" "Quit"; do 
        case "$REPLY" in
        1) spikeIn_org="hs1";; 2) spikeIn_org="hg38";; 3) spikeIn_org="mm39";;
        4) spikeIn_org="dm6";; 5) spikeIn_org="canFam6";; 6) spikeIn_org=false;;
        $((${#options[@]}+1))) echo "Goodbye!"; exit 0; ;;
        *) echo "Invalid option. Try another one.";continue;;
        esac
        break
    done
    echo

    #-------------------------------------------

### bowtie2

    cd ${organism}/genome


    if [ -f "bowtie_$organism.1.bt2" ]; then
        echo "Bowtie2 index already exists, moving on..."
    else
        echo "Downloding the reference genome for $organism and running bowtie2-build - this might take some time!"
        # download the reference genome from NCBI, rename the file "reference_genome.fa.gz" and unzip the .gz FASTA file:
        if wget --no-check-certificate -cO genome_$organism.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/${organism}/bigZips/${organism}.fa.gz; then
            gunzip genome_$organism.fa.gz
        else
            read -p "The genome could not be found using the golden path - please enter the URL to the desired genome with the .fa.gz filetype (or press enter to use local genome): " URL
            if [ ! -z "${URL}" ]; then 
                wget --no-check-certificate -cO genome_$organism.fa.gz $URL
                gunzip genome_$organism.fa.gz 
            else 
                echo "assuming the genome is avialible in the folder named 'genome' with the name formated as 'genome_${organism}.fa'"
                sleep 5
            fi

        fi
        bowtie2-build genome_$organism.fa bowtie_$organism
    fi

    cd ../../

    #-------------------------------------------

### download samples

    if [ $fastq_source == "Input_SRR" ]; then     #Input your own SRR numbers from GEO
        cd ${organism}/data

        if [ $cvar1 -gt 1 ] ; then
            for nr in ${fasta_numbers[@]}; do
                fasterq-dump ${nr} -p
            done

            loop1=0; loop2=1; loop3=0;
            while [ $loop2 -le $cvar2 ] ; do
                cat ${fasta_numbers[$loop1]}.fastq ${fasta_numbers[$((loop1+1))]}.fastq > ~-/${organism}/samples/${sample_names[$loop3]}.fastq
                loop1=$((loop1+2))
                loop2=$((loop2+1))
                loop3=$((loop3+1))
            done
            echo "Samples have been downloaded and replicates have been concatinated"
        else
            loop1=1; loop2=0;
            while [ $loop1 -le ${#fasta_numbers[@]} ] ; do
                for nr in ${fasta_numbers[@]}; do
                    fasterq-dump ${nr}
                    mv ${nr}.fastq ~-/${organism}/samples/${sample_names[$loop2]}.fastq
                    loop1=$((loop1+1))
                    loop2=$((loop2+1))
                done
            done
            echo "Samples have been downloaded"
        fi

    #--------------------------------------------------

    elif [ $fastq_source == "Load_files" ]; then    #Use your own fastq-files
        cd ${organism}/data

        case ${answer2:0:1} in
            y|Y ) cat *.fastq | fastx_barcode_splitter.pl --bcfile barcodes.txt --bol --mismatches 2 --partial 2 --prefix ../samples/ --suffix ".fastq" --quiet;;
            * ) mv *.fastq ~-/${organism}/samples/;; 
        esac

        cd ../samples

        #--------------------------------------------------

        case ${answer:0:1} in
            y|Y )
                case ${answer4:0:1} in
                y|Y )
                    loop3=0
                    for ((i=1;i<=cvar2;i++)); do
                        cat *sample_${i}*.fastq > ~-/${organism}/samples/${sample_names[$loop3]}.fastq
                        loop3=$((loop3+1))
                    done
                ;;
                * )
                    loop3=0;
                    for ((i=1;i<=cvar2;i++)); do
                        cp sample_${i}.fastq ~-/${organism}/samples/${sample_names[$loop3]}.fastq
                        loop3=$((loop3+1))
                    done;;
                esac;;
            * )
                echo "move the files and come back later"
            ;;
        esac
    fi

    cd ../../

    #-----------------------------------------------------------------------

    case ${ans_QC:0:1} in
        y|Y )
            # Make a list of your samples, without the .fastq extension

            sample_list=()  # Declare an empty array

            for file in ${organism}/samples/*.fastq; do
                filename=$(basename "$file")
                filename=${filename%.*}
                sample_list+=("$filename")  # Append the filename to the array
            done

            cd ${organism}/data/quality_reports

            for x in "${sample_list[@]}"; do 

                # Output quality report
                output_report="quality_report_${x}.txt"

                # Run fastq_quality_stats
                fastx_quality_stats -i "../../samples/${x}.fastq" -o "$output_report"

                echo "Quality report generated: $output_report"
                mv $output_report ../data/quality_reports 
                
            done 
            cd ../../../
                ;;

        * ) ;; 
    esac

    #-----------------------------------------------------------------------

### convert to bigWig

    # Make a list of your samples, without the .fastq extension

    sample_list=()  # Declare an empty array

    for file in ${organism}/samples/*; do
        filename=$(basename "$file")
        filename=${filename%.*}
        sample_list+=("$filename")  # Append the filename to the array
    done

    cd ${organism}/genome

    if [ ! -f chrm_sizes_${organism}.txt ]; then
        if fetchChromSizes ${organism} > chrm_sizes_$organism.temp; then
            grep -v "_" chrm_sizes_$organism.temp > chrm_sizes2_$organism.temp
            grep -v "chrM" chrm_sizes2_$organism.temp > chrm_sizes_$organism.txt
            sort -V -o chrm_sizes_${organism}.txt chrm_sizes_${organism}.txt
            rm *.temp
        else
            echo; echo "Could not find a list of chromosome sizes - create a list in txt.format with the chromosomes and the size of each chromosome, name it chrm_sizes_${organism} and put it in the 'genome' folder before you start the conversion to bigWig."
            exit
        fi
    else
        echo
    fi
    cd ../../

    #-----------------------------------------------

    project_root=$(pwd)

    mkdir -p ${organism}/bigWig
    cd ${organism}/bigWig

    #-----------------------------------------------

    #Add the path to the genome version that was previously built with bowtie2
    bowtie="${project_root}/${organism}/genome/bowtie_${organism}"
    chromSizes="${project_root}/${organism}/genome/chrm_sizes_${organism}.txt"
    spike_in_seq="${project_root}/${spikeIn_org}/genome/bowtie_${spikeIn_org}" 

    for x in "${sample_list[@]}"; do 

    # ---------------[PART 1]-----------------------------

    ### Convert the file to a fasta file
        if [[ ${project_root}/${organism}/samples/${x}.fastq == *.fastq ]]; then
            echo; echo "Converting the samples to .fa files"
            sed -n '1~4s/^@/>/p;2~4p' ${project_root}/${organism}/samples/${x}.fastq > ${x}.fa; fi

    ### Remove the adapter sequence:
        if [ -n "$adapter_sequence" ]; then
            echo "removing adapter sequence from ${x}:"
            fastx_clipper -v -n -a ${adapter_sequence} -l 12 -i ${x}.fa -o ${x}_clipped.fa 
        else
            rsync -ah --progress ${x}.fa ${x}_clipped.fa; fi

    ### Remove the UMI:
        if [ "$UMI_status" = true ]; then
            echo "removing PCR duplicates and collapsing reads from ${x}:" 
            fastx_collapser -v -i ${x}_clipped.fa -o ${x}_clipped_UMI.fa 
            mv ${x}_clipped_UMI.fa ${x}_clipped.fa
        fi
    
    ### Remove the 7bp barcode:
        if [ "$barcode_status" = true ]; then    
            echo trimming the ${barcode_length} bp molecular barcodes on left from ${x} :
            fastx_trimmer -f ${barcode_length} -v -i ${x}_clipped.fa -o ${x}_clipped_barcode.fa  
            mv ${x}_clipped_barcode.fa ${x}_clipped.fa
        fi
            
    ### Create a reverese complement of the sequences:
        echo "making reverse complement of the clipped sequences of ${x}"
        fastx_reverse_complement -i ${x}_clipped.fa -o ${x}_RC.fa 
        echo

        echo "aligning reads from ${x} to the $organism genome with bowtie2"
        bowtie2 --end-to-end -p 7 -x ${bowtie} --un ${x}_hgRemovedTemp.fa -f -U ${x}_RC.fa | awk '$2 != 4 {print}' | sed '/XS:/d' | samtools view -S -b '-' > ${x}.bam
        bowtie2 --end-to-end -p 7 -x ${bowtie} --un ${x}_hgRemovedTemp.fa -f -U ${x}_RC.fa | awk '$2 != 4 {print}' | samtools view -S -b '-' > ${x}_allMap.bam
        ## aligns reads to hg using bowtie2. --end-to-end required the whole read to align. I.e. no soft clipping allowed, maintains the nucleotide-resolution of PRO-seq. -p 7 uses 7 cores, -x ${bowtie} provides the bowtie2 indexed genome, and -f -U ${x}_clipped_RC.fa provides the input file.
        ## the output sam file is piped to awk which removes the read that has no reported alignments (flag 4 in 2nd column of sam file), and the sed '/XS:/d' command gets rid of multimapped reads. --un saves unaligned reads to a separate fasta file that we will use to map the spike-in reads.
        ## Samtools view command converts sam file to bam file, -S is for input is SAM, -b is for output BAM.
        echo

        fastx_clipper -l 40 -i ${x}_hgRemovedTemp.fa -o ${x}_hgRemovedTemp_min40nt.fa   ### Removes short reads before spike-in mapping. The short reads are the most like ones to cross-map to distinct genomes.

        if [ "$spikeIn_org" != false ]; then
            if [ -f "${project_root}/${spikeIn_org}/genome/bowtie_${spikeIn_org}.1.bt2" ]; then
                echo "Bowtie2 index for ${spikeIn_org} already exists, moving on..."
            else
                cd ../../

                mkdir -p ${spikeIn_org}/genome
                cd ${spikeIn_org}/genome
                echo "Downloding the reference genome for ${spikeIn_org} and running bowtie2-build - this might take some time!"
                # download the reference genome from NCBI, rename the file "reference_genome.fa.gz" and unzip the .gz FASTA file:
                wget -q --no-check-certificate -cO - https://hgdownload.soe.ucsc.edu/goldenPath/${spikeIn_org}/bigZips/${spikeIn_org}.fa.gz > genome_${spikeIn_org}.fa.gz
                
                fetchChromSizes ${spikeIn_org} > chrm_sizes_${spikeIn_org}.temp
                grep -v "_" chrm_sizes_${spikeIn_org}.temp > chrm_sizes2_${spikeIn_org}.temp
                grep -v "chrM" chrm_sizes2_${spikeIn_org}.temp > chrm_sizes_${spikeIn_org}.txt
                sort -V -o chrm_sizes_${spikeIn_org}.txt chrm_sizes_${spikeIn_org}.txt
                
                gunzip genome_${spikeIn_org}.fa.gz
                bowtie2-build genome_${spikeIn_org}.fa bowtie_${spikeIn_org}

                rm *.temp
                cd ${project_root}/${organism}/bigWig
            fi
        fi

        rm *Temp*

    # ---------------[PART 2]-----------------------------

    ### the bam file is converted to bed and piped to be sorted by chromosome and then start position.
        echo "converting bam to bed file and sorting it"
        bamToBed -i ${x}.bam > ${x}.bed 
        sort -k1,1V -k2,2n ${x}.bed -o ${x}.bed
        bamToBed -i ${x}_allMap.bam > ${x}_allMap.bed
        sort -k1,1V -k2,2n ${x}_allMap.bed -o ${x}_allMap.bed
    
        if [ "$spikeIn_org" != false ]; then 
            bamToBed -i ${x}_${spikeIn_org}_spikeIn.bam > ${x}_${spikeIn_org}_spikeIn.bed
            sort -k1,1V -k2,2n ${x}_${spikeIn_org}_spikeIn.bed -o ${x}_${spikeIn_org}_spikeIn.bed; fi



        if [ "$spikeIn_org" != false ]; then
            nSpikeIn=$(grep -c ^ ${x}_${spikeIn_org}_spikeIn.bed)
            ## n is calculated by summing the reads, which is equivalent to the number of lines in the dm6 bed file. Since we only retained uniquely mapping reads (see the bowtie2 line), only reads that uniquely mapped to dm6 are counted here.
            echo; echo "Number of reads mapping to spikeIn dm6 genome is: $nSpikeIn"; fi
            # prints the number of reads that uniquely map to spike-in genome

    # ---------------[PART 3]-----------------------------
   
        echo; echo "generating non-normalized bedgraph files of ${x}. Retains only the three_prime_nt of each read"
        awk '$6 == "+"' ${x}.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizes} > ${x}_unnorm_pl.bedgraph
        awk '$6 == "-"' ${x}.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizes} > ${x}_unnorm_m.bedgraph
        awk '{$4=$4*-1; print}' ${x}_unnorm_m.bedgraph > ${x}_unnorm_mn.bedgraph    # turn the number 1 to -1

        echo; echo "generating non-normalized bedgraph files of allMapping reads of ${x}. Retains only the three_prime_nt of each read"
        awk '$6 == "+"' ${x}_allMap.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizes} > ${x}_allMap_unnorm_pl.temp
        awk '$6 == "-"' ${x}_allMap.bed | genomeCoverageBed -i stdin -3 -bg -g ${chromSizes} > ${x}_allMap_unnorm_m.temp
        awk '{$4=$4*-1; print}' ${x}_allMap_unnorm_m.temp > ${x}_allMap_unnorm_mn.temp    # turn the number 1 to -1

        awk -F'\t' '!/_/{print}' ${x}_allMap_unnorm_pl.temp > ${x}_allMap_unnorm_pl.bedgraph
        awk -F'\t' '!/_/{print}' ${x}_allMap_unnorm_mn.temp > ${x}_allMap_unnorm_mn.bedgraph

        echo; echo "making bigwig from non-normalized bedgraph files of ${x}"
        bedGraphToBigWig ${x}_unnorm_pl.bedgraph ${chromSizes} ${x}_unnorm_pl.bigWig
        bedGraphToBigWig ${x}_unnorm_mn.bedgraph ${chromSizes} ${x}_unnorm_mn.bigWig
        ## makes bigwig file from the bedgraph files

        echo; echo "making bigwig from non-normalized bedgraph files of allMapping reads of ${x}"
        bedGraphToBigWig ${x}_allMap_unnorm_pl.bedgraph ${chromSizes} ${x}_allMap_unnorm_pl.bigWig
        bedGraphToBigWig ${x}_allMap_unnorm_mn.bedgraph ${chromSizes} ${x}_allMap_unnorm_mn.bigWig
        ## makes bigwig file from the bedgraph files

        cSample=$(grep -c ^ ${x}.bed)
        ## cSample is the total number of uniquely mappig reads
    
        echo $cSample | awk '{ c="'$cSample'"; printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*1000000)/c}' ${x}_unnorm_pl.bedgraph  > ${x}_CPMnorm_pl.bedgraph
        echo $cSample | awk '{ c="'$cSample'"; printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*1000000)/c}' ${x}_unnorm_mn.bedgraph  > ${x}_CPMnorm_mn.bedgraph
 
        echo $cSample | awk '{ c="'$cSample'"; printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*1000000)/c}' ${x}_allMap_unnorm_pl.bedgraph  > ${x}_allMap_CPMnorm_pl.bedgraph
        echo $cSample | awk '{ c="'$cSample'"; printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*1000000)/c}' ${x}_allMap_unnorm_mn.bedgraph  > ${x}_allMap_CPMnorm_mn.bedgraph
        ## generates the coverage normalized bedgraph files by dividing the fourth column of bedgraph file with total number of mapped reads and multiplying the resulting number with 1M (RPM - Reads Per Million)

        echo; echo "making bigwig from CPM-normalized bedgraph files of ${x}"
        bedGraphToBigWig ${x}_CPMnorm_pl.bedgraph ${chromSizes} ${x}_CPMnorm_pl.bigWig
        bedGraphToBigWig ${x}_CPMnorm_mn.bedgraph ${chromSizes} ${x}_CPMnorm_mn.bigWig
        ## makes bigwig file from the bedgraph files

        echo; echo "making bigwig from CPM-normalized bedgraph files of allMapping reads of ${x}"
        bedGraphToBigWig ${x}_allMap_CPMnorm_pl.bedgraph ${chromSizes} ${x}_allMap_CPMnorm_pl.bigWig
        bedGraphToBigWig ${x}_allMap_CPMnorm_mn.bedgraph ${chromSizes} ${x}_allMap_CPMnorm_mn.bigWig
        ## makes bigwig file from the bedgraph files

        rm *.temp

    done
    cd ../../


    #-----------------------------------------------------------------------

### Map functional genomic regions

        # Check if the gene list exists
    if test -e "${organism}/genome/genes.txt"; then
        echo
    else
        echo; echo "Downloading a gene list from GitHub..."

        if [ $organism == "hs1" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_hs1.txt > ${organism}/genome/genes.txt
        elif [ $organism == "hg38" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_hg38.txt > ${organism}/genome/genes.txt
        elif [ $organism == "mm39" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_mm39.txt > ${organism}/genome/genes.txt
        elif [ $organism == "dm6" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_dm6.txt > ${organism}/genome/genes.txt
        elif [ $organism == "canFam6" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_canFam6.txt > ${organism}/genome/genes.txt
        elif [ $organism == "TAIR" ]; then
            wget -q --no-check-certificate -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/genes/genes_TAIR.txt > ${organism}/genome/genes.txt
        else
            read -p "Please enter a URL to a gene list with the .gff.gz filetype: " URL_genes
            wget -q --no-check-certificate -cO genes_${organism}.gff.gz $URL_genes

            #-----------------------------------------

            gunzip genes_${organism}.gff.gz

            #-----------------------------------------

            # extract all rows containing genes:
            awk '$3 == "gene" {print $1, $4, $5, $7, $9}' genes_${organism}.gff > extracted_genes.txt

            # extract all gene names:
            sed -Ee 's/(.*Name=)([^;]*)(;gbkey.*)/\2/' -Ee 's/(.*Name=)([^;]*)(;description.*)/\2/' -Ee 's/(.*ID=gene-Dmel_)([^;]*)(;Dbxref.*)/\2/' extracted_genes.txt > modified_file.txt

            # remove all data except chr, txstart, txend and strand:
            cut -d' ' -f 1-4 extracted_genes.txt > temp.txt; paste -d' ' temp.txt modified_file.txt > genes0.temp

            #-----------------------------------------

            # extract old chromosome names as a list
            awk '{print $1}' genes0.temp > chromosomes.temp
            chromosomes_old=($(cat chromosomes.temp | tr '[:space:]' '[\n*]' | uniq))

            #-----------------------------------------

            # create list of new chromosome lists dependant on type of organism
            echo "Input the name of each chromosome for the genome ${organism}"
            chromosomes_new=()
            while IFS= read -r -p "Next item (end with an empty line): " line; do
                [[ $line ]] || break  # break if line is empty
                chromosomes_new+=("$line")
            done

            #-----------------------------------------

            # determine number of chromosomes in list
            array_length=${#chromosomes_new[@]}

            # rename chromosomes in gene-file from old to new chromosomes
            for i in $(seq 0 $((array_length-1))); do
            sed "s/${chromosomes_old[${i}]}/${chromosomes_new[${i}]}/g" genes0.temp > temporary.temp && mv temporary.temp genes0.temp
            done

            # add titles to the different columns in the file and change all seperators to tabs
            echo "chr txStart txEnd strand geneName" | cat - genes0.temp > genes0_temp.temp && mv genes0_temp.temp genes0.temp
            sed 's/ /\t/g' genes0.temp > temp.temp && mv temp.temp genes_txt.temp

            pattern=$(IFS="|"; echo "${chromosomes_new[*]}")
            grep -E "$pattern" genes_txt.temp > genes.txt
            mv genes.txt ${organism}/genome/

            # remove temporary files
            rm *.temp genes_${organism}.gff temp.txt extracted_genes.txt modified_file.txt 
        fi

        #-----------------------------------------

        echo "List of genes created.."
        echo
    fi

    #-----------------------------------------------------------------------

    #Make a list of your samples, without the .fastq extension

    cd ${organism}/samples
    declare -a sample_list
    for file in *.fastq; do
        sample_list=("${sample_list[@]}" "$file")
    done
    sample_list=("${sample_list[@]%.*}")
    cd ../../

    #---------------------------------------------------

    for sample in ${sample_list[@]}; do 

        echo "Creating functional genomics .bed-file for sample $sample"

        mkdir -p ${organism}/analysis/functionalGenomics_${sample}          # create a folder for output files

        #---------------------------------------------------

        wget -q -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/scripts/R1.R > R1.R
        wget -q -cO - https://raw.githubusercontent.com/SerhatAktay/TrackTx/master/scripts/R2.R > R2.R

        Rscript R1.R $organism $sample --save   # R_script that adds coordinates for different features such as TSS, CPS, genebody etc.

        cat ${organism}/bigWig/${sample}_allMap_unnorm_mn.bedgraph | tr ' ' '\t' > ${organism}/bigWig/${sample}_allMap_unnorm_mn_temp.bedgraph

        bedtools slop -i ${organism}/bigWig/${sample}_allMap_unnorm_pl.bedgraph -g ${organism}/genome/chrm_sizes_${organism}.txt -b 10 | bedtools merge -d 1 -c 4 -o max > ${organism}/analysis/functionalGenomics_${sample}/pl_bedgraph.bed
        bedtools slop -i ${organism}/bigWig/${sample}_allMap_unnorm_mn_temp.bedgraph -g ${organism}/genome/chrm_sizes_${organism}.txt -b 10 | bedtools merge -d 1 -c 4 -o max > ${organism}/analysis/functionalGenomics_${sample}/mn_bedgraph.bed
        bedtools intersect -a ${organism}/analysis/functionalGenomics_${sample}/pl_bedgraph.bed -b ${organism}/analysis/functionalGenomics_${sample}/mn_bedgraph.bed > ${organism}/analysis/functionalGenomics_${sample}/divergent_transcription.bed

        #-----------------------------------------------

        bedtools intersect -v -a ${organism}/analysis/functionalGenomics_${sample}/divergent_transcription.bed -b ${organism}/analysis/functionalGenomics_${sample}/refGenes_TSSpm500.txt > ${organism}/analysis/functionalGenomics_${sample}/enhancers.bed
        bedtools intersect -wa -a ${organism}/analysis/functionalGenomics_${sample}/refGenes_TSSpm500.txt -b ${organism}/analysis/functionalGenomics_${sample}/divergent_transcription.bed > ${organism}/analysis/functionalGenomics_${sample}/activeGenes.bed
        Rscript R2.R $organism $sample --save

        rm R1.R R2.R

        #-----------------------------------------------

        echo; echo "retaining 3' most coordinate of the bed file."

        cd ${organism}/analysis/functionalGenomics_${sample}

        project_root="$(dirname $(dirname $(dirname $(realpath $0))) )"
        awk '$6 == "+"' ${project_root}/bigWig/${sample}.bed > tempPL.bed
        awk '$6 == "-"' ${project_root}/bigWig/${sample}.bed > tempMN.bed
        awk '{$2 = $3; print}' tempPL.bed > tempPL_3p.bed
        awk '{$3 = $2; print}' tempMN.bed > tempMN_3p.bed

        cat tempPL_3p.bed tempMN_3p.bed | tr ' '  '\t' > PROseq_3pnt.bed

        rm *temp*

        #-----------------------------------------------

        echo; echo "creating .bed files"

        ##counting engaged Pol II at promoter-proximal regions
        bedtools intersect -s -u -wa -a PROseq_3pnt.bed -b ppPolII.txt > PROseq_ppPolII.bed
        bedtools intersect -v -a     PROseq_3pnt.bed -b ppPolII.txt > ppRemoved.bed

        ##counting engaged Pol II at the sites of divergent transcription
        bedtools intersect -S -u -wa -a ppRemoved.bed -b divTx.txt > PROseq_ppDiv.bed
        bedtools intersect -v -a     ppRemoved.bed -b divTx.txt > ppdivRemoved.bed

        ##counting engaged Pol II at enhancers
        bedtools intersect -u -wa -a ppdivRemoved.bed -b enhancers.bed > PROseq_enhancers.bed
        bedtools intersect -v -a     ppdivRemoved.bed -b enhancers.bed > ppdivEnhRemoved.bed

        ##counting engaged Pol II at CPS
        bedtools intersect -s -u -wa -a ppdivEnhRemoved.bed -b CPS.txt > PROseq_CPS.bed
        bedtools intersect -v -a     ppdivEnhRemoved.bed -b CPS.txt > ppdivEnhCPSRemoved.bed

        ##counting engaged Pol II at GB
        bedtools intersect -s -u -wa -a ppdivEnhCPSRemoved.bed -b geneBody.txt > PROseq_GB.bed
        bedtools intersect -v -a     ppdivEnhCPSRemoved.bed -b geneBody.txt > ppdivEnhCPSgbRemoved.bed

        ##counting engaged Pol II at short genes
        bedtools intersect -s -u -wa -a ppdivEnhCPSgbRemoved.bed -b shortGenes.txt > PROseq_SG.bed
        bedtools intersect -v -a     ppdivEnhCPSgbRemoved.bed -b shortGenes.txt > ppdivEnhCPSgbSGRemoved.bed

        ##counting engaged Pol II at termination windows
        bedtools intersect -s -u -wa -a ppdivEnhCPSgbSGRemoved.bed -b TW.txt > PROseq_TW.bed
        bedtools intersect -v -a     ppdivEnhCPSgbRemoved.bed -b TW.txt > PROseq_noGene_noEnh.bed

        rm *Removed.bed

        echo; echo ".bed files created"

        #-----------------------------------------------

        echo; echo "here comes count data:"

        ## Total count of active sites of transcription (uniquely mapping reads) in the dataset:
        echo total reads:
        totC=$(awk 'END { print NR }' PROseq_3pnt.bed)
        echo ${totC}
        ## Engaged Pol II at transcribed enhancers:
        echo enhancers:
        enhC=$(awk 'END { print NR }' PROseq_enhancers.bed)
        echo ${enhC}

        ## Engaged Pol II at divergent transcripts:
        echo divergentTx:
        divC=$(awk 'END { print NR }' PROseq_ppDiv.bed)
        echo ${divC}

        ## Engaged Pol II at promoters:
        echo promoter:
        ppC=$(awk 'END { print NR }' PROseq_ppPolII.bed)
        echo ${ppC}

        ## Engaged Pol II at gene bodies:
        echo gene body:
        GBC=$(awk 'END { print NR }' PROseq_GB.bed)
        echo ${GBC}

        ## Engaged Pol II at the cleavage and polyadenylation regions:
        echo CPS:
        CPSC=$(awk 'END { print NR }' PROseq_CPS.bed)
        echo ${CPSC}

        ## Engaged Pol II at short genes:
        echo short genes:
        SGC=$(awk 'END { print NR }' PROseq_SG.bed)
        echo ${SGC}

        ## Engaged Pol II at the termination windows:
        echo termination window:
        twC=$(awk 'END { print NR }' PROseq_TW.bed)
        echo ${twC}

        ## Engaged Pol II at sites that did not localize to any of the identified functional regions:
        echo noGenesNoEnhancers:
        noGdnoE=$(awk 'END { print NR }' PROseq_noGene_noEnh.bed)
        echo ${noGdnoE}

        #-----------------------------------------------

        echo; echo "creating a bed-file named functionalGenomicRegions for further analysis"

        bedtools intersect -c -wa -a ppPolII.txt -b PROseq_ppPolII.bed > ppPolCounts.tmp
        bedtools intersect -c -wa -a divTx.txt -b PROseq_ppDiv.bed > ppDivCounts.tmp
        bedtools intersect -c -wa -a enhancers.bed -b PROseq_enhancers.bed > enhancerCounts.tmp
        bedtools intersect -c -wa -a geneBody.txt -b PROseq_GB.bed > geneBodyCounts.tmp
        bedtools intersect -c -wa -a CPS.txt -b PROseq_CPS.bed > CPSCounts.tmp
        bedtools intersect -c -wa -a shortGenes.txt -b PROseq_SG.bed > SGCounts.tmp
        bedtools intersect -c -wa -a TW.txt -b PROseq_TW.bed > TerminationWinCounts.tmp

        cut -f 1-4,6- -d " " enhancerCounts.tmp > enhancerCounts.tmp

        #-----------------------------------------------

        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="243,132,0"; print }' ppPolCounts.tmp > ppPolCounts.bed                         # orange
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="178,59,212"; print }' ppDivCounts.tmp > ppDivCounts.bed                        # purple
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="115,212,122"; print }' enhancerCounts.tmp > enhancerCounts.bed                 # green
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="0,0,0"; print }' geneBodyCounts.tmp > geneBodyCounts.bed                       # black
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="103,200,249"; print }' CPSCounts.tmp > CPSCounts.bed                           # blue
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="253,218,13"; print }' SGCounts.tmp > SGCounts.bed                              # yellow
        awk -F '\t' -v OFS='\t' '{ $(NF+1) ="255,54,98"; print }' TerminationWinCounts.tmp > TerminationWinCounts.bed       # red

        #-----------------------------------------------

        cat ppPolCounts.bed ppDivCounts.bed enhancerCounts.bed geneBodyCounts.bed CPSCounts.bed SGCounts.bed TerminationWinCounts.bed > catRegions.temp

        awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $7 "\t" $8}' catRegions.temp > catRegions2.temp
        awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $4 "\t" $5 "\t" $2 "\t" $3 "\t" $7}' catRegions2.temp > catRegions3.temp
        awk '!seen[$1,$2,$3]++' catRegions3.temp | sortBed > catRegions4.temp
        #awk -F '\t' -v OFS='\t' '{ $(NF+1) ="."; print }' catRegions3.temp > catRegions4.temp

        touch headerLine.txt
        echo track name="functional_genomic_regions_$sample" itemRgb="On" >> headerLine.txt

        cat headerLine.txt catRegions4.temp > functionalGenomicRegions_${sample}.bed
        awk -v c=${totC} '{ $4 /= c; print }' OFS='\t' functionalGenomicRegions_${sample}.bed > functionalGenomicRegions_${sample}_norm.bed
        awk -i inplace -v OFS='\t' '{ $4 *= 1000000; print }' functionalGenomicRegions_${sample}_norm.bed

        mv functionalGenomicRegions_${sample}.bed ${project_root}/analysis/functionalGenomicRegions_${sample}.bed
        mv functionalGenomicRegions_${sample}_norm.bed ${project_root}/analysis/functionalGenomicRegions_${sample}_norm.bed

        #-----------------------------------------------

        rm *.temp
        rm *.tmp
        rm headerLine.txt
        cd ../../../

        #-----------------------------------------------

    done

    echo "Analysis of functional genomic regions done"
    echo

    #---------------------------------------------------

    cd ../../
}

#-----------------------------------------------------------------------


while true; do
    options=("Download genome and build bowtie2 alignment" "Download/load experimental data" "Perform quality control" 
         "Convert fastq files" "Map functional genomic regions" "Do step 1-5 all in one"
         "Compare gene expression using DEseq2" "Map Sat III repeats" "Quit")
    echo "Choose an option: "
    select opt in "${options[@]}"; do
        case $REPLY in
            1) genome; break ;;
            2) sample; break ;;
            3) check_quality; break ;;
            4) convert; break ;; 
            5) functionalgenomicregions; break ;;
            6) download_and_convert; break;;
            7) CompGeneExpression; break ;;
            8) map_repeats; break ;;
            9) break 2 ;;
            *) echo "Invalid option. Try another one." >&2
        esac
    done
done

echo "Bye bye!"

#-----------------------------------------------------------------------

# Deactivate Conda environment
if [ "$conda_environment_active" = true ]; then
    conda deactivate
fi

#-----------------------------------------------------------------------
