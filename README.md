# TrackTx - The Integrative Data Processing Pipeline for Tracking Transcription

## Instructions
The pipeline is designed to run in a terminal window. Follow the instructions below. 

>**Terminal commands to run the pipeline:**
>1. Clone: git clone https://github.com/SerhatAktay/TrackTx.git
>2. Change directory: cd TrackTx
>3. Start pipeline: ./TrackTx.sh

### Prerequisites

#### Software
This pipeline uses bash and R and the required packages will all be installed automatically using Conda and it's recommended that you install [miniConda](https://docs.conda.io/projects/miniconda/en/latest/) before running the pipeline. 

#### Data
The pipeline downloads all data needed, including reference genomes and experimental data (from GEO). If you want to use your own reference genome and/or your own locally stored data you can do so by following the instructions in the pipeline. 

### Recommended usage
It's recommended to run the pipeline using a virtual machine or similar that will lessen the burden on your computer. It's also recommended to run the pipeline using a window manager such as [screen](https://linux.die.net/man/1/screen) so that you don't need to keep the computer on and connected at all times. It takes up to 10 hours to run two human samples so allocate enough time. 

## Contact
If you run into errors, bugs, or have questions either use the "Issues" section on this repo or contact me at <serhat.aktay@scilifelab.se>. 
