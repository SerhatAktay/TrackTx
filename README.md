# TrackTx - The Integrative Data Processing Pipeline for Tracking Transcription

## Instructions
To run the script from a terminal window download the file **TrackTx.sh** and run in terminal by using the command ./TrackTx.sh after you are in the correct folder **or**
download and put the **script.command** in a location where you want to run the script from. Double-click the file to start the pipeline and follow the instructions in the terminal window. 

### Prerequisites

#### Software
This pipeline uses bash and R and the required packages will all be installed automatically using Conda. If you do not have conda installed on your system you can do so by visiting "https://docs.conda.io/projects/miniconda/en/latest/" and installing miniconda. 

#### Data
The pipeline downloads all data needed, including reference genomes and experimental data (from GEO). If you want to use your own reference genome and/or your own locally stored data you can do so by following the instructions in the pipeline. 

### Recommended usage
It's recommended to run the pipeline using a virtual machine or similar that will lessen the burden on your computer. It's also recommended to run the pipeline in a "screen" so that you don't need to keep the computer on and connected at all times. It takes up to 10 hours to run two human samples so allocate enough time. 

