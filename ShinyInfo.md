# Info on Python Shiny and Shiny.io

#### - look into pipes in R
#### - look into dockers and containers 
####  https://code.visualstudio.com/docs/containers/quickstart-python

## If you want to download data from SRA by NCBI run following in terminal if computer is Mac arm64
### conda create -n ngstools
### conda activate ngstools
### conda config --env --set subdir osx-64
### config --env --add channels defaults
### conda config --env --add channels bioconda
### conda config --env --add channels conda-forge
### conda install sra-tools
### vdb-config  -i 
#### You want to enable the "Remote Access" option on the Main screen.Proceed to the "Cache" tab where you will want to enable "local file-caching" and you want to set the "Location of user-repository".The repository directory needs to be set to an empty folder. This is the folder where prefetch will deposit the files. Go to your cloud provider tab and accept to "report cloud instance identity".
## Do this every time you want to try to download data from here 

