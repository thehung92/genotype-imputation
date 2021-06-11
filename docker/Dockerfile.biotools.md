# load container from biotools and add required file
FROM thehung92phuyen/biotools:v3.5
LABEL MAINTAINER="thehung, thehung92phuyen@gmail.com"
# update general tools in existing conda environment
SHELL ["conda", "run", "-n", "biotools", "/bin/bash", "-c"]

```bash
# copy minimac3 and minimac4 to container
# docker run --rm -v binary-files:/mnt/ -it thehung92phuyen/biotools:v3.5 bash
cp -r /mnt/* /opt/
ln -s /opt/Minimac3/bin/Minimac3 /usr/local/bin/minimac3
ln -s /opt/Minimac4/bin/minimac4 /usr/local/bin/minimac4
#
apt-get install -y tree
#
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y shapeit4
# detach container
# docker commit -m 'add shapeit4, minimac3, minimac4' 63fa4c9d529a thehung92phuyen/biotools:v4.0
```