#

## docker command to build images from dockerfile without context dir
```bash
# cd [path/contain/Dockerfile]
docker build -t test:alpha - < Dockerfile.build
```

## docker command to extract content of images
```bash
docker container create --name extract test:alpha
docker container cp extract:/opt/ ./
docker container rm -f extract
```
## docker command to test build binary on debian
```bash
docker run --rm -v ${PWD}/opt/:/opt/ -it debian:9 bash
```
## docker command to build new images from existing one
```bash
docker build -t test:beta - < Dockerfile.build2
docker run --rm -it test:beta bash
```
## docker command to build a specific stage in dockerfile
```bash
docker build --target builder -t test:alpha - < Dockerfile.build
docker build --target builder2 -t test:beta - < Dockerfile.build
docker build --target final -t test:gamma - < Dockerfile.build
docker run --rm -it test:gamma bash
```
## docker commit a succesfull build from container
```bash
docker commit -m 'binary for minimac4' build2 test:beta
```
## create volume to save binary file
```bash
docker volume create binary-files
docker run --rm -v binary-files:/mnt/ -it test:gamma bash
# cp -r /opt/Minimac* /mnt/
# exit the container then inspect volume with busybox
docker run --rm -v binary-files:/mnt/ -it busybox
```