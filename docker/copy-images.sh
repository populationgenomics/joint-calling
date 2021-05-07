#!/bin/bash

set -ex

gcloud config set project cpg-common
images=$(cat images.txt)
dest_repo="australia-southeast1-docker.pkg.dev/cpg-common/joint-calling"

for source in ${images}
do
  image_name=$(basename "${source}")
  if [ -e Dockerfile ] ; then trap "Remove Dockefile from this folder" ; fi
  echo "FROM ${source}" > Dockerfile
  dest="${dest_repo}/${image_name}"
  gcloud builds submit --tag "${dest}" .
  if [ -e Dockerfile ] ; then rm Dockerfile ; fi
done

set +ex
