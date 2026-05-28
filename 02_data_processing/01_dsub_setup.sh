# This file is meant to be sourced in the terminal or jupyter notebook. It sets up dsub and defines 
# a function for running dsub with reasonable defaults for AoU RWB. Dsub is a command-line tool 
# for running batch jobs on Google Cloud. It provides a simple interface for running jobs in paralle,
# selectting appropriate VMs, and managing resources such as CPUs, memory, and storage.

# Setup dsub
dsub --version

# optional
pip3 install --upgrade dsub

# Verify dsub installation and configuration options
dsub --help

# Setup aou_dsub function
%%writefile ~/aou_dsub.bash

#!/bin/bash

# --[ Parameters ]--
# any valid dsub parameter flag

#--[ Returns ]--
# the job id of the job created by dsub

#--[ Details ]--
# The first seven parameters below should always be those values when running on AoU RWB.
# Feel free to change the values for --user, --regions, --logging, and --image if you like.

# Note that we insert some job data into the logging path.
# https://github.com/DataBiosphere/dsub/blob/main/docs/logging.md#inserting-job-data

function aou_dsub () {

  # Get a shorter username to leave more characters for the job name.
  local DSUB_USER_NAME="$(echo "${OWNER_EMAIL}" | cut -d@ -f1)"
  
  dsub \
      --provider google-batch \
      --user-project "${GOOGLE_PROJECT}" \
      --project "${GOOGLE_PROJECT}" \
      --network "global/networks/network" \
      --subnetwork "regions/us-central1/subnetworks/subnetwork" \
      --service-account "$(gcloud config get-value account)" \
      --use-private-address \
      --user "${DSUB_USER_NAME}" \
      --image "${ARTIFACT_REGISTRY_DOCKER_REPO}/ubuntu:latest" \
      --regions us-central1 \
      --logging "${WORKSPACE_BUCKET}/dsub/logs/{job-name}/{user-id}/$(date +'%Y%m%d/%H%M%S')/{job-id}-{task-id}-{task-attempt}.log" \
      "$@"
}
