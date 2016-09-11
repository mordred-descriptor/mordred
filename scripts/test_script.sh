#!/bin/bash

set -e
source ./scripts/add_path.sh

nosetests mordred -q
