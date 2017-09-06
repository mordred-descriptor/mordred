#!/bin/bash

set -e
source ./scripts/add_path.sh

python -m mordred.tests -q --with-coverage
