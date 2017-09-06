#!/bin/bash

set -e
source ./scripts/add_path.sh

if [[ -n "$COVERAGE" ]]; then
    python -m mordred.tests -q --with-coverage
else
    python -m mordred.tests -q
fi
