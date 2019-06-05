#!/bin/bash

set -e
source ./extra/ci/common.sh

PYTHON=python

conda list

if [[ -n "$COVERAGE" ]]; then
    info $PYTHON -m mordred.tests -q --with-coverage
else
    info $PYTHON -m mordred.tests -q
fi

echo "test README.rst" >&2
info $PYTHON -m doctest README.rst

for example in `find examples -name '*.py'`; do
    echo "test $example" >&2
    PYTHONPATH=. info $PYTHON $example > /dev/null
done

info $PYTHON setup.py flake8
