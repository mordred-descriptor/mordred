#!/bin/bash

set -e

# install conda

if [[ "$TRAVIS_OS_NAME" == osx ]]; then
    export OS_NAME=MacOSX
elif [[ "$TRAVIS_OS_NAME" == linux ]]; then
    export OS_NAME=Linux
elif [[ -n "$APPVEYOR" ]]; then
    export OS_NAME=Windows
fi
echo "detected OS: $OS_NAME"

if [[ -n "$TRAVIS_OS_NAME" ]]; then
    echo "download & install conda"
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-${OS_NAME}-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
fi

source ./scripts/conda.sh

# install requirements

CONDA_REQ_FILE=requirements-conda.txt
PIP_REQ_FILE=requirements-pip.txt

cat ./scripts/requirements-conda.txt ./scripts/requirements.txt > $CONDA_REQ_FILE

[[ "$PYTHON_VERSION" < 3.4 ]] && echo enum34 >> $CONDA_REQ_FILE

[[ -n "$COVERAGE" ]] && echo coveralls >> $PIP_REQ_FILE
[[ -n "$LINT" ]] && cat ./scripts/requirements-flake8.txt >> $PIP_REQ_FILE

if [[ -n "$DOCUMENTATION" ]]; then
    echo sphinx >> $CONDA_REQ_FILE
    echo sphinx_rtd_theme >> $CONDA_REQ_FILE
    echo sphinxcontrib-bibtex >> $PIP_REQ_FILE
fi

echo "=== conda install ==="
cat $CONDA_REQ_FILE
echo "=== pip   install ==="
cat $PIP_REQ_FILE
echo "====================="

conda install python=$PYTHON_VERSION --file $CONDA_REQ_FILE

if [[ -f $PIP_REQ_FILE ]]; then
    pip install -r $PIP_REQ_FILE
fi
