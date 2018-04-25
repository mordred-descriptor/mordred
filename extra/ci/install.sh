#!/bin/bash

set -e

REQUIREMENTS=./extra/requirements
CONDA_REQ_FILE=requirements-conda.txt
PIP_REQ_FILE=requirements-pip.txt

# install conda

if [[ "$TRAVIS_OS_NAME" == osx ]]; then
    export OS_NAME=MacOSX
elif [[ "$TRAVIS_OS_NAME" == linux ]]; then
    export OS_NAME=Linux
elif [[ -n "$APPVEYOR" ]]; then
    export OS_NAME=Windows
fi

if [[ -n "$TRAVIS_OS_NAME" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-${OS_NAME}-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
fi

source ./extra/ci/conda.sh

# setup conda

hash -r

conda config --set always_yes yes --set changeps1 no
conda config --add channels rdkit --add channels mordred-descriptor

# install requirements

cat $REQUIREMENTS/requirements-conda.txt $REQUIREMENTS/requirements.txt > $CONDA_REQ_FILE

if [[ "$PYTHON_VERSION" < 3.4 ]]; then
    echo enum34 >> $CONDA_REQ_FILE
fi

if [[ -n "$COVERAGE" ]]; then
    echo coveralls >> $PIP_REQ_FILE
fi

if [[ -n "$LINT" ]]; then
    cat $REQUIREMENTS/requirements-flake8.txt >> $PIP_REQ_FILE
fi

if [[ -n "$DOCUMENTATION" ]]; then
    cat $REQUIREMENTS/requirements-documentation.txt >> $PIP_REQ_FILE
fi

conda install python=$PYTHON_VERSION --file $CONDA_REQ_FILE

if [[ -f $PIP_REQ_FILE ]]; then
    pip install -r $PIP_REQ_FILE
fi
