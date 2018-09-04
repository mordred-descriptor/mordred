#!/bin/bash
set -e

source ./extra/ci/common.sh

# install conda
if [[ -n "$TRAVIS_OS_NAME" ]]; then
    if [[ "$TRAVIS_OS_NAME" == osx ]]; then
        export OS_NAME=MacOSX
    elif [[ "$TRAVIS_OS_NAME" == linux ]]; then
        export OS_NAME=Linux
    fi
    if [[ ! -f "miniconda.sh" ]]; then
        info wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-${OS_NAME}-x86_64.sh -O miniconda.sh
    fi
    rm -rf $HOME/miniconda
    info bash miniconda.sh -b -p $HOME/miniconda
elif [[ -n "$APPVEYOR" ]]; then
    export OS_NAME=Windows
fi

# setup conda
hash -r

info conda config --set always_yes yes --set changeps1 no
info conda config --add channels rdkit --add channels mordred-descriptor
info conda update -y --all

info conda install python=$PYTHON_VERSION

RDKIT="rdkit==$(python ./extra/ci/get-rdkit-version.py ./extra/requirements/rdkit-versions.txt $OS_NAME $PYTHON_VERSION)"
info conda install $RDKIT --file ./extra/requirements/requirements-conda.txt

info pip install -r ./extra/requirements/requirements-pip.txt
