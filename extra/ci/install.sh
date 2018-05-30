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
    info wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-${OS_NAME}-x86_64.sh -O miniconda.sh
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

# install requirements
info pip install pipenv

pipenv lock -r > requirements.txt
info python ./extra/ci/scrub-requirements.py requirements.txt

banner "requirements.txt start"
cat requirements.txt
banner "requirements.txt  end "

RDKIT="rdkit==$(python ./extra/ci/get-rdkit-version.py $OS_NAME $PYTHON_VERSION)"
info conda install $RDKIT --file requirements.txt --file ./extra/ci/requirements-conda.txt

pipenv lock -r --dev > requirements.txt
banner "requirements.txt(dev) start"
cat requirements.txt
banner "requirements.txt(dev)  end "

info pip install -r requirements.txt
