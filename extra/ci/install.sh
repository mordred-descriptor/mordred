#!/bin/bash
set -e

source ./extra/ci/common.sh

EXTENSION=sh
if [[ "$OS_NAME" == Windows ]]; then
    EXTENSION=exe
fi

MINICONDA_INSTALLER=Miniconda3-latest-${OS_NAME}-x86_64.$EXTENSION

if [[ ! -f "$MINICONDA_INSTALLER" ]]; then
    info wget -oq https://repo.continuum.io/miniconda/$MINICONDA_INSTALLER
else
    rm -rf $HOME/miniconda
fi

if [[ "$OS_NAME" == Windows ]]; then
    cmd.exe /C "$MINICONDA_INSTALLER /InstallationType=JustMe /RegisterPython=0 /S /D=%UserProfile%\\miniconda"
else
    info bash $MINICONDA_INSTALLER -b -p $HOME/miniconda
fi

# setup conda
hash -r

info conda config --set always_yes yes --set changeps1 no
info conda config --add channels rdkit
info conda update -y --all

info conda install python=$PYTHON_VERSION

RDKIT="rdkit==$(python ./extra/ci/get-rdkit-version.py ./extra/requirements/rdkit-versions.txt $OS_NAME $PYTHON_VERSION)"
info conda install $RDKIT --file ./extra/requirements/requirements-conda.txt

info pip install -r ./extra/requirements/requirements-pip.txt
