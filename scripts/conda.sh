#!/bin/bash

echo "add conda PATH"
if [[ -n "$APPVEYOR" ]]; then
    export PATH="$MINICONDA:$MINICONDA/Scripts:$PATH"
else
    export PATH="$HOME/miniconda/bin:$PATH"
fi

hash -r

conda config --set always_yes yes --set changeps1 no
conda config --add channels rdkit --add channels mordred-descriptor
