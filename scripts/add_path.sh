#!/bin/bash

if [[ -n "$APPVEYOR" ]]; then
    export PATH="$MINICONDA:$MINICONDA/Scripts:$PATH"
else
    export PATH="$HOME/miniconda/bin:$PATH"
fi
