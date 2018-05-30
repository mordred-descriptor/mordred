#!/bin/bash

if [[ -n "$APPVEYOR" ]]; then
    export PATH="$MINICONDA:$MINICONDA/Scripts:$PATH"
else
    export PATH="$HOME/miniconda/bin:$PATH"
fi

info() {
    echo ">>> $@"
    "$@"
}

banner() {
    line="=================="
    echo "$line $@ $line"
}
