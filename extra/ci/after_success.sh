#!/bin/bash

set -e

source ./extra/ci/common.sh

if [[ -n "$COVERAGE" ]]; then
    info coveralls
fi

if [[ -z "$TRAVIS_TAG" && -z "$APPVEYOR_REPO_TAG_NAME" ]]; then
    LABEL=dev
    echo $(cat mordred/_version.txt).post1.dev1 > mordred/_version.txt
else
    LABEL=main
fi

info conda build . --no-test

OUTPUT=`conda build . --output --python $PYTHON_VERSION`
if [[ -n "$ANACONDA_CLOUD_TOKEN" ]]; then
    if [[ -n "$APPVEYOR" ]]; then
        cmd /c "anaconda -t $ANACONDA_CLOUD_TOKEN upload --label $LABEL --force $OUTPUT"
    else
        anaconda -t $ANACONDA_CLOUD_TOKEN upload --label $LABEL --force $OUTPUT
    fi
fi

# documentation
if [[ -f ~/.ssh/id_rsa && "$TRAVIS_PULL_REQUEST" == false && -n "$DOCUMENTATION" && "$TRAVIS_OS_NAME" == linux ]]; then
    eval $(ssh-agent -s)
    ssh-add
    ssh-add -l
    echo "$SSH_AGENT_PID"

    cd docs
    info make html

    info git clone -b gh-pages $DOC_REMOTE gh-pages
    info rm -r gh-pages/$TRAVIS_BRANCH || true
    info cp -r _build/html gh-pages/$TRAVIS_BRANCH

    cd gh-pages
    info git add .
    info git commit -m "update documentation to mordred-descriptor/mordred@$TRAVIS_COMMIT"
    info git push origin gh-pages
fi
