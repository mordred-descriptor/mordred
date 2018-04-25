#!/bin/bash

set -e

if [[ ! -f ~/.ssh/id_rsa ]]; then
    exit 0
fi

source ./scripts/add_path.sh

[[ -n "$COVERAGE" ]] && coveralls

if [[ -z "$TRAVIS_TAG" && -z "$APPVEYOR_REPO_TAG_NAME" ]]; then
    LABEL=dev
    echo $(cat mordred/_version.txt).post1.dev1 > mordred/_version.txt
else
    LABEL=main
fi

conda build . --no-test

OUTPUT=`conda build . --output --python $PYTHON_VERSION`
if [[ -n "$APPVEYOR" ]]; then
    cmd /c "anaconda -t $ANACONDA_CLOUD_TOKEN upload --label $LABEL --force $OUTPUT"
else
    anaconda -t $ANACONDA_CLOUD_TOKEN upload --label $LABEL --force $OUTPUT
fi

# documentation
if [[ "$TRAVIS_PULL_REQUEST" == false && -n "$DOCUMENTATION" && "$TRAVIS_OS_NAME" == linux ]]; then
    eval $(ssh-agent -s)
    ssh-add
    ssh-add -l
    echo "$SSH_AGENT_PID"

    cd docs
    make html

    git clone -b gh-pages $DOC_REMOTE gh-pages
    rm -r gh-pages/$TRAVIS_BRANCH || true
    cp -r _build/html gh-pages/$TRAVIS_BRANCH

    cd gh-pages
    git add .
    git commit -m "update documentation to mordred-descriptor/mordred@$TRAVIS_COMMIT"
    git push origin gh-pages
fi
