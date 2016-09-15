#!/bin/bash

set -e
source ./scripts/add_path.sh

[[ -n "$COVERAGE" ]] && coveralls

if [[ -z "$TRAVIS_TAG" && -z "$APPVEYOR_REPO_TAG_NAME" ]]; then
    echo $(cat mordred/_version.txt).post1.dev1 > mordred/_version.txt
fi

conda build . --no-test

OUTPUT=`conda build . --output --python $PYTHON_VERSION`
if [[ -n "$APPVEYOR" ]]; then
    cmd /c "anaconda -t $ANACONDA_CLOUD_TOKEN upload --label main --force $OUTPUT"
else
    anaconda -t $ANACONDA_CLOUD_TOKEN upload --label main --force $OUTPUT
fi

# documentation
if [[ "$TRAVIS_PULL_REQUEST" == "false" && -n "$DOCUMENTATION" && ( -n "$TRAVIS_TAG" || "$TRAVIS_BRANCH" == master) ]]; then
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
