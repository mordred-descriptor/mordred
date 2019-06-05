#!/bin/bash

set -e

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo this is PR
    exit 0
fi

mkdir -p ~/.ssh
chmod 700 ~/.ssh

openssl aes-256-cbc -K $encrypted_9f3357583363_key -iv $encrypted_9f3357583363_iv -in ./extra/id_rsa.enc -out ~/.ssh/id_rsa -d
chmod 600 ~/.ssh/id_rsa

echo "StrictHostKeyChecking no" >> ~/.ssh/config

git config --global user.email "philopon.dependence@gmail.com"
git config --global user.name "philopon"
