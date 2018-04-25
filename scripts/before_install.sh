#!/bin/bash

set -e

if [[ -z "$encrypted_9f3357583363_key" || -z "$encrypted_9f3357583363_iv" ]]; then
    exit 0
fi

openssl aes-256-cbc -K $encrypted_9f3357583363_key -iv $encrypted_9f3357583363_iv -in ./scripts/id_rsa.enc -out ~/.ssh/id_rsa -d
chmod 600 ~/.ssh/id_rsa

echo "StrictHostKeyChecking no" >> ~/.ssh/config

git config --global user.email "philopon.dependence@gmail.com"
git config --global user.name "philopon"
