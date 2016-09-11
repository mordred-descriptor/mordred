#!/bin/bash

set -e

eval "$(ssh-agent -s)"
ssh-add -D

openssl aes-256-cbc -K $encrypted_8d665ddaa7bd_key -iv $encrypted_8d665ddaa7bd_iv -in .id_rsa.enc -out ~/.ssh/.id_rsa -d
chmod 600 ~/.ssh/.id_rsa
ssh-add ~/.ssh/.id_rsa


echo "StrictHostKeyChecking no" >> ~/.ssh/config

git config --global user.email "philopon.dependence@gmail.com"
git config --global user.name "philopon"
