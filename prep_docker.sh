#!/bin/bash

TAG=commandlinegirl/sashimi:0.1
#git commit --amend
#git push -f
docker build -t $TAG docker/ --build-arg CACHEBUST=$(date +%s)
docker push $TAG
