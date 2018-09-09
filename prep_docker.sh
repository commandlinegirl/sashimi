#!/bin/bash

TAG=commandlinegirl/sashimi:0.3
#git commit --amend
#git push -f
docker build -t $TAG docker/ --build-arg CACHEBUST=$(date +%s)
docker push $TAG
