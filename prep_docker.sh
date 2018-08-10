#!/bin/bash

TAG=commandlinegirl/sashimi:0.2
#git commit --amend
#git push -f
docker build -t $TAG docker/ --build-arg CACHEBUST=$(date +%s)
docker push $TAG
