#!/bin/bash

TAG=commandlinegirl/sashimi:hom_range_test
#git commit --amend
#git push -f
docker build -t $TAG docker/ --build-arg CACHEBUST=$(date +%s)
docker push $TAG
