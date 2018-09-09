#!/bin/bash

VERSION=0.3
TAG=commandlinegirl/sashimi:$VERSION
docker build -t $TAG docker/ --build-arg CACHEBUST=$(date +%s)
docker tag commandlinegirl/sashimi:$VERSION commandlinegirl/sashimi:latest
docker push $TAG
