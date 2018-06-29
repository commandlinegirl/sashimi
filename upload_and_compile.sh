#!/usr/bin/env bash

dx upload wdl/sashimi.wdl
dx run workflow_importer -iprimary_descriptor=sashimi.wdl
