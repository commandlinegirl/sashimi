#!/usr/bin/env bash

dx upload wdl/sashimi_simple.wdl
dx run workflow_importer -iprimary_descriptor=sashimi.wdl
