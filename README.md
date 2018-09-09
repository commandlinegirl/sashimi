# Sashimi

## Description

...

## Requirements

Docker is required to run sashimi.

If you want to run standalone scripts instead, install:
* numpy
* scipy
* pandas
* pybedtools

for example with pip:
```
pip install \
    numpy \
    scipy \
    pandas \
    pybedtools
```

## Running sashimi

### In a docker container

```
docker run commandlinegirl/sashimi:0.3 /opt/sashimi.py -h
```

### Using the sashimi script directly

```
python sashimi.py -h
```
