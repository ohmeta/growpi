#!/bin/sh
rm -rf build
rm -rf growpi.egg-info
rm -rf dist
python3 setup.py sdist bdist_wheel
twine upload dist/*
