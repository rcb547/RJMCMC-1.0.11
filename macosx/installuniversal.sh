#!/bin/bash

export PREFIX=/usr/local

pushd install
cp -vr * $PREFIX
popd

