#!/bin/bash

set -e

dpkg-buildpackage -rfakeroot
