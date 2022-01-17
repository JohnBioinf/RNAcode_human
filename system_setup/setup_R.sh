#!/bin/bash

set -e

R_PATH="$(python3 -c "import sys, json;\
	print(json.load(sys.stdin)['R_PATH'])" \
	< "./parameters.json")"
if [ -z "$R_PATH" ]; then
	echo "set variables"
	exit 1
fi

if [ ! -d "$R_PATH" ]; then
	mkdir "$R_PATH"
fi

if [ ! -d "$R_PATH/libs" ]; then
	mkdir "$R_PATH/libs"
fi

export R_LIBS_USER="$R_PATH/libs"
export PATH="$R_PATH/R/bin/:$PATH"

# script_path="$(pwd)"
# cd "$R_PATH"
# wget "http://cran.rstudio.com/src/base/R-4/R-4.1.1.tar.gz"
# tar xvf "R-4.1.1.tar.gz"
# cd R-4.1.1
# ./configure --prefix="$R_PATH/R" --with-readline=no --with-x=no
# make && make install
# cd "$script_path"

Rscript ./install_packages.R
