#!/bin/sh
set -eu
cd "$(dirname "$0")"
mkdir -p .buildutils/m4
autoreconf -fvi
