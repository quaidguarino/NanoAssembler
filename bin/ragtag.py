#!/bin/sh
DIR=$(cd "$(dirname "$0")" && pwd)
exec "$DIR/../python/bin/python3" "$DIR/../python/bin/ragtag.py" "$@"
