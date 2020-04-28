#!/bin/bash

# Handle minio stored optionally centrally. Use if on PATH.
if command -v minio > /dev/null 2>&1; then
  MINIO=minio
else
  MINIO=./minio
fi

MINIO_ACCESS_KEY=hello MINIO_SECRET_KEY=k1tty-pass $MINIO server ./stash
