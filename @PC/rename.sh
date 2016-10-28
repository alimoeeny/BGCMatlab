#!/bin/sh
for f in xx*; do
  #d="$(head -1 "$f" | sed 's/function .*= \(.*\)(.*$/\1.m/')"
  d="$(head -1 "$f" | sed 's/function \(.*\)(.*$/\1.m/')"
  if [ ! -f "$d" ]; then
    mv $f $d
    echo "moving '$f' to '$d'"
  else
    echo "File '$d' already exists! skipping '$f'"
  fi
done
