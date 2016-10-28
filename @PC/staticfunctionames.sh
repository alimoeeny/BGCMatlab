#!/bin/sh
for f in *; do
  #d="$(head -1 "$f" | sed 's/function .*= \(.*\)(.*$/\1.m/')"
  #d="$(head -1 "$f" | sed 's/function \(.*\)(.*$/\1.m/')"
  d="$(head -1 "$f")"
    echo "$d"
done
