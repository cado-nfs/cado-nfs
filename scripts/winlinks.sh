#! /usr/bin/env bash

# This script looks for symbolic links stored in the Git repository,
# tests for each symlink whether it is actually a file containing the 
# name of the linked-to file (as Git produces on Windows when checking 
# out a symlink) and, if so, replaces the symlink by a copy of the file

git ls-files -s | awk '/120000/{print $4}' > symlinked.files
for F in $(< symlinked.files)
do
  TARGET="`dirname "$F"`"/"`cat "$F" | head -n 1`"
  if ! test -f "$TARGET"
  then
    echo "$F does not seem to be a borked symlink file"
  else
    rm "$F"
    cp "$TARGET" "$F"
    git update-index --assume-unchanged "$F"
  fi
done
rm symlinked.files
