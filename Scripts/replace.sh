#! /bin/bash

## Small utility script used to search and replace text in all project files.

old=${1//|/\\|}
new=${2//|/\\|}
if [ -z "${old}" ]; then
    echo "usage: $0 <old> <new>" 1>&2
    exit 1
fi

if [[ `uname` == Darwin ]]; then
_replace()
{
  sed -i '' "s|$2|$3|g" "$1"
}
else
_replace()
{
  sed -i'' "s|$2|$3|g" "$1"
}
fi

for f in $(find "$PWD" -mindepth 1 -maxdepth 1 -type f ! -name '.*'); do
  _replace "$f" "$old" "$new"
done
for d in $(find "$PWD" -mindepth 1 -maxdepth 1 -type d ! -name .git ! -name ThirdParty); do
  for f in $(find "$d" -type f ! -name .git); do
    _replace "$f" "$old" "$new"
  done
done
