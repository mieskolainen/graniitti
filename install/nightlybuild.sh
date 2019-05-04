#!/bin/sh
#
# Run with: source ./install/nightlybuild.sh

echo '* GRANIITTI automated nightly build script *'

# ------------------------------------------------
# Set the manually here
NAME="GRANIITTI"
VERSION="0.35"
TYPE="beta"
DATE="28/04/2019"
UPDATE="Algorithm updates"

BUILD_DATE=`date`
echo $BUILD_DATE
# ------------------------------------------------

# Fetch latest version, compile
git pull origin master && source ./install/setenv.sh && make superclean && make -j4

# Test if compilation went fine
FILE=./bin/gr

if test -f "$FILE"; then
  echo "Build successful"
  BUILD_SUCCESS="true"
  cp ./install/img/passing.svg ./install/img/build-status.svg
else
  echo "Build not successful"
  BUILD_SUCCESS="false"
  cp ./install/img/failing.svg ./install/img/build-status.svg
fi

# Output message
MSG="{\n  \"name\": \"$NAME\",\n \"version\": \"$VERSION\",\n \"type\": \"$TYPE\",\n \"date\": \"$DATE\",\n \"update\": \"$UPDATE\",\n \"build-success\": $BUILD_SUCCESS,\n \"build-date\": \"$BUILD_DATE\"\n}"

echo -e $MSG > VERSION.json

# Push via SSH
git remote set-url origin git+ssh://git@github.com/mieskolainen/GRANIITTI.git
git add -A
git commit -m "nightly build"
git push origin master

# Change back to https
git remote set-url origin https://github.com/mieskolainen/GRANIITTI.git
