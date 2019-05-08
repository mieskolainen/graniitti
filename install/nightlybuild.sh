#!/bin/sh
#
# Run with: source ./install/nightlybuild.sh

echo '* GRANIITTI automated nightly build script *'

# ------------------------------------------------
# Set the manually here
NAME="GRANIITTI"
VERSION="0.38"
TYPE="beta"
DATE="08/05/2019"
UPDATE="spinor/vector/global tests"

BUILD_DATE=`date -R`
echo $BUILD_DATE
# ------------------------------------------------

# Fetch latest version, compile
git pull origin master
source ./install/setenv.sh
make superclean
make -j4 ROOT14=TRUE

# Test if compilation went fine
FILE=./bin/gr

if test -f "$FILE"; then
  echo "Build successful"
  BUILD_PASSING="true"
  cp ./install/img/passing.svg ./install/img/build-status.svg
else
  echo "Build not successful"
  BUILD_PASSING="false"
  cp ./install/img/failing.svg ./install/img/build-status.svg
fi

# Output message
MSG="{\n \"name\": \"$NAME\",\n \"version\": \"$VERSION\",\n \"type\": \"$TYPE\",\n \"date\": \"$DATE\",\n \"update\": \"$UPDATE\",\n \"build_passing\": $BUILD_PASSING,\n \"build_date\": \"$BUILD_DATE\"\n}\n"

printf "$MSG" > VERSION.json

# Push via SSH
git remote set-url origin git+ssh://git@github.com/mieskolainen/GRANIITTI.git
git add -A
git commit -m "nightly build"
git push origin master

# Change back to https
git remote set-url origin https://github.com/mieskolainen/GRANIITTI.git
