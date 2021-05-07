#!/bin/sh
#
# Run with: source ./install/nightlybuild.sh

echo '* GRANIITTI automated nightly build script *'

# Fetch latest version, compile
git pull origin master
source ./install/setenv.sh
make superclean
make -j4 TEST=TRUE

# Test if compilation went fine
FILE=./bin/gr

# Run program so VERSION.json is constructed
./bin/gr

if test -f "$FILE"; then
  echo "Build successful"
  BUILD_PASSING="true"
  cp ./install/img/passing.svg ./install/img/build-status.svg
else
  echo "Build not successful"
  BUILD_PASSING="false"
  cp ./install/img/failing.svg ./install/img/build-status.svg
fi

# Build timestamp
BUILD_DATE=`date -R`
echo $BUILD_DATE

# Output message
MSG="{\n \"name\": \"GRANIITTI\",\n \"build_passing\": $BUILD_PASSING,\n \"build_date\": \"$BUILD_DATE\"\n}\n"
printf "$MSG" > ./develop/BUILD.json

# Push via SSH
git remote set-url origin git+ssh://git@github.com/mieskolainen/graniitti.git
git add -A
git commit -m "nightly build"
git push origin master

# Change back to https
git remote set-url origin https://github.com/mieskolainen/graniitti.git
