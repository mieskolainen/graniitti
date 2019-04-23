#!/usr/bin/env bash
#
# clang-tidy and clang-format based automatic style processing
#
# Use with: ./format.sh sourcefile.cc
# -------------------------------------------------------------------

# -------------------------------------------------------------------

##TIDY="clang-tidy"
##
##if [ -f /usr/local/opt/llvm/bin/clang-tidy ]; then
##    TIDY="/usr/local/opt/llvm/bin/clang-tidy"
##fi


##echo "clang-tidy processing "$*" ..."
##sleep 1
##
##$TIDY \
##    -fix \
##    -fix-errors \
##    -header-filter=.* \
##    --checks=readability-braces-around-statements,misc-macro-parentheses \
##    $*
##    #\
##    #-- -I.
##echo "clang-tidy done!"

# -------------------------------------------------------------------

FORMAT="clang-format"

if [ -f /usr/local/bin/clang-format ]; then
    FORMAT="/usr/local/bin/clang-format"
fi

echo "clang-format processing "$*" ..."
sleep 1

# -sort-includes
$FORMAT -i $*

echo "clang-format done!"
