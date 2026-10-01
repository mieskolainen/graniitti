#!/bin/sh
#
# Clear all output

echo "Cleaning all output folders completely (all results gone) -- Are you sure? (write yes or no)"
read yn
case $yn in
    Yes|yes|Y|y )
        echo "Deleting ..."
        rm tmp/* -f -r
        rm runs/* -f -r
        rm figs/* -f -r
        rm output/* -f -r
        rm vgrid/* -f -r
        rm eikonal/* -f -r
        rm sudakov/* -f -r
        rm nuclear/* -f -r
        ;;
    No|no|N|n ) 
        echo "Aborted."
        ;;
    * ) 
        echo "Invalid input. Aborted."
        ;;
esac
