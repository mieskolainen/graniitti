# Code to pdf with syntax highlight and codelines
#
# Run with: source ./develop/code2pdf.sh

COMMAND="--verbose --color=1 -C -Ecpp -fCourier8 -o -"

# cc files (one can combine output of commands with (command1; command2; ... commandN)
(find ./src/ -name "M*.cc" & \
 echo ./src/Program/gr.cc) | xargs enscript $COMMAND | ps2pdf - code_cc.pdf

# h files
find ./include/ -name "M*.h" | xargs enscript $COMMAND | ps2pdf - code_h.pdf

# combine
#gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None -sOutputFile=merged.pdf code_h.pdf code_cc.pdf
