# Code to pdf with syntax highlight and codelines
#
# Run with: source ./develop/code2pdf.sh

# makefile
COMMAND="--verbose --color=1 -C -fCourier8 -o -"

find Makefile | xargs enscript -Emakefile --header='$n|$%/$=|' $COMMAND | ps2pdf - code_makefile.pdf

# MadGraph to Graniitti conversion
find ./develop/MG2GRA/MG2GRA.py | xargs enscript -Epython --header='$n|$%/$=|' $COMMAND | ps2pdf - code_MG2GRA.pdf

# -----------------------------

# cc files (one can combine output of commands with (command1; command2; ... commandN)
(find ./src/ -name "M*.cc" & \
 echo ./src/Program/gr.cc) | xargs enscript -Ecpp --header='$n|$%/$=|' $COMMAND | ps2pdf - code_cc.pdf

# h files
find ./include/ -name "M*.h" | xargs enscript -Ecpp --header='$n|$%/$=|' $COMMAND | ps2pdf - code_h.pdf

# combine
#gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None -sOutputFile=merged.pdf code_h.pdf code_cc.pdf



