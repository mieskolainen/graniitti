# Code to pdf with syntax highlight and codelines
#
# Run with: source ./develop/code2pdf.sh

find ./src/ -name "M*.cc" | xargs enscript --color=1 -C -Ecpp -fCourier8 -o - | ps2pdf - code_cc.pdf
find ./include/ -name "M*.h" | xargs enscript --color=1 -C -Ecpp -fCourier8 -o - | ps2pdf - code_h.pdf
