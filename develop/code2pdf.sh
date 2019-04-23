find . -name "*.cc" | xargs enscript --color=1 -C -Ecpp -fCourier8 -o - | ps2pdf - code.pdf
