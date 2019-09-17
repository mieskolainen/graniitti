#/bin/bash
#
# Energy dependence scanning test

# Path
P=./tests/run_minbias

# CMS-energies
E=1.995262E+01,2.241326E+01,2.517735E+01,2.828233E+01,3.177022E+01,3.568825E+01,4.008947E+01,4.503346E+01,5.058717E+01,5.682579E+01,6.383377E+01,7.170601E+01,8.054908E+01,9.048272E+01,1.016414E+02,1.141762E+02,1.282569E+02,1.440741E+02,1.618419E+02,1.818009E+02,2.042213E+02,2.294067E+02,2.576980E+02,2.894784E+02,3.251780E+02,3.652803E+02,4.103282E+02,4.609315E+02,5.177754E+02,5.816296E+02,6.533585E+02,7.339333E+02,8.244449E+02,9.261187E+02,1.040331E+03,1.168629E+03,1.312749E+03,1.474643E+03,1.656502E+03,1.860788E+03,2.090268E+03,2.348049E+03,2.637619E+03,2.962901E+03,3.328298E+03,3.738757E+03,4.199836E+03,4.717777E+03,5.299592E+03,5.953159E+03,6.687326E+03,7.512035E+03,8.438449E+03,9.479112E+03,1.064811E+04,1.196128E+04,1.343640E+04,1.509343E+04,1.695481E+04,1.904575E+04,2.139454E+04,2.403301E+04,2.699685E+04,3.032621E+04,3.406616E+04,3.826734E+04,4.298662E+04,4.828791E+04,5.424297E+04,6.093243E+04,6.844686E+04,7.688800E+04,8.637014E+04,9.702166E+04,1.089868E+05,1.224274E+05,1.375257E+05,1.544859E+05,1.735377E+05,1.949391E+05,2.189798E+05,2.459853E+05,2.763212E+05,3.103982E+05,3.486778E+05,3.916781E+05,4.399814E+05,4.942417E+05,5.551936E+05,6.236623E+05,7.005749E+05,7.869726E+05,8.840252E+05,9.930468E+05,1.115513E+06,1.253083E+06,1.407618E+06,1.581211E+06,1.776213E+06,1.995262E+06

# Scan
./bin/xscan -i $P/sd.json,$P/dd.json,$P/sd_xi_005.json,$P/dd_xi_005.json -l true -e $E

