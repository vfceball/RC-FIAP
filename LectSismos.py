inFileID = open(inFilename, 'r')
#set GMdirection 1;              # ground-motion direction
eqnms_listX = open('Sismos/nmsfileX.txt', 'r')
#set eqnms_listY     [read [open $nmsfileY "r"]];
set dts_list        [read [open "Sismos/dtsfile.txt" "r"]];
set durs_list       [read [open "Sismos/dursfile.txt" "r"]];
set nrecs           [llength $dts_list];
#puts "eqnms_listX=$eqnms_listX"
#set sisdir "Sismos"
source ReadSMDfile.tcl;     # procedure for reading GM file and converting it to proper format

for {set i 1} {$i <= $nrecs} {incr i 1} {

    #set IM_log [open $sisdir/IM_${i}.txt "w"];

    # Load the info about the record
    set EQnameX [lindex $eqnms_listX $i-1];     # Get the name of the record1
    #set EQnameY [lindex $eqnms_listY $i-1];     # Get the name of the record2
    set dt  [lindex $dts_list $i-1];    # Current dt
    set dur [lindex $durs_list $i-1];   # Current duration



    #set GMfile "H-e12140" ;         # ground-motion filenames
    set GMfact 1.;             # ground-motion scaling factor

    # set up ground-motion-analysis parameters
    #    set DtAnalysis  [expr 0.005*$sec];   # time-step Dt for lateral analysis
    #    set TmaxAnalysis    [expr 35.0*$sec];   # maximum duration of ground-motion analysis -- should be 50*$sec
    #

    #  ---------------------------------    perform Dynamic Ground-Motion Analysis
    # the following commands are unique to the Uniform Earthquake excitation
    #set IDloadTag 400;  # for uniformSupport excitation
    # Uniform EXCITATION: acceleration input
    set inFile Sismos/$EQnameX.at2
    set outFile Sismos/$EQnameX.g3;  # set variable holding new filename (PEER files have .at2/dt2 extension)
    ReadSMDFile $inFile $outFile dt;        # call procedure to convert the ground-motion file
#    set GMfatt [expr $g*$GMfact];       # data in input file is in g Unifts -- ACCELERATION TH
#    set AccelSeries "Series -dt $dt -filePath $outFile -factor  $GMfatt";   # time series information
#
#    pattern UniformExcitation  $IDloadTag  $GMdirection -accel  $AccelSeries  ;     # create Unifform excitation
#
#    set Nsteps [expr int($TmaxAnalysis/$DtAnalysis)];
#    set ok [analyze $Nsteps $DtAnalysis];           # actually perform analysis; returns ok=0 if analysis was successful
}

