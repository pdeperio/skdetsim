#!/bin/csh -f
#set DS_T_TUNE = 1
set DS_T_TUNE = 0
set posIDe = $1
set wav    = $2
set yooabs = $3
set yooray = $4
set yoomie = $5
set DATE = $6
set OUTDIR = $7
set BDIR = $8
set NEVT = $9
set tailnamev = $10
set nphotons = $11
set sigbms = $12
set tailname = $tailnamev
if ( $tailname == -1 ) set tailname =
set trun = $13

set BYRAN1 = 341551
set BYRAN2 = 234
set BYRAN3 = 38948
#BDIR: Base DIR storing MCs : it shoud exist.
#OUTDIR will be set as  $OUTDIR/$DATE/data/pos$posIDe/$wav
#DATE skdetsim version yrmn
#card dir, file_dir

#set SKDETSIM = /home/byang/skrep/lsrsim.$DATE/$MACHINE/skdetsim
set SKDETSIM = ./skdetsim
set scr = $0
if ( ! -e $SKDETSIM ) then
  sh -c 'echo $0 : Error. NOT exisiting $1 1>&2' $scr $SKDETSIM
  exit(11)
endif
if ( ! -d $BDIR ) then
  sh -c 'echo $0 : Error. NOT exisiting $1 1>&2' $scr $BDIR
  exit(12)
endif

switch ("$wav")
case "337":
case "365":
case "400":
case "420":
case "409.6":
	breaksw
default:
	sh -c 'echo $0 : wrong wavelength $1  1>&2' $scr $wav
	exit(14)
endsw

#the following input parameters of each injector and wavelength must be tuned!!
switch ("$posIDe")
case "0":
	set x =  -35.3; set y =  777.7; set z = 1802.7; set cos1 =  0.0138591; set cos2 =  0.0381394; set cos3 = -0.999176;
	switch ("$wav")
	case "337":
		set nphoton = 2800
		set sigbm = 2100; breaksw
	case "365":
	 	set nphoton = 2800
		set sigbm = 3300; breaksw
	case "400":
	 	set nphoton = 7800
		set sigbm = 2500; breaksw
	case "420":
	 	set nphoton = 6000
		set sigbm = 2350; breaksw
	case "409.6":
	 	set nphoton = 5000
		set sigbm = 2000; breaksw
	endsw
	breaksw
case "1":
	set x =  -70.7; set y = -777.7; set z = 1802.7; set cos1 =  0.0126454; set cos2 =  0.0230219; set cos3 = -0.999654;
	switch ("$wav")
	case "337":
	 	set nphoton = 4500
		set sigbm = 2900; breaksw
	case "365":
	 	set nphoton = 7152
		set sigbm = 4000; breaksw
	case "400":
	 	set nphoton = 9500
		set sigbm = 3300; breaksw
	case "420":
	 	set nphoton = 7250
		set sigbm = 3300; breaksw
	case "409.6":
	 	set nphoton = 7250
		set sigbm = 2900; breaksw
	endsw
	breaksw
case "2":
	set x = 1490.73; set y =  768.14; set z = 1232.25; set cos1 = -0.880855; set cos2 = -0.473377; set cos3 =  0.00319347;
	switch ("$wav")
	case "337":
	 	set nphoton = 3690
		set sigbm = 3700; breaksw
	case "365":
	 	set nphoton = 5248
		set sigbm = 5200; breaksw
	case "400":
	 	set nphoton = 7860
		set sigbm = 4200; breaksw
	case "420":
	 	set nphoton = 6000
		set sigbm = 3800; breaksw
	case "409.6":
	 	set nphoton = 5400
		set sigbm = 4000; breaksw
	endsw
	breaksw
case "3":
	set x = 1490.73; set y =  768.14; set z =  595.95; set cos1 = -0.87501; set cos2 = -0.484104; set cos3 =  0.00120348;
	switch ("$wav")
	case "337":
	 	set nphoton = 5000
		set sigbm = 3450; breaksw
	case "365":
	 	set nphoton = 5902
		set sigbm = 4200; breaksw
	case "400":
	 	set nphoton = 9000
		set sigbm = 3600; breaksw
	case "420":
	 	set nphoton = 6900
		set sigbm = 3150; breaksw
	case "409.6":
	 	set nphoton = 8000
		set sigbm = 3300; breaksw
	endsw
	breaksw
case  "4":
	set x = 1490.73; set y =  768.14; set z =  -40.35; set cos1 = -0.886625; set cos2 = -0.462162; set cos3 = -0.0174184;
	switch ("$wav")
	case "337":
	 	set nphoton = 6700
		set sigbm = 3500; breaksw
	case "365":
	 	set nphoton = 4997
		set sigbm = 4000; breaksw
	case "400":
	 	set nphoton = 5100
		set sigbm = 3500; breaksw
	case "420":
	 	set nphoton = 7100
		set sigbm = 3200; breaksw
	case "409.6":
	 	set nphoton = 7750
		set sigbm = 3300; breaksw
	endsw
	breaksw
case "5":
	set x = 1490.73; set y =  768.14; set z = -605.95; set cos1 = -0.876555; set cos2 = -0.481151; set cos3 =  0.0121669;
	switch ("$wav")
	case "337":
	 	set nphoton = 5000
		set sigbm = 3800; breaksw
	case "365":
	 	set nphoton = 5018
		set sigbm = 4500; breaksw
	case "400":
	 	set nphoton = 9000
		set sigbm = 3800; breaksw
	case "420":
	 	set nphoton = 6900
		set sigbm = 3600; breaksw
	case "409.6":
	 	set nphoton = 7170
		set sigbm = 3500; breaksw
	endsw
	breaksw
case "6":
	set x =  1490.73; set y =   768.14; set z = -1242.25; set cos1 = -0.867773; set cos2 = -0.496951; set cos3 =  0.00304752;
	switch ("$wav")
	case "337":
	 	set nphoton = 4640
		set sigbm = 4400; breaksw
	case "365":
	 	set nphoton = 4985
		set sigbm = 4500; breaksw
	case "400":
	 	set nphoton = 8890
		set sigbm = 4200; breaksw
	case "420":
	 	set nphoton = 7000
		set sigbm = 4100; breaksw
	case "409.6":
	 	set nphoton = 6610
		set sigbm = 4200; breaksw
	endsw
	breaksw
case "7":
	set x =   -70.7; set y =   777.7; set z = -1802.7; set cos1 = -0.00894029; set cos2 = -0.00246342; set cos3 =  0.999956;
	switch ("$wav")
	case "337":
	 	set nphoton = 4330
		set sigbm = 4300; breaksw
	case "365":
	 	set nphoton = 4564
		set sigbm = 4100; breaksw
	case "400":
	 	set nphoton = 8000
		set sigbm = 4200; breaksw
	case "420":
	 	set nphoton = 6250
		set sigbm = 4100; breaksw
	case "409.6":
	 	set nphoton = 6500
		set sigbm = 4150; breaksw
	endsw
	breaksw
default :
	sh -c 'echo $0 : wrong injector ID $1 1>&2' $scr $posIDe
	exit(15)
endsw

#echo "nphoton" $nphoton
if ( $sigbms != -1 ) then
 @ sigbm = $sigbms
endif
if ( $nphotons != -1 ) then
 @ nphoton = $nphotons
endif

if ( $DS_T_TUNE != 0 ) then
  set DSOFFSET = CDSOFFSET
  set DSTRGDS  = (DS-TRGDS `grep $trun dst$wav.table | awk '{print $2}'`)
  set DSTRGOFF = CDSTRGOFF
endif


set DIR  =  $OUTDIR/$DATE/data/pos$posIDe/$wav
#checking and making out dir
if ( ! -d $DIR ) then
  if ( ! -d $OUTDIR/$DATE ) then
    mkdir $OUTDIR/$DATE
  endif
  if ( ! -d $OUTDIR/$DATE/data ) then
    mkdir $OUTDIR/$DATE/data
  endif
  if ( ! -d $OUTDIR/$DATE/data/pos$posIDe ) then
    mkdir $OUTDIR/$DATE/data/pos$posIDe
  endif
  mkdir $DIR
endif

if ( $yooabs != -1 || $yooray != -1 || $yoomie != -1 ) then
  if ( $yooabs != -1 ) then
    set abs = $yooabs
  else
    set abs = dabs
  endif
  if ( $yooray != -1 ) then
    set ray = $yooray
  else
    set ray = dray
  endif
  if ( $yoomie != -1 ) then
    set mie = $yoomie
  else
    set mie = dmie
  endif
  set out =  $DIR/{$wav}.mc.pos{$posIDe}.{$abs}_{$ray}_{$mie}{$tailname}.dat
  set card = $BDIR/{$wav}.pos{$posIDe}.{$abs}_{$ray}_{$mie}{$tailname}.card.$$.`hostname`
else
  set out =  $DIR/{$wav}.mc.pos{$posIDe}.default{$tailname}.dat
  set card = $BDIR/{$wav}.pos{$posIDe}.default{$tailname}.card.$$.`hostname`
endif

echo  ''
echo  out    : $out
echo  card   : $card
echo  yooabs : $yooabs
echo  yooray : $yooray
echo  yoomie : $yoomie
echo  nphoton: $nphoton
echo  sigbm  : $sigbm
echo  tail   : $tailname
echo  x, y, z : $x $y $z
echo  cos1 cos2 cos3 : $cos1 $cos2 $cos3
echo  ''

if ( ! -e xsneut.dat ) then
  ln -s $CERN_ROOT/lib/xsneut95.dat xsneut.dat
endif
if ( ! -e chetc.dat ) then
  ln -s /home/atmpd/const/chetc.dat .
endif

#if ( -e $out && ! -z $out ) then
  #sh -c 'echo $0 : Already $1 exists. Please check. 1>&2' $scr $out
  #exit(17)
#endif

## SUBRUN number is input from enviromment variable.
setenv MC_SUBNUM 0

#card writing start
#cat <<! >! $card

  sed -e 's/NEVENT/'$NEVT'/; s/BYRAN1/'$BYRAN1'/; s/BYRAN2/'$BYRAN2'/; s/BYRAN3/'$BYRAN3'/; s/VTXX/'$x'/; s/VTXY/'$y'/; s/VTXZ/'$z'/; s/DCOSX/'$cos1'/; s/DCOSY/'$cos2'/; s/DCOSZ/'$cos3'/' laser_sk4_example.card > $card


#card over

$SKDETSIM  $card $out $nphoton $wav $sigbm $yooabs $yooray $yoomie

unsetenv MC_SUBNUM

if ( ! -e $out || -z $out ) then
  sh -c 'echo $0 ": Error. NOT exisiting" $1 1>&2' $scr $OUT
  exit(17)
endif

exit(0)
