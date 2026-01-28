#!/usr/bin/perl -w 
use List::Util qw(min max) ;

# Compute the statistical similarity between the SMILES in the input file
#
# Also output the 10 compounds that are most similar to each other (10 NN)
#
# usage: recsim.pl [-options] list_of_smiles.smi [alternative_freq_mtx]
#
# options: -v  verbose output 
#          -d  debugging output. Also enables verbose output
#
# if alternative_freq_mtx is provided, those will be used instead of
# atpairfrq{1-4}.mtx
#
# Output is written to standard out in space separated format:
#
# (symmetic) matrix of simililarities between the compounds
# 10 NN compounds
#
#
# limitations:
#
# SMILES containing lower case letters (e.g. c1ccccc1) are not allowed
# isotopes (e.g. [13C]), except [2H] should be converted to their normal element
# inorganic elements and counter ions (e.g. Fe, Mg, Na, K) should be removed
#
#
$sep     = " " ; # separator (white space)
$MNAME   = "atpairfrq"; 
$option  = "xxxx" ;
$query   = "xxxx" ;
$mollist = "xxxx" ;
$verbose = 0 ;
$debug   = 0 ;
$matrix  = 0 ;
$ncomp   = 0 ;
$ncompmax = 250 ; # limit of accepted number of compounds
$reference = "xxxx" ;

# atom types
#
my %taty =("ZZ","0", "C4","1", "C3","2", "CO","3", "C2","4", "CP","5", "C+","6", "CZ","7", "CA","8", "CB","9", "CC","10", "CD","11", "CE","12", "N3","13", "N2","14", "N1","15", "NA","16", "NH","17", "NP","18", "NB","19", "NZ","20", "NO","21", "NC","22", "O2","23", "O1","24", "OF","25", "OC","26", "OE","27", "ON","28", "B3","29", "B4","30", "SI","31", "P","32", "P5","33", "S2","34", "S+","35", "SO","36", "S4","37", "SA","38", "F","39", "CL","40", "BR","41", "I","42", "H","43", "HO","44", "HN","45", "HV","46", "HB","47", "HX","48"); 
$maxt = 48 ; # maximum number of atom types
#$maxth = 38 ; # number without hydrogens and halogens as they cannot 
#               occur as middle atoms in angles and torsions
$maxt2 = $maxt * $maxt ;
$maxt3 = $maxt2 * $maxt ;

if (defined($ARGV[0])) {
  $option = $ARGV[0] ;     # the first argument of the command line
} 
else {
  print "No arguments given. Aborting\n" ;
  die ;
}

if (defined($ARGV[1])) {
  $mollist = $ARGV[1] ;
}

if (defined($ARGV[2])) {
  $MNAME = $ARGV[2] ;
}


if (substr($option,0,1) eq "-") { # options requested
  $option = substr($option,1) ;   # remove the dash
  $verbose = 0 ;
  $debug   = 0 ;
  $matrix  = 0 ;
  for ($i=0 ; $i < length($option) ; $i++) {
    if (substr($option,$i,1) eq "v") {$verbose = 1}
    if (substr($option,$i,1) eq "d") {
      $debug = 1 ;
      $verbose = 1 ;
    }
  }
  $sumo = $verbose + $debug ;
  if ($sumo > 0) {print "requested options: \n"} 
  if ($debug > 0) {
    print "d debugging output\n" ;
  }
  if ($verbose > 0) {
    print "v verbose output\n" ;
  }
}
else {  # no options
  if (defined($ARGV[1])) {
    $MNAME = $ARGV[1] ;
  }
  $mollist = $option ;
}

# print "$mollist $MNAME\n" ;

if (-e "$mollist") {
  if ($verbose > 0) {     
    print "list of SMILES $mollist found\n" ;
  } 
}
else { 
  if ($mollist eq "xxxx") {
    print "No SMILES given. Aborting\n" ; 
    die ;
  }
  print "list of SMILES $mollist not found\n" ; 
  die ;
}

if ($verbose > 0) {     
  print "checking if frequency matrices are present:\n" ;
}
for ($i=1 ; $i <= 4 ; $i++) {
  $tnam = $MNAME.$i.'.mtx' ;
  if (-e "$tnam") {
    if ($verbose > 0) {     
      print "matrix $tnam found\n" ;
    }
  }
  else {
    print "matrix $tnam not found. Aborting\n" ; 
    die ;
  }
}


open(OTHERS,$mollist) or die "error opening $mollist : $!\n"; 

while(defined($ci = <OTHERS>)) {
  @cols = split(/[\s\t]+/,$ci); # split at spaces or tabs

  $smi = " "  ;     # compound SMILES
  $cid = "" ;       # compound name

  if ($cols[0] ne "") {
    $smi = $cols[0] ;
  }

  if (defined($cols[1])) {
    if ($cols[1] ne "") {
      $cid = $cols[1] ;
    }
  }

  if (defined($cols[2])) {
    if ($cols[2] ne "") {      # splitted name, e.g. acetic ester
      if ($cols[2] =~ /[A-Z]/) {} # e.g. CHIRAL NATURAL WITHDRAWN
      else {
        $cid = $cid.$cols[2] ;
      }
    }
  }

  if ($verbose > 0) {     
    print "$smi $cid\n" ;
  }


# check for unallowed lower case letters
#
  if($smi =~ m/[cno]/) {
    print "$smi $cid lower case letters not allowed! Aborting\n" ;
    die ;
  }

# print "$smi $cid\n" ;

  $ncomp++ ;
  $smil[$ncomp] = $smi ; # SMILES
  $smin[$ncomp] = $cid ; # name of SMILES

}

close (OTHERS) or die "error closing $mollist : $!\n";

print "number of compounds $ncomp\n" ;

if ($ncomp > $ncompmax) {
  print "count of compounds is above a reasonable size ($ncompmax) !\n" ;
  print "aborting\n" ;
  die ;
}

# read in propability matrices
# $p11[i]  # propability of having atom type i   
# $p12[i]  # propability of for atom type j being 1st neighbor of atom type i
# index i is assigned as j >= i
# $p13[i]  # propability of for atom type k being 2nd neighbor of atom type i
# $p14[i]  # propability of for atom type l being 3rd neighbor of atom type i
#
for ($j=1 ; $j <= $maxt ; $j++) {
  $p11[$j] = 0.000000001 ;  # init with default
}
$tnam = $MNAME.'1'.'.mtx' ;

open(MATX,$tnam) or die "error opening $tnam : $!\n"; 

while(defined($ck = <MATX>)) {
  @cols = split(/[\s\t]+/,$ck); # split at spaces or tabs
  if (defined($cols[0]) & defined($cols[1])) {
    $ttt = $taty{$cols[0]} ; # convert string to numerical atom type
    $p11[$ttt] = $cols[1] ;  # propability of atom type
  }
}
close (MATX) or die "error closing $tnam : $!\n";

for ($j=1 ; $j <= $maxt ; $j++) {
  for ($k=$j ; $k <= $maxt ; $k++) {
    $dij = ($j-1) * $maxt + $k ; # index
    $p12[$dij] = 0.000000001 ;  # init with default
  }
}
$tnam = $MNAME.'2'.'.mtx' ;
open(MATX,$tnam) or die "error opening $tnam : $!\n"; 
while(defined($ck = <MATX>)) {
  @cols = split(/[\s\t]+/,$ck); # split at spaces or tabs
  if (defined($cols[0]) & defined($cols[1])) {
     @cls = split(/-/,$cols[0]); # split at dash sign
     $ttt = $taty{$cls[0]} ; # convert string to numerical atom type
     $uuu = $taty{$cls[1]} ; # 
     $pp = exp($cols[1]) ; # convert back from logarithmic value 
     if ($uuu >= $ttt) {
       $dij = ($ttt-1) * $maxt + $uuu ; # j>=i
     }
     else {
       $dij = ($uuu-1) * $maxt + $ttt ; # i>j
     }
     $p12[$dij] = $pp ;  # propability of i-j
#    print "cls0 $cls[0] ttt $ttt  cls1 $cls[1] uuu $uuu  value $cols[1] p12 $pp index $dij \n" ;
   }
 }
close (MATX) or die "error closing $tnam : $!\n";

for ($j=1 ; $j <= $maxt ; $j++) {
  for ($k=1 ; $k <= $maxt ; $k++) { # just to be sure to init all array elements
    for ($l=1 ; $l <= $maxt ; $l++) { 
      $dij = ($j-1) * $maxt +($k-1) * $maxt2 + $l ; # index
      $p13[$dij] = 0.000000001 ;  # init with default
    }
  }
}
$tnam = $MNAME.'3'.'.mtx' ;
open(MATX,$tnam) or die "error opening $tnam : $!\n"; 
while(defined($ck = <MATX>)) {
  $cols[0] = substr($ck,0,8) ;
  $cols[1] = substr($ck,8,8) ;
  if (defined($cols[0]) & defined($cols[1])) {
    @cls = split(/-/,$cols[0]); # split at dash sign
    if (substr($cls[0],1,1) eq " ") {$cls[0] = substr($cls[0],0,1) } # e.g. "F " as atom i
    if (substr($cls[2],1,1) eq " ") {$cls[2] = substr($cls[2],0,1) } # e.g. "F " as atom k
    $ttt = $taty{$cls[0]} ; # convert string to numerical atom type
    $uuu = $taty{$cls[1]} ; 
    $vvv = $taty{$cls[2]} ; 
    $pp = exp($cols[1]) ; # convert back from logarithmic value 
    if ($vvv >= $ttt) {
      $dij = ($ttt-1)*$maxt + ($uuu-1)*$maxt2 + $vvv ; # k>=i
    }
    else {
      $dij = ($vvv-1)*$maxt + ($uuu-1)*$maxt2 + $ttt ; # i>k
    }
    $p13[$dij] = $pp ;  # propability of i-j-k
#   print "cls0 $cls[0] ttt $ttt  cls1 $cls[1] uuu $uuu  cls2 $cls[2] vvv $vvv  value $cols[1]  pp13 $pp index $dij \n" ;
  }
}
close (MATX) or die "error closing $tnam : $!\n";

for ($j=1 ; $j <= $maxt ; $j++) { # ijkl = lkji, l>=i
  for ($k=1 ; $k <= $maxt ; $k++) {
    for ($l=1 ; $l <= $maxt ; $l++) { 
      for ($m=1 ; $m <= $maxt ; $m++) { # 
        $dij = ($j-1) * $maxt +($k-1) * $maxt2 + ($l-1)*$maxt3 + $m ; # index
        $p14[$dij] = 0.000000001 ;  # init with default
      }
    }
  }
}
$tnam = $MNAME.'4'.'.mtx' ;
open(MATX,$tnam) or die "error opening $tnam : $!\n"; 
while(defined($ck = <MATX>)) {
  $cols[0] = substr($ck,0,11) ;
  $cols[1] = substr($ck,11,8) ;
  if (defined($cols[0]) & defined($cols[1])) {
    @cls = split(/-/,$cols[0]); # split at dash sign
    if (substr($cls[0],1,1) eq " ") {$cls[0] = substr($cls[0],0,1) } # e.g. "F " as atom i
    if (substr($cls[3],1,1) eq " ") {$cls[3] = substr($cls[3],0,1) } # e.g. "F " as atom l
    $ttt = $taty{$cls[0]} ; # convert string to numerical atom type
    $uuu = $taty{$cls[1]} ; 
    $vvv = $taty{$cls[2]} ; 
    $www = $taty{$cls[3]} ; 
    $pp = exp($cols[1]) ; # convert back from logarithmic value 
    if ($www >= $ttt) { 
      $dij = ($ttt-1)*$maxt + ($uuu-1)*$maxt2 + ($vvv-1)*$maxt3 + $www ; # l>=i ijkl
    }
    else {
      $dij = ($www-1)*$maxt + ($vvv-1)*$maxt2 + ($uuu-1)*$maxt3 + $ttt ; # i<l lkji
    }
    $p14[$dij] = $pp ;  # propability of i-j-k-l
#    print "cls0 $cls[0] ttt $ttt  cls1 $cls[1] uuu $uuu  cls2 $cls[2] vvv $vvv  cls3 $cls[3] www $www  value $cols[1]  pp14 $pp  index $dij \n" ;
  }
}
close (MATX) or die "error closing $tnam : $!\n";

# loop iteratively over compounds to deterime 1 to 1 similarity
for ($mi=1 ; $mi <= $ncomp ; $mi++) {
  $smi = $smil[$mi] ; # query molecule
  $query = $smi ;
  smi2hin() ; # retrive information of query 
# copy array information of query compound
  $natq = $natoms ; # number of atoms in query compound
  for ($i=1 ; $i <= $natq ; $i++) {
#   $typq[$i]     = $typ[$i] ;     # atom type of atom i  (not used)
    $nsubstq[$i]  = $nsubst[$i] ;  # count of substituents of atom i
    $nanglesq[$i] = $nangles[$i] ; # count of angles involving atom i as first atom
    $ntorsq[$i]   = $ntors[$i] ;   # count of torsions involving atom i as first atom
  }
  proparray() ; # retrive information of query    
  for ($j=1 ; $j <= $natq ; $j++) {
    $qfq11[$j] = $qf11[$j] ; # $qf11[i] => $qfq11[i]
    for ($k=1 ; $k <= $nsubstq[$j] ; $k++) {
      $qfq12[$j][$k] = $qf12[$j][$k] ; # $qf12[i][j] => $qfq12[i][j]
    }
    for ($l=1 ; $l <= $nanglesq[$j] ; $l++) {
      $qfq13[$j][$l] = $qf13[$j][$l] ; # $qf13[i][j] => $qfq13[i][j]
    }
    for ($m=1 ; $m <= $ntorsq[$j] ; $m++) {
      $qfq14[$j][$m] = $qf14[$j][$m] ; # $qf14[i][j] => $qfq14[i][j]
    }
  }
# query finished
  for ($ni=$mi ; $ni <= $ncomp ; $ni++) {
    $smi = $smil[$ni] ; #  reference molecule
    $reference = $smi ;
#   print "query $mi   reference $ni" ;
    if ($mi != $ni) { # otherwise skip computation if query and reference are identical 
    smi2hin() ; # retrive information of reference
#   information of reference molecule is used directly 
#   $natoms     # number of atoms
#   $typ[i]     # atom type of atom i "NA"
#   $nsubst[i]  # number of 1st neighbors of atom i (including hydrogens) = number of bonds
#   $neib[i][j] # list of 1st neighbors of atom i
#   $tn12[i][j] # list of atom types of 1st neighbors, total of $nsubst[i] "NA-CA"
#   $nangles[i] # number of angles involving atom i (as first atom)
#   $tn13[i][j] # list of atom types of 2nd neighbors, total of $nangles[i] "NA-CA-CA"
#   $ntors[i]   # number of torsions involving atom i (as first atom)
#   $tn14[i][j] # list of atom types of 3rd neighbors, total of $ntors[i] "NA-CA-CA-BR"

    proparray() ; # retrive information of reference

#   now compute similarity

#   atom to atom matrix
    for ($j=1 ; $j <= $natq ; $j++) {
      for ($k=1 ; $k <= $natoms ; $k++) {
        $s11[$j][$k] = ($qfq11[$j] + $qf11[$k]) / (2.0*max($qfq11[$j],$qf11[$k])) ;
#       print "j $j  k $k  s11 $s11[$j][$k]\n";
        $scm[$j][$k] = $s11[$j][$k] ;  # holds the sum of all 4 matrices
      }
    }

#   first neighbors to frist neighbors matrix
    for ($j=1 ; $j <= $natq ; $j++) {
      $sq12 = 0 ;
      for ($jj=1 ; $jj <= $nsubstq[$j] ; $jj++) {
        $sq12 = $sq12 + $qfq12[$j][$jj] ;
      }
      for ($k=1 ; $k <= $natoms ; $k++) {
        $sr12 = 0 ;
        for ($kk=1 ; $kk <= $nsubst[$k] ; $kk++) {
          $sr12 = $sr12 + $qf12[$k][$kk] ;
        }
        $s12[$j][$k] = ($sq12 + $sr12) / (2.0*max($sq12,$sr12)) ;
#       print "j $j  k $k  s12 $s12[$j][$k]\n";
        $scm[$j][$k] = $scm[$j][$k] + $s12[$j][$k] ;
      }
    }

#   second neighbors to second neighbors matrix
    for ($j=1 ; $j <= $natq ; $j++) {
      $sq13 = 0 ;
      for ($jj=1 ; $jj <= $nanglesq[$j] ; $jj++) {
        $sq13 = $sq13 + $qfq13[$j][$jj] ;
      }
      for ($k=1 ; $k <= $natoms ; $k++) {
        $sr13 = 0 ;
        for ($kk=1 ; $kk <= $nangles[$k] ; $kk++) {
          $sr13 = $sr13 + $qf13[$k][$kk] ;
        }
        $s13[$j][$k] = ($sq13 + $sr13) / (2.0*max($sq13,$sr13)) ;
#       print "j $j  k $k  s13 $s13[$j][$k]\n";
        $scm[$j][$k] = $scm[$j][$k] + $s13[$j][$k] ;
      }
    }

#   third neighbors to third neighbors matrix
    for ($j=1 ; $j <= $natq ; $j++) {
      $sq14 = 0 ;
      for ($jj=1 ; $jj <= $ntorsq[$j] ; $jj++) {
        $sq14 = $sq14 + $qfq14[$j][$jj] ;
      }
      for ($k=1 ; $k <= $natoms ; $k++) {
        $sr14 = 0 ;
        for ($kk=1 ; $kk <= $ntors[$k] ; $kk++) {
          $sr14 = $sr14 + $qf14[$k][$kk] ;
        }
        $s14[$j][$k] = ($sq14 + $sr14) / (2.0*max($sq14,$sr14,0.000001)) ;
#       print "j $j  k $k  s14 $s14[$j][$k]\n";
        $scm[$j][$k] = $scm[$j][$k] + $s14[$j][$k] ;
      }
    }

#   compute similarity
    $simt  = 0 ; # total similarity 
    $matat = 0 ; # count of matched atoms 
    for ($j=1 ; $j <= $natq ; $j++) { 
      $matq[$j] = 0 ;
      for ($k=1 ; $k <= $natoms ; $k++) {
        $pas[$j][$k] = 0 ; # pairwise atom to atom assignment 
      }
    }
    for ($k=1 ; $k <= $natoms ; $k++) {
      $matr[$k] = 0 ;
    }
    for ($j=1 ; $j <= $natq ; $j++) {
      $matq[$j] = 0 ;
    }
    if ($natoms > $natq) { # query molecule is smaller
      while ($matat < $natq) {
        $ttt = 0 ;
        for ($j=1 ; $j <= $natq ; $j++) { # get overall maximum
          if ($matq[$j] > 0) {next} # query atom already matched
          for ($k=1 ; $k <= $natoms ; $k++) {
            if ($matr[$k] > 0) {next} # reference atom already matched
            if ($scm[$j][$k] > $ttt) { # maximum pairwise value
              $ttt = $scm[$j][$k] ;
	      $ttj = $j ; # query atom
              $ttk = $k ; # matched reference atom
            }
          }
        }   
        $simt = $simt + $ttt ; 
        $pas[$ttj][$ttk] = $ttt / 4.0 ; # match score of atoms 
        $matq[$ttj] = 1 ;
        $matr[$ttk] = 1 ;
#       print "j $ttj  ttk $ttk  pas $pas[$ttj][$ttk]\n" ;
#       if ($matrix > 0) {print "matat $matat  query $ttj  reference $ttk $pas[$ttj][$ttk]\n"} 
        $matat++ ;
      } # end while
    }
    else { 
      while ($matat < $natoms) {
        $ttt = 0 ;
        for ($j=1 ; $j <= $natoms ; $j++) { # get overall maximum
          if ($matr[$j] > 0) {next} # refernce atom already matched
          for ($k=1 ; $k <= $natq ; $k++) {
            if ($matq[$k] > 0) {next} # query atom already matched
            if ($scm[$k][$j] > $ttt) { # maximum pairwise value
              $ttt = $scm[$k][$j] ;
              $ttk = $k ; # matched query atom
              $ttj = $j ; # reference atom
            }
          }
        }   
        $simt = $simt + $ttt ; 
        $pas[$ttk][$ttj] = $ttt / 4.0 ; # match score of atoms 
        $matq[$ttk] = 1 ;
        $matr[$ttj] = 1 ;
        $matat++ ;
      } # end while
    }
  
    $simt = $simt / 4.0 ; # normalize (4 matrices)

#   similarity <= 1 
#   $simtha = 0.5 * $simt * (1.0 / $natq + 1.0 / $natoms) ; # harmonic mean
    
    $minn = min($natoms,$natq) ; 
    $maxn = max($natoms,$natq) ;
	  
#   $simtha = $simt / ($minn + 0.5 * log($maxn - $minn + 1)) ; # log scaling
    $simtha = $simt / ($minn + log($maxn - $minn + 1)) ; # log scaling2

    }
    else {  # query and reference molecule are identical
      $simtha   = 1.0 ; 
    }

#   printf ("%7.4f",$simtha) ;
#   print "\n" ;

#   store similarity matrix elements
    $simqq[$ni][$mi] = $simtha ;
    $simqq[$mi][$ni] = $simtha ;

  } # end inner loop
} # end outer loop

# print similarity matrix
print "similarity" ;
for ($mi=1 ; $mi <= $ncomp ; $mi++) {
  print "$sep$smin[$mi]" ;
}
print "\n" ;
for ($mi=1 ; $mi <= $ncomp ; $mi++) {
  print "$smil[$mi]" ;
  for ($ni=1 ; $ni <= $ncomp ; $ni++) {
    printf ("%7.4f",$simqq[$mi][$ni]) ;
  }
  print "\n" ;
}

print "\n" ;

# determine 10 NN as sum of ranks or scores, respectively
#
for ($mi=1 ; $mi <= $ncomp ; $mi++) {
  $tmax = -9.9 ; # maximum similarity in this row
  $tmmm = 1 ;    # molecule that has maximum similarity
  $tsum[$mi] = 0.0 ;  # sum of similarities in this row
  $ttmn[$mi] = 0 ; # count of molecule index
  for ($ni=1 ; $ni <= $ncomp ; $ni++) {
    if ($ni == $mi) {next} # skip idenitical molecules
    if ($simqq[$mi][$ni] > $tmax) {
      $tmax = $simqq[$mi][$ni] ;
      $tmmm = $ni ;
    }
    $tsum[$mi] = $tsum[$mi] + $simqq[$mi][$ni] ;
  }
  $ttmm[$mi] = $tmmm ; # molecule that has maximum similarity in this row
  $ttms[$mi] = $tmax ; # maximum similarity in this row
  $tavs[$mi] = $tsum[$mi] / ($ncomp - 1) ; # average similarity of molecule
}

$ttav = 0 ;
for ($mi=1 ; $mi <= $ncomp ; $mi++) {
  $ttmn[$ttmm[$mi]]++ ; # count of occurence as maximum
  $ttav = $ttav + $tavs[$mi] ;
  print "$mi av $tavs[$mi] sum $tsum[$mi]  max $ttmm[$mi] $ttms[$mi]\n" ;
}
$ttav = $ttav / $ncomp ; # total average of similarity between all molecules

$sse = 0 ;
for ($mi=1 ; $mi <= $ncomp ; $mi++) {
  $sse = $sse + ($tavs[$mi] - $ttav) * ($tavs[$mi] - $ttav) ;
}
$sse = $sse / ($ncomp - 1) ; # variance
$sse = sqrt($sse) ; # standard deviation
print "average similarity $ttav  standard deviation $sse\n" ;


# determine cluster using DBSCAN
#
$factor = int($ncomp / 10.0) +1 ;
$minpts = max($factor,4) ; # at least 4 neighbors in a cluster
#$eps = 0.9400 ; # maximum distance, minimum similarity
$eps = min($ttav + 2.5 * $sse, 0.98) ; # maximum distance, minimum similarity

$nclust = 0 ; # number of clusters

for ($i=1; $i<=$ncomp ; $i++) {
  if (defined($label[$i])) {
    next ; # already labelled 
  }
  else {
    $nn = 0 ; # number of neighbors 
    for ($j=1; $j<=$ncomp ; $j++) { # find neighbors of i
      if ($j == $i) {next}
      if ($simqq[$i][$j] >= $eps) {$nn++}
    }
    if ($nn < $minpts) {
      $label[$i] = "N" ; # point is noise (not directly reachable)
      next ; 
    }
    $nclust++ ;
    $label[$i] = $nclust ; 
  }
}

print "\n" ;
print "# cluster compound max-count\n" ;
for ($mi=1 ; $mi <= $ncomp ; $mi++) {
  print "$mi $label[$mi] $smin[$mi] $ttmn[$mi]\n" ;
}


exit(0) ;

# end of main program
# 

sub smi2hin {
# 
# setup connectivities, atomtypes, and neighborlists
#
# input: $smi  SMILES string
#
# determines:
# $natoms     # number of atoms
# $elem[i]    # element of atom i
# $typ[i]     # atom type of atom i "NA"
# $nsubst[i]  # number of 1st neighbors of atom i (including hydrogens) = number of bonds
# $neib[i][j] # list of 1st neighbors of atom i

# $tn12[i]    # list of atom types of 1st neighbors, total of $nsubst[i] "NA-CA"
# $nangles[i] # number of angles involving atom i (as first atom)
# $tn13[i]    # list of atom types of 2nd neighbors, total of $nangles[i] "NA-CA-CA"
# $ntors[i]   # number of torsions involving atom i (as first atom)
# $tn14[i]    # list of atom types of 3rd neighbors, total of $ntors[i] "NA-CA-CA-BR"

$dative = 1 ;  # dative bonds
$neutra = 1 ;  # neutralize acidic oxygens
$aromat = 1 ;  # force aromatic bond information in .hin file
$natoms = 0           ;  # total number of atoms
$nheavy = 0           ;  # number of non-hydrogen atoms
$nhydr = 0            ;  # number of hydrogen atoms

# substitute [NH] to N
    $smi =~ s/\[NH\]/N/g ;

# substitute [2H] to H
    $smi =~ s/\[2H\]/H/g ;

#   print "$smi\n" ;

#   now convert SMILES to .hin
#   interprete SMILES: get each character build up connection matrix
    $natoms = 0      ;  # total number of atoms
    $nheavy = 0      ;  # number of non-hydrogen atoms
    $nhydr = 0       ;  # number of hydrogen atoms
    $single = "s"    ;
    $double = "d"    ;
    $triple = "t"    ;
    $aromatic = "a"  ;
    $pos = 0         ;  # current position / atom  in SMILES string
    $bat = 0         ;  # last atom position (also from where branching occurs)
    $cm1 = " "       ; 

# init arrays
    for ($i = 0; $i <= length($smi) ; $i++) {
       $neighbors[$i] = "" ;
       $el[$i] = "" ;    # atom position in SMILES string
       $elem[$i] = "" ;  # atoms numbered consequtively
       $emap[$i] = -1 ;  # map conseq. atoms position to pos. in SMILES
       $chrg[$i] = 0 ;
       $nsub[$i] = 0 ;   # number of substituents of each atom
       $bondlist[$i] = "" ;
    }
# number of ring systems that are open at the same time
# is limited here to 9 !
    for ($i = 0; $i <= 9  ; $i++) { 
       $ring[$i] = -1 ;
    }

    for ($j=1 ; $j <= length($smi) ; $j++) {
      $c1 = substr($smi,$j-1,1) ;   # one character (e.g. C, N, O)
      if ($j < length($smi)) {
        $c2 = substr($smi,$j,1) ;
      }
      else { 
        $c2 = "" ;                  # next character (e.g. Cl, Br, Si)
      }


# atom or branching or bond or ring opening/closing ?
# element
# still to do: treat stereoinformation @ and @@
      if ($c1 eq "@") {           # skip @ and @@
        $pos++ ;
        $cm2 = $cm1 ; # pre-preceeding character
        $cm1 = $c1 ; # preceeding character
        next ;
      }
      if ($c1 eq "H" && $cm1 eq "@") {
        $pos++ ;
        $cm2 = $cm1 ; # pre-preceeding character
        $cm1 = $c1 ; # preceeding character
        next ;
      }

# still to do: treat cis/trans information / and \


      if ($c1 =~ /[a-z]/) {next}  # skip lower case letters  
      if ($c1 eq "[" && $pos == 0) {next}  # if SMILES starts with [
      if ($c1 =~ /[A-Z]/) {       # element (always starts with upper case letter)
        $nheavy++;
        if ($c2 =~ /[a-z]/) {
          $el[$pos] = $c1.$c2 ;
        }
        else {
          $el[$pos] = $c1 ;
        }
#       $neighbors[$pos] = "$pos" ;  # obsolte !!!
        if ($pos > 0) {  # bonded to previous atom
          $neighbors[$bat] = $neighbors[$bat]." ".$pos ; 
          $neighbors[$pos] = $neighbors[$pos]." ".$bat ; 
          $nsub[$pos]++ ;
          $nsub[$bat]++ ;
          if ($c2 =~ /\+/) { $chrg[$pos] =  1 }
          if ($c2 =~ /\-/) { $chrg[$pos] = -1 }
          if ($cm1 eq "[") {$cm1 = $cm2} # cases like C=[N+] 
          if ($cm1 eq "=") {
            $bondlist[$bat] = $bondlist[$bat]." ".$double ; 
            $bondlist[$pos] = $bondlist[$pos]." ".$double ; 
          }
          elsif ($cm1 eq "#") {
            $bondlist[$bat] = $bondlist[$bat]." ".$triple ; 
            $bondlist[$pos] = $bondlist[$pos]." ".$triple ; 
          }
          elsif ($cm1 eq ":") {
            $bondlist[$bat] = $bondlist[$bat]." ".$aromatic ; 
            $bondlist[$pos] = $bondlist[$pos]." ".$aromatic ; 
          }
          else {
            $bondlist[$bat] = $bondlist[$bat]." ".$single ; 
            $bondlist[$pos] = $bondlist[$pos]." ".$single ; 
          }
        }
        $bat = $pos ;
#       if ($c2 =~ /\+/) { $chrg[$nheavy] =  1 }
#       if ($c2 =~ /\-/) { $chrg[$nheavy] = -1 }
        if ($c2 =~ /[a-z]/) { $pos++ }
      }

# open branch
# the branching atom is put on the stack
      if ($c1 =~ /[(]/) {
        push(@blist,$bat) ;
##      print "$c1 $bat\n" ;
      }

# close branch
# the latest branching atom is retrieved from the stack
      if ($c1 =~ /[)]/) {
        $bat = pop(@blist) ;
##      print "$c1 $bat\n" ;
      }

# open/close ring
      if ($c1 =~ /[1-9]/) {
        if ($ring[$c1] == -1) {  # ring opening
          $ring[$c1] = $bat ;
          $nsub[$bat]++ ;
#         print " open ring system $c1 at position $bat\n" ;
        }
        else {                   # ring closure, connect atoms
          $neighbors[$bat] = $neighbors[$bat]." ".$ring[$c1] ; 
          $neighbors[$ring[$c1]] = $neighbors[$ring[$c1]]." ".$bat ; 
          $nsub[$bat]++ ;
#         print "close ring system $c1 between $bat and $ring[$c1]\n" ;
##        print "bat $bat  ring $ring[$c1]\n" ;
          if ($cm1 eq "[") {$cm1 = $cm2} # cases like C=[N+] 
          if ($cm1 eq "=") {
            $bondlist[$bat] = $bondlist[$bat]." ".$double ; 
            $bondlist[$ring[$c1]] = $bondlist[$ring[$c1]]." ".$double ; 
          }
          elsif ($cm1 eq "#") {
            $bondlist[$bat] = $bondlist[$bat]." ".$triple ; 
            $bondlist[$ring[$c1]] = $bondlist[$ring[$c1]]." ".$triple ; 
          }
          elsif ($cm1 eq ":") {
            $bondlist[$bat] = $bondlist[$bat]." ".$aromatic ; 
            $bondlist[$ring[$c1]] = $bondlist[$ring[$c1]]." ".$aromatic ; 
          }
          else {
            $bondlist[$bat] = $bondlist[$bat]." ".$single ; 
            $bondlist[$ring[$c1]] = $bondlist[$ring[$c1]]." ".$single ; 
          }

          $ring[$c1] = -1 ; # reset to allow re-use of this ring number
        }
      }

# explicit bonds
#     if ($c1 =~ /[=#:]/) {
##       print "explicit bond $c1\n" ;
#     }

      $cm2 = $cm1 ; # pre-preceeding character
      $cm1 = $c1 ; # preceeding character
      $pos++ ;

   }

# now derive atoms sequentially 
# setup connectivity matrix $neib[][], add hydrogens
# $elem[] $typ[] $charge[] $typ[] $nsubst[] $neib[][] $bond[][]

for ($i = 0; $i < $pos ; $i++) {
  if ($el[$i] =~ /[A-Z]/) {
    $natoms++ ; 
    if ($el[$i] ne "H") {
      $nheavy++ ;
    }
    else {
      $nhydr++ ;
    }
    $emap[$i] = $natoms ;    # atom position in SMILES to conseq. number
    $elem[$natoms] = $el[$i] ;
    $typ[$natoms] = "**" ;
    $charge[$natoms] = $chrg[$i] ;
    $nsubst[$natoms] = $nsub[$i] ;
#   split bond and neighbor lists
    @coln = split(/\s+/,$neighbors[$i]) ;
    @colb = split(/\s+/,$bondlist[$i]) ;
    for ($j = 1; $j <= $nsubst[$natoms] ; $j++) {
##     print "$j $coln[$j] $colb[$j]\n" ;
      $neib[$natoms][$j] = $coln[$j] ;
      $bond[$natoms][$j] = $colb[$j] ;
    }
  }
}

# now correct atom positions in neighbor list 
for ($i=1 ; $i <= $natoms ; $i++) {
  for ($j = 1; $j <= $nsubst[$i] ; $j++) {
    $neib[$i][$j] = $emap[$neib[$i][$j]] ;
  }
}


# Determine valencies and add hydrogens
$newato = $natoms ;
for ($i=1 ; $i <= $natoms ; $i++) {
  $neb = 0 ;  # number of bonds 
  $val = 0 ;  # number of valencies
  if ($nsubst[$i] > 0) {
    for ($k=1 ; $k <= $nsubst[$i] ; $k++) {
      if ($bond[$i][$k] eq "s" ) {
        $neb = $neb + 1 ;
      }
      elsif ($bond[$i][$k] eq "d") {
        $neb = $neb + 2 ;
      }
      elsif ($bond[$i][$k] eq "t") {
        $neb = $neb + 3 ;
      }
      else {   # aromatic
        $neb = $neb + 1 ;
      }
    }
  }
  if ($elem[$i] eq "B" ) {
    if ($charge[$i] == -1) {
      $val = 4 ;
      $typ[$i] = "B4" ;
    }
    else {
      $val = 3 ;
      $typ[$i] = "B3" ;
    }
  }
  if ($elem[$i] eq "C" ) {
    if ($charge[$i] == -1) {
      $val = 3 ;
    }
    else {
      $val = 4 ;
    }
  }
  if ($elem[$i] eq "Si") {
    $val = 4 ;
    $typ[$i] = "SI" ;
  }
  if ($elem[$i] eq "O" ) {
    if ($charge[$i] == -1) {
      $val = 1 ;
#     $val = 0 ;
    }
    elsif ($charge[$i] == 1) {
      $val = 3 ;
    }
    else {
      $val = 2 ;
    }
    if ($neb == 1 && $dative > 0) {  # no hydrogens added to [N+][O-] 
      $nn = $neib[$i][1] ;
      if ($elem[$nn] eq "N" && $nsubst[$nn] == 3) {
        $val = 1 ;
        $typ[$i] = "O1" ;
#       print "atom $i $elem[$i] $typ[$i] $charge[$i] $nsubst[$i]\n" ;
      }
    }
    if ($charge[$i] == -1 && $neutra > 0 && $typ[$i] eq "**") {
      $nn = $elem[$neib[$i][1]] ;
      if ($nn eq "C" || $nn eq "S" || $nn eq "P" ) {
        $charge[$i] = 0 ;   # neutralize acidic oxygens
        $val = 2 ;
      }
    }
  }
  if ($elem[$i] eq "F" ) {
    $val = 1 ;
    $typ[$i] = "F" ;
  }
  if ($elem[$i] eq "Cl") {
    $val = 1 ;
    $typ[$i] = "CL" ;
  }
  if ($elem[$i] eq "Br") {
    $val = 1 ;
    $typ[$i] = "BR" ;
  }
  if ($elem[$i] eq "I" ) {
    $val = 1 ;
    $typ[$i] = "I" ;
  }

#   Nitrogen, Sulfur and Phoshorus are treated individually
  if ($elem[$i] eq "N" ) {
    if ($neb > 3) {
      $val = $neb ;
      $charge[$i] = 1 ;
    }
#   if ($neb == 2) {      # wouldn't add hydogens to -N-
#     $charge[$i] = -1 ;
#     $val = 2 ;
#   }
    if ($charge[$i] == 0) {
      $val = 3 ;
    }
    if ($dative > 0) {   # convert dative to directional bonds
# dative nitro groups to double bonds   [N+]([O-])=O  to N(=O)=O
      $oneib   = 0 ; # number of oxygens as substituents
      $onei[1] = 0 ; 
      $onei[2] = 0 ; 
      $onei[3] = 0 ; 
      $nneib   = 0 ; # number of nitrogens as substituents
      $nnei[1] = 0 ;
      $nnei[2] = 0 ;
      $nnei[3] = 0 ;
      for ($k=1 ; $k <= $nsubst[$i] ; $k++) {
        $nn = $neib[$i][$k] ;
        if ($elem[$nn] eq "O") {
          $oneib++ ;
          $onei[$oneib] = $nn ;
        }
        if ($elem[$nn] eq "N") {
          $nneib++ ;
          $nnei[$oneib] = $nn ;
        }
      }
      if ($nsubst[$i] == 3) {
        if ($oneib == 2 || $oneib == 3) {  # NO2 
          $charge[$i] = 1 ;
          $typ[$i] = "NO" ;   # MM3 assigns NO
          for ($k=1 ; $k <= $oneib ; $k++) {
            $typ[$onei[$k]] = "O1" ;
#   print "atom $onei[$k] $charge[$onei[$k]] $nsubst[$onei[$k]]\n" ;
          }
          $val = $neb ;
        }
        if ($oneib == 1) {  # N+O-
          if ($neb == 3) {  
            $typ[$i] = "N2" ;
          }
          else {
            $typ[$i] = "N1" ;
          }
          $charge[$i] = 1 ;
          $typ[$onei[1]] = "O1" ;
          $charge[$onei[1]] = -1 ;
          $val = $neb ;
        }
      }
    }
  }

  if ($elem[$i] eq "S") {
    if ($nsubst[$i] < 2) {
      if ($neb == 1) {
#       $val = 1 ;  # didn't add hydrogen to SX2
        $val = 2 ;  # same as for -S-
      }
      else {
        $val = 0 ;  # =S
      }
    }
    elsif ($nsubst[$i] == 2) {
      if ($neb == 3) {
        $val = 1 ;
      }
      else {
        $val = 0 ;
      }
    }
    elsif ($nsubst[$i] == 3) {
      if ($neb == 5) {
        $val = 1 ;
      }
      else {
        $val = 0 ;
      }
    }
    elsif ($nsubst[$i] == 4) {
      $val = 0 ;
    }
    else {
      $val = 0 ;
    }
  }

  if ($elem[$i] eq "P") {
    if ($neb < 4) {   
      if ($nsubst[$i] < 3) {   #
        $val = 3 ;   # PX3
        $typ[$i] = "P" ;
      }
      else {
        $val = 5 ;   # PX5
        $typ[$i] = "P5" ;
      }
    }
    elsif ($neb == 4 && $nsubst[$i] == 4) {  # [P+]X4 
      $val = 4 ;
      $charge[$i] = 1 ; 
      $typ[$i] = "P5" ;
    }
    else {  # PX5
      $val = 5 ;
      $typ[$i] = "P5" ;
    }
  }


  if ($neb > 0) {
    $val = $val - $neb ;
  }
  else {
    $val = 0 ;
  }

#   $val is the number of hydrogens to be added
#  print "atom $i $elem[$i] nsub $nsubst[$i] bonds $neb  addhyd $val\n" ;

# add hydrogens
  if ($val > 0) {
    for ($j=1 ; $j <= $val ; $j++) {
      $newato++ ;               # number of new hydrogen atom
      $elem[$newato] = "H" ;
      $typ[$newato] = "**" ;
      $charge[$newato] = 0 ;
      $nsubst[$newato] = 1 ;
      $neib[$newato][1] = $i ;
      $bond[$newato][1] = "s" ;
      $nsubst[$i]++ ;           # update heavy atom
      $neib[$i][$nsubst[$i]] = $newato ;
      $bond[$i][$nsubst[$i]] = "s" ;
    }
  }

}
$natoms = $newato ; # update number of atoms

# generate topological bond matrix
for ($i=1 ; $i <= $natoms ; $i++) {
  for ($j=1 ; $j <= $natoms ; $j++) {
    $nbo[$i][$j] = 0 ;
  }
}
for ($i=1 ; $i <= $natoms ; $i++) {
  $nn = $nsubst[$i] ;
  for ($j=1 ; $j <= $nn ; $j++) {
    $kk = $neib[$i][$j] ;
    $nbo[$i][$kk] = $bond[$i][$j] ; # entries: s,d,t,a
  }
}

# ring perception is done according to Figuera's approach
# see J.Chem.Inf.Comput.Sci. 36 (1996) pp.986-991
# init arrays
  $nr3 = 0 ;  # number of 3-membered rings
  $nr4 = 0 ;  # number of 4-membered rings
  $nr5 = 0 ;  # number of 5-membered rings
  $nr6 = 0 ;  # number of 6-membered rings
  $nra = 0 ;  # number of aromatic rings
for ($i=1 ; $i <= $natoms ; $i++) {
  $ring3[$i] = 0 ;  # atom is part of a 3-membered ring
  $ring4[$i] = 0 ;  # atom is part of a 4-membered ring
  $ring5[$i] = 0 ;  # atom is part of a 5-membered ring
  $ring6[$i] = 0 ;  # atom is part of a 6-membered ring
  # Hier noch moeglich groesse Ringe (7,8) fuer 3D-Generierung
  $aring[$i] = 0 ;  # atom is part of a Hueckel aromatic ring
  $nring3[$i][0] = 0 ; # number and atoms of 3-membered rings
  $nring4[$i][0] = 0 ; # number and atoms 
  $nring5[$i][0] = 0 ; # number and atoms 
  $nring6[$i][0] = 0 ; # number and atoms 
  $atl[$i]   = 0 ;  # atom is not terminating
  $atpath[$i][0]  = 0 ;  # path of each atom 
  $atplen[$i] = 0     ;  # length of each atom path
}

# 1st generate list of non-terminating atoms (speed up ring search)
for ($i=1 ; $i <= $natoms ; $i++) {
  if ($elem[$i] eq "H") {next}  
  if ($elem[$i] eq "F") {next}  
  if ($elem[$i] eq "Cl") {next}  
  if ($elem[$i] eq "Br") {next}  
  if ($elem[$i] eq "I") {next}  
  if ($nsubst[$i] <= 1) {next}
  $hneib = 0 ;
  for ($j=1 ; $j <=$nsubst[$i] ; $j++) {
    $kk = $neib[$i][$j] ;
    if ($elem[$kk] eq "H") {$hneib++} 
  }
  # -CH3, -NH2, -OH
  if ($elem[$i] eq "C" && $hneib >= 3) {next}
  if ($elem[$i] eq "N" && $hneib >= 2) {next}
  if ($elem[$i] eq "O" && $hneib >= 1) {next}
  $atl[$i] = 1 ; # not terminating atom or group
} 

# put all non-terminating atoms into the queue
#
##print "\nnon-terminating atoms are " ;
$atqul = 0 ; # length of atom queue array
for ($i=1 ; $i <= $natoms ; $i++) {
  if ($atl[$i] == 1) {
    $atqul++ ;
    push(@ntqu,$i) ;
##    print "$i " ;
  }
}

for ($rn=1 ; $rn <= $atqul ; $rn++) { # loop over all non-terminating atoms
  @atqu = () ;   # non-terminating atoms (array of type queue)
  @source = () ; # preceeding atom (array of type queue)
  $t1 = shift(@ntqu) ;
  push(@atqu,$t1) ;
  for ($t2=1 ; $t2 <= $natoms ; $t2++) {
    $atpath[$t2][0]  = 0 ;  # path of each atom 
    $atplen[$t2] = 0     ;  # length of each atom path
  }
  $frsrc = 0 ;
  $nxtr  = 0 ;
  for ($i=1 ; $i <= $atqul ; $i++) {  # breadth-first-search loop
    $front = shift(@atqu) ;   # step 1; get front node and its source
##   print "current front node is atom $front\n" ;
    if ($i == 1) { # init path of first atom
      $atpath[$front][0] = $front ;
      $atplen[$front] = 1 ;
    }
    else {$frsrc = shift(@source)} 
##   print "current source atom is $frsrc\n" ;
    for ($j=1 ; $j <= $nsubst[$front] ; $j++) { # step 2
      $m = $neib[$front][$j] ;
      if ($atl[$m] == 0) {next} # skip terminating atoms
      if ($m == $frsrc) {next}  # atom connected to front node and =! source
      if ($atpath[$m][0] == 0) {     
##       print "compute path for atom $m is " ; # path empty, compute path
        for ($k=0 ; $k < $atplen[$front] ; $k++) { 
          $atpath[$m][$k] = $atpath[$front][$k] ;
        } 
        $atpath[$m][$atplen[$front]] = $m ; # add atom m to path
        $atplen[$m] = $atplen[$front] + 1;
        for ($k=0 ; $k < $atplen[$m] ; $k++) {
##         print " $atpath[$m][$k]" ;
        }
##       print "\n" ;
        push(@atqu,$m) ;       # add atom $m to atom queue
        push(@source,$front) ; # add front node of atom $m to queue
##       print "adding atom $m to queue\n" ;
##       print "adding atom $front to source queue\n" ;
      }
      else { # compute intersection of path[front node] and path[atom m]
##       print "path to atom $m exits\n" ;
        $isectl = 0 ; # count of intersecting atoms
        for ($k1=0 ; $k1 < $atplen[$m] ; $k1++) {
          for ($k2=0 ; $k2 < $atplen[$front] ; $k2++) {
            if ($atpath[$m][$k1] == $atpath[$front][$k2]) { 
              $isectl++;
            }
          }
        }
        if ($isectl == 1) { # one atom in common => ring closure
#          print "ring closure between atoms $front and $m\n" ;
          # compute new ring set containing atoms of path[m]+path[front node]
          @newlt = () ; # init with empty array
          for ($k1=0 ; $k1 < $atplen[$m] ; $k1++) {
            push(@newlt, $atpath[$m][$k1]) ;
          }
          for ($k2=0 ; $k2 < $atplen[$front] ; $k2++) {
            push(@newlt, $atpath[$front][$k2]) ;
          }
          @newl1 = sort {$a <=> $b} @newlt ; # sort common atom path
##         print "common path is @newl1\n" ;
#         remove double atoms
          $rsize = 0 ;
          $atold = 0 ;
##          print "ring contains " ;
          for ($k1=0 ; $k1 < $atplen[$m]+$atplen[$front] ; $k1++) {
            $tt = shift(@newl1) ;
            if ($tt != $atold) {
              $newring[$rsize] = $tt ;
              $rsize++ ;
            }
            $atold = $tt ;
          }
          for ($k=0 ; $k < $rsize ; $k++) {
##            print " $newring[$k]"
          }
#         print "\n" ;       
          # assign atoms to corresponding ring systems
          # and check if this ring was already found 
          if ($rsize < 7) {
            if ($rsize == 3) {
              for ($k1=0 ; $k1 < 3 ; $k1++) {
                $tt = $newring[$k1] ;
                $ring3[$tt] = 1 ;
              }
              # loop over all existing 3-membered rings
              $atold = 0 ;
              for ($k1=1 ; $k1 <= $nr3 ; $k1++) { 
                $atold = 0 ;
                for ($l1=0 ; $l1 < 3 ; $l1++) {
                  $tt = $newring[$l1] ;
                  for ($l2=0 ; $l2 < 3 ; $l2++) {
                    if ($tt == $nring3[$k1][$l2]) {$atold++} 
                  }
                }
                if ($atold >= 3) {last} # ring already present
              }
              if ($atold < 3) { # ring not present yet
                $nr3++ ;
##                print "3-membered ring no $nr3 contains" ;
                for ($l1=0 ; $l1 < 3 ; $l1++) {
                  $nring3[$nr3][$l1] = $newring[$l1] ;
##                  print " $newring[$l1]" ;
                }
##                print "\n" ;
              }
            }
            if ($rsize == 4) {
              for ($k1=0 ; $k1 < 4 ; $k1++) {
                $tt = $newring[$k1] ;
                $ring4[$tt] = 1 ;
              }
              # loop over all existing 4-membered rings
              $atold = 0 ;
              for ($k1=1 ; $k1 <= $nr4 ; $k1++) { 
                $atold = 0 ;
                for ($l1=0 ; $l1 < 4 ; $l1++) {
                  $tt = $newring[$l1] ;
                  for ($l2=0 ; $l2 < 4 ; $l2++) {
                    if ($tt == $nring4[$k1][$l2]) {$atold++} 
                  }
                }
                if ($atold >= 4) {last} # ring already present
              }
              if ($atold < 4) { # ring not present yet
                $nr4++ ;
##                print "4-membered ring no $nr4 contains" ;
                for ($l1=0 ; $l1 < 4 ; $l1++) {
                  $nring4[$nr4][$l1] = $newring[$l1] ;
##                  print " $newring[$l1]" ;
                }
##                print "\n" ;
              }
            }
            if ($rsize == 5) {
              for ($k1=0 ; $k1 < 5 ; $k1++) {
                $tt = $newring[$k1] ;
                $ring5[$tt] = 1 ;
              }
              # loop over all existing 5-membered rings
              $atold = 0 ;
              for ($k1=1 ; $k1 <= $nr5 ; $k1++) { 
                $atold = 0 ;
                for ($l1=0 ; $l1 < 5 ; $l1++) {
                  $tt = $newring[$l1] ;
                  for ($l2=0 ; $l2 < 5 ; $l2++) {
                    if ($tt == $nring5[$k1][$l2]) {$atold++} 
                  }
                }
                if ($atold >= 5) {last} # ring already present
              }
              if ($atold < 5) { # ring not present yet
                $nr5++ ;
##              print "5-membered ring no $nr5 contains" ;    #remove pr
                for ($l1=0 ; $l1 < 5 ; $l1++) {
                  $nring5[$nr5][$l1] = $newring[$l1] ;
##                print " $newring[$l1]" ;                    #remove pr
                }
##              print "\n" ;                                  #remove pr
              }
            }
            if ($rsize == 6) {
              for ($k1=0 ; $k1 < 6 ; $k1++) {
                $tt = $newring[$k1] ;
                $ring6[$tt] = 1 ;
              }
              # loop over all existing 6-membered rings
              $atold = 0 ;
              for ($k1=1 ; $k1 <= $nr6 ; $k1++) { 
                $atold = 0 ;
                for ($l1=0 ; $l1 < 6 ; $l1++) {
                  $tt = $newring[$l1] ;
                  for ($l2=0 ; $l2 < 6 ; $l2++) {
                    if ($tt == $nring6[$k1][$l2]) {$atold++} 
                  }
                }
                if ($atold >= 6) {last} # ring already present
              }
              if ($atold < 6) { # ring not present yet
                $nr6++ ;
##                print "6-membered ring no $nr6 contains" ;
                for ($l1=0 ; $l1 < 6 ; $l1++) {
                  $nring6[$nr6][$l1] = $newring[$l1] ;
##                  print " $newring[$l1]" ;
                }
##                print "\n" ;
              }
            }
          }
          $nxtr = 1 ;
          last ; # start with next atom as root
        }
        if ($nxtr == 1) {last} 
      }
      if ($nxtr == 1) {last} 
    }
    if ($nxtr == 1) {last} # start from next non-terminating atom
  }
# print "###\n" ;
} # end of loop over all non-terminating atoms

# assign aromatic rings  
# $aring[$i] = 1 ;  # atom is part of a Hueckel aromatic ring
#
# Hueckel aromatic rings: only 5- and 6-membered rings 
# check for  
# 6-membered ring: 3 double bonds of C and N, but not exocyclic
# 5-membered ring: 2 double bonds of C and N and
#                    O with two neighbors or
#                    S with two neighbors
#
# This approach does not work for case such as
#         
#   //\ //\   Here, only the right ring is recognizes as
#   |  |  ||  being aromatic, since the left one possess
#   \\/ \\/   two exocyclic double bonds.
#     

for ($i=1 ; $i <= $nr6 ; $i++) { # check all 6-membered rings
  $dbc = 0 ; # double bonds of carbons
  $dbn = 0 ; # double bonds of nitrogens
  for ($j=0 ; $j < 6 ; $j++) { 
    $tt = $nring6[$i][$j] ;
    if ($elem[$tt] eq "O") {last} # these elements do not occur 
    if ($elem[$tt] eq "S") {last} # in 6-membered aromatic rings
    if ($elem[$tt] eq "P") {last}
    if ($elem[$tt] eq "Si") {last}
    if ($elem[$tt] eq "C") {
      if ($nsubst[$tt] != 3) {
        last ;
      }
      else { 
        for ($k1=1 ; $k1 <= 3 ; $k1++) {
          $nn = $neib[$tt][$k1] ; # double bond to neighbor in the same ring
          # double bonds are counted twice this way !
          if ($bond[$tt][$k1] eq "d" || $bond[$tt][$k1] eq "a") { 
            for ($k2=0 ; $k2 < 6 ; $k2++) {
              if ($nn == $nring6[$i][$k2]) {$dbc++}
            }
          }
        }
      }
    }
    if ($elem[$tt] eq "N") {
      if ($nsubst[$tt] != 2) {
        last ;
      }
      else { 
        for ($k1=1 ; $k1 <= 2 ; $k1++) {
          $nn = $neib[$tt][$k1] ; # double bond to neighbor in the same ring
          # double bonds are counted twice this way !
          if ($bond[$tt][$k1] eq "d" || $bond[$tt][$k1] eq "a") { 
            for ($k2=0 ; $k2 < 6 ; $k2++) {
              if ($nn == $nring6[$i][$k2]) {$dbn++}
            }
          }
        }
      }
    }
  }
# print "6-membered ring no $i contains $dbc C= and $dbn N= bonds\n" ;
  if ($dbc+$dbn == 6) { # aromatic ring  (double bonds counted twice)
    $nra++ ;
    for ($l1=0 ; $l1 < 6 ; $l1++) {
      $tt = $nring6[$i][$l1] ;
      $aring[$tt] = 1 ;  
    }
    if ($aromat == 1) {
      for ($l1=0 ; $l1 < 6 ; $l1++) { # update bonding information
        $tt = $nring6[$i][$l1] ;
        for ($k1=1 ; $k1 <= $nsubst[$tt] ; $k1++) {
          $nn = $neib[$tt][$k1] ;
          for ($k2=0 ; $k2 < 6 ; $k2++) {
            if ($nn == $nring6[$i][$k2]) {
              $nbo[$tt][$nring6[$i][$k2]] = "a" ;
              $bond[$tt][$k1] = "a" ;
            }
          }
        }
      }
    }
  }
}

for ($i=1 ; $i <= $nr5 ; $i++) { # check all 5-membered rings
  $dbc = 0 ; # double bonds of carbons
  $dbn = 0 ; # double bonds of nitrogens
  $dbo = 0 ; # count of oxygens and sulfurs with two substituents
  for ($j=0 ; $j < 5 ; $j++) { 
    $tt = $nring5[$i][$j] ;
    if ($elem[$tt] eq "P")  {last} # these elements do not occur 
    if ($elem[$tt] eq "Si") {last} # in 6-membered aromatic rings
    if ($elem[$tt] eq "C") {
      if ($nsubst[$tt] != 3) {
        last ;
      }
      else { 
        for ($k1=1 ; $k1 <= 3 ; $k1++) {
          $nn = $neib[$tt][$k1] ; # double bond to neighbor in the same ring
          # double bonds are counted twice this way !
          if ($bond[$tt][$k1] eq "d" || $bond[$tt][$k1] eq "a") { 
            for ($k2=0 ; $k2 < 5 ; $k2++) {
              if ($nn == $nring5[$i][$k2]) {$dbc++}
            }
          }
        }
      }
    }
    if ($elem[$tt] eq "N") {
      if ($nsubst[$tt] != 2) {
        last ;
      }
      else { 
        for ($k1=1 ; $k1 <= 2 ; $k1++) {
          $nn = $neib[$tt][$k1] ; # double bond to neighbor in the same ring
          # double bonds are counted twice this way !
          if ($bond[$tt][$k1] eq "d" || $bond[$tt][$k1] eq "a") { 
            for ($k2=0 ; $k2 < 5 ; $k2++) {
              if ($nn == $nring5[$i][$k2]) {$dbn++}
            }
          }
        }
      }
    }
    if ($elem[$tt] eq "O" || $elem[$tt] eq "S") {
      if ($nsubst[$tt] == 2) {$dbo++}  
    }
  }
# print "5-membered ring no $i contains $dbc C= and $dbn N= bonds\n" ;
# print "                     and $dbo oxygens and sulfurs\n" ;
  if ($dbc+$dbn == 4 && $dbo == 1) { # aromatic ring  
    $nra++ ;
    for ($l1=0 ; $l1 < 5 ; $l1++) {
      $tt = $nring5[$i][$l1] ;
      $aring[$tt] = 1 ;  
    }
    if ($aromat == 1) {
      for ($l1=0 ; $l1 < 5 ; $l1++) { # update bonding information
        $tt = $nring5[$i][$l1] ;
        for ($k1=1 ; $k1 <= $nsubst[$tt] ; $k1++) {
          $nn = $neib[$tt][$k1] ;
          for ($k2=0 ; $k2 < 5 ; $k2++) {
            if ($nn == $nring5[$i][$k2]) {
              $nbo[$tt][$nring5[$i][$k2]] = "a" ;
              $bond[$tt][$k1] = "a" ;
            }
          }
        }
      }
    }
  }
}

# check ring assignment   
#print "\natom  3ring  4ring  5ring  6ring  aromat\n" ;
#for ($k=1 ; $k <= $natoms ; $k++) {
#  print "$k\t$ring3[$k]\t$ring4[$k]\t$ring5[$k]\t$ring6[$k]\t$aring[$k]\n" ;
#}
#print "\n" ;


for ($i=1 ; $i <= $natoms ; $i++) {
  if ($typ[$i] ne "**") {next}      # type already assigned
  $nn = $nsubst[$i] ;               # number of neighbors
  # The sequence is ordered in decreasing frequency of elements 
  # found in molecules to speed up things somewhat

  if ($elem[$i] eq "C" && $nn > 0) {                        # carbon
    if ($nn == 4) {                                         # sp3
      $typ[$i] = "C4" ;
      if ($ring4[$i] == 1) {$typ[$i] = "CB"} # cyclobutane
      if ($ring3[$i] == 1) {$typ[$i] = "CP"} # cyclopropane
      # atom type for 3-membered ring superseeds larger rings
      next ;
    }
    if ($nn == 3) {                                         # sp2
      $typ[$i] = "C3" ;
      if ($ring4[$i] == 1) {$typ[$i] = "CC"} # cyclobutene
      if ($ring3[$i] == 1) {$typ[$i] = "CZ"} # cyclopropene
      $oneib = 0 ;         # >C=O ?
      $nneib = 0 ;
      for ($j=1 ; $j <= $nn ; $j++) {
        $kk = $neib[$i][$j] ;
        if ($elem[$kk] eq "O" && $bond[$i][$j] eq "d") {
          $oneib++ ;
        }
        if ($elem[$kk] eq "N") {$nneib++}
      }
      if ($oneib > 0) {
        $typ[$i] = "CO" ;
        if ($ring4[$i] == 1) {               # cyclobutanone
          if ($nneib < 1) {$typ[$i] = "CD"}  # except lactim
        }
        if ($ring3[$i] == 1) {               # cyclopropanone
          if ($nneib < 1) {$typ[$i] = "CE"}  # except lactim
        }
        next ;
      }
      if ($aring[$i] == 1) {  # aromatic 
        $typ[$i] = "CA" ;
        next ;
      }
    }
    if ($nn == 2) {                                         # sp
      $typ[$i] = "C2" ;
      next ;
    }
    if ($nn == 1) {   # sp isonitrile [N+]#[C-]
      $typ[$i] = "C2" ;
      next ;
    }
    if ($charge[$i] == 1) {  # carbonium
      $typ[$i] = "C+" ;
      next ;
    }
    if ($typ[$i] eq "**") {  # default
#     $typ[$i] = "C?" ;   # C? is not defined 
      $typ[$i] = "C4" ;
      next ;
    } 
  }

  if ($elem[$i] eq "O" && $nn > 0) {                        # oxygen
    if ($nn == 2) { 
      $typ[$i] = "O2" ;
      if ($ring3[$i] == 1) {
        $typ[$i] = "OE" ;      # expoxide
        next ;
      }
      if ($ring5[$i] == 1) {
        $oneib = 0;
        for ($j=1 ; $j <= $nn ; $j++) {
          $kk = $neib[$i][$j] ;
          if ($elem[$kk] eq "C" && $nsubst[$kk] == 3) {
            $oneib++ ;
          }
        }
#       check if aromatic ring is furane (OF) or contains N or S (O2)
        if ($oneib == 2 && $aring[$i] == 1) {  # Car-O-Car
          for ($j=1 ; $j <= $nr5 ; $j++) {
            $m5 = 0 ;
            for ($k=0 ; $k < 5 ; $k++) {
              $kk = $nring5[$j][$k] ;
              if ($kk == $i) {$m5++}
              if ($kk == $neib[$i][1]) {$m5++}
              if ($kk == $neib[$i][2]) {$m5++}
            }
            if ($m5 == 3) {  # O and its neighbors are in ring $j
              $m5 = $j ;
              last ;
            }
          }
          $nns = 0 ;
          for ($j=0 ; $j < 5 ; $j++) {
            $kk = $nring5[$m5][$j] ;
            if ($elem[$kk] eq "N" || $elem[$kk] eq "S") {$nns++}
          }
          if ($nns == 0) {
            $typ[$i] = "OF" ;  # furane
          }
          else {
            $typ[$i] = "O2" ;  # other heteroaromatic 5-ring
          }
          next ;
        }
      }
    }
    if ($nn == 1) { 
      $kk = $neib[$i][1] ;
      if ($elem[$kk] eq "C" || $elem[$kk] eq "S" || $elem[$kk] eq "P") {
        if ($bond[$i][1] eq "d") {              # acidic oxygens
          $typ[$i] = "O1" ; # carbonyl oxygen
        }
        else {
          $typ[$i] = "OC" ; # carboxylate oxygen
        }
        next ;
      }
      if ($elem[$kk] eq "N") {
        if ($bond[$i][1] eq "s") {
          $typ[$i] = "ON" ;   # amine oxide oxygen
        }
        else {
          $typ[$i] = "O1" ;   # N=O, N(=O)(=O) 
        }
        next ;
      }
    }
    if ($typ[$i] eq "**") {  # default
#     $typ[$i] = "O?" ; # O? is not defined
      $typ[$i] = "O2" ;
      next ;
    } 
  }

  if ($elem[$i] eq "S" && $nn > 0) {                        # sulfur
    if ($nn == 1) {
      $typ[$i] = "S2" ;
      next ;
    }
    if ($nn == 2) {
      $typ[$i] = "S2" ;
      $oneib = 0;
      for ($j=1 ; $j <= $nn ; $j++) {
        $kk = $neib[$i][$j] ;
        if ($elem[$kk] eq "C" && $nsubst[$kk] == 3) {
          $oneib++ ;
        }
      }
      if ($oneib == 2 && $aring[$i] == 1) {  # Car-S-Car
        $typ[$i] = "SA" ;  # thiophen
        next ;
      }
    }
    if ($nn == 3) {
      $oneib = 0;
      for ($j=1 ; $j <= $nn ; $j++) {
        $kk = $neib[$i][$j] ;
        if ($elem[$kk] eq "O" && $nsubst[$kk] == 1) {
          $oneib++ ;
        }
      }
      if ($oneib > 0) {  
        $typ[$i] = "SO" ;  # >S=O
        next ;
      }
    }
    if ($nn == 4) {
      $typ[$i] = "S4" ;  # >S(=O)(=O)
      next ;
    }
    if ($typ[$i] eq "**") {  # default
#     $typ[$i] = "S?" ; # S? is not defined
      $typ[$i] = "S4" ;
      next ;
    } 
  }

  if ($elem[$i] eq "N" && $nn > 0) {                        # nitrogen
    if ($nn == 4) {
#     $typ[$i] = "N4" ;  # ammonium  
      $typ[$i] = "NH" ;  # ammonium  
      next ;
    }
    if ($nn == 1) {
      $typ[$i] = "N1" ;  # sp1
      next ;
    }
    if ($nn == 2 || $nn ==3) {
      if ($nn == 3) {$typ[$i] = "N3"}  # default 
      if ($nn == 2) {$typ[$i] = "N2"}  # default
      $hneib = 0 ;
      $oneib = 0 ;
      $nneib = 0 ;
      $cneib = 0 ;
      $sneib = 0 ;
      for ($j=1 ; $j <= 3 ; $j++) { 
        $hnei[$j] = 0 ; 
        $onei[$j] = 0 ; 
        $nnei[$j] = 0 ; 
        $cnei[$j] = 0 ; 
        $snei[$j] = 0 ;
      }
      for ($j=1 ; $j <= $nn ; $j++) { # determine kind of substituents
        $kk = $neib[$i][$j] ;
        if ($elem[$kk] eq "H") {
          $hneib++ ;
          $hnei[$hneib] = $kk ;
        }
        if ($elem[$kk] eq "O") {
          $oneib++ ;
          $onei[$oneib] = $kk ;
        }
        if ($elem[$kk] eq "N") {
          $nneib++ ;
          $nnei[$nneib] = $kk ;
        }
        if ($elem[$kk] eq "C") {
          $cneib++ ;
          $cnei[$cneib] = $kk ;
        }
        if ($elem[$kk] eq "S") {
          $sneib++ ;
          $snei[$sneib] = $kk ;
        }
      }
      if ($nn == 3) { 
        if ($hneib > 0) {
          if ($cneib == 1) {
            if ($nsubst[$cnei[1]] > 3) {
              $typ[$i] = "N3" ; # nonplanar amine or amide
              next ;
            }
            else {
              $typ[$i] = "N2" ; # planar amine or amide
              next ;
            }
          }
          if ($cneib > 1) {
            if ($cneib == 2 && $hneib == 1) {
              if ($nsubst[$cnei[1]] == 3 || $nsubst[$cnei[2]] == 3) {
                $typ[$i] = "N2" ; # planar amine 
#               check for pyrrole ring
                if ($ring5[$i] == 1) {
                  for ($j=1 ; $j <= $nr5 ; $j++) {
                    $m5 = 0 ;
                    for ($k=0 ; $k < 5 ; $k++) {
                      $kk = $nring5[$j][$k] ;
                      if ($kk == $i) {$m5++}
                      if ($kk == $cnei[1]) {$m5++}
                      if ($kk == $cnei[2]) {$m5++}
                    }
                    if ($m5 == 3) {  # N and its C neighbors are in ring $j
                      $m5 = $j ;
                    last ;
                    }
                  }
                  $nos = 0 ;
                  for ($j=0 ; $j < 5 ; $j++) {
                    $kk = $nring5[$m5][$j] ;
                    if ($elem[$kk] eq "O" || $elem[$kk] eq "S") {$nos++}
                  }
                  if ($nos == 0) {
                    $typ[$i] = "NP" ;   # pyrrole
                    next ;
                  }
                }
              }
              next ;
            }
            for ($k=1 ; $k <= $cneib ; $k++) { # determine kind of substituents
              $tt = $cnei[$k] ;
              if ($nsubst[$tt] == 3) {
                for ($k1=1 ; $k1 <= 3 ; $k1++) {
                  $oo = $neib[$tt][$k1] ;
                  if ($elem[$oo] eq "O" && $bond[$tt][$k1] eq "d") {
                    $typ[$i] = "N2" ; # peptide
                    next ;
                  }
                }
              }
            }
          }
          if ($sneib > 0) {
            if ($nsubst[$snei[1]] > 2) {
              $typ[$i] = "N2" ;
              next ;
            }
            if ($sneib == 2) {
              if ($nsubst[$snei[1]] == 2 && $nsubst[$snei[2]] == 2) {
                $typ[$i] = "N3" ;
                next ;
              }
              else {
                $typ[$i] = "N2" ;
                next ;
              }
            }
          }
          if ($nneib > 0) {
            $typ[$i] = "N2" ; 
            next ;
          }

        }
        else {
          if ($cneib > 0) {
            $typ[$i] = "N2" ; 
            next ;
          }
        }
      }
      if ($nn == 2) { 
        if ($oneib == 2) {     # NO2
          $typ[$i] = "NO" ;    # nitro 
          next ;
        }   
        if ($bond[$i][1] eq "t" && $bond[$i][2] eq "s") {
          $typ[$i] = "N1" ;    #   #[N-]  e.g. isonitrile
          next ;
        }
        if ($bond[$i][1] eq "s" && $bond[$i][2] eq "t") {
          $typ[$i] = "N1" ;    #   [N-]#[C+]  isonitrile
          next ;
        }
        if ($nneib == 2) {     
          if ($bond[$i][1] eq "d" && $bond[$i][2] eq "d") {
            # middle atom in -N=N=N 
            $typ[$i] = "N1" ;  # Hyperchem assigns N1
            next ; 
          }
          else {
            $typ[$i] = "NA" ;    # middle atom in -N=N-N= 
            next ;
          }
        }
        if ($cneib == 2) {     
          if ($bond[$i][1] eq "d" && $bond[$i][2] eq "s") {
            if ($nsubst[$cnei[2]] == 4) {
              $typ[$i] = "NC" ;    # imine C=N-C
              next ;
            }
            elsif ($nsubst[$cnei[2]] == 3) {
              $typ[$i] = "N2" ;    # C=N-Car 
              next ;
            }
            else {
              $typ[$i] = "NA" ;    
              next ;
            }
          }
          if ($bond[$i][1] eq "s" && $bond[$i][2] eq "d") {
            if ($nsubst[$cnei[1]] == 4) {
              $typ[$i] = "NC" ;    # imine C-N=C
              next ;
            }
            elsif ($nsubst[$cnei[1]] == 3) {
              $typ[$i] = "N2" ;    # Car-N=C    
              next ;
            }
            else {
              $typ[$i] = "NA" ;    
              next ;
            }
          }

#         check for pyrrole if explicit double bonds are present
          if ($ring5[$i] == 1) {
#   print "nsub N $nsubst[$i] C1 $nsubst[$cnei[1]] C2 $nsubst[$cnei[2]] \n" ;
            if ($nsubst[$cnei[1]] == 3 && $nsubst[$cnei[2]] == 3 &&                             $nsubst[$i] == 3) {
              for ($j=1 ; $j <= $nr5 ; $j++) {
                $m5 = 0 ;
                for ($k=0 ; $k < 5 ; $k++) {
                  $kk = $nring5[$j][$k] ;
                  if ($kk == $i) {$m5++}
                  if ($kk == $cnei[1]) {$m5++}
                  if ($kk == $cnei[2]) {$m5++}
                }
                if ($m5 == 3) {  # N and its C neighbors are in ring $j
                  $m5 = $j ;
                last ;
                }
              }
              $nos = 0 ;
              for ($j=0 ; $j < 5 ; $j++) {
                $kk = $nring5[$m5][$j] ;
                if ($elem[$kk] eq "O" || $elem[$kk] eq "S") {$nos++}
              }
              if ($nos == 0) {
                $typ[$i] = "NP" ;   # pyrrole
                next ;
              }
            }
          }

          if ($aring[$cnei[1]] == 1 && $aring[$cnei[2]] == 1) {
            if ($ring5[$i] == 1) {
#             check if pyrolle (NP) or other heteroaromatic ring (NA)
              for ($j=1 ; $j <= $nr5 ; $j++) {
                $m5 = 0 ;
                for ($k=0 ; $k < 5 ; $k++) {
                  $kk = $nring5[$j][$k] ;
                  if ($kk == $i) {$m5++}
                  if ($kk == $cnei[1]) {$m5++}
                  if ($kk == $cnei[2]) {$m5++}
                }
                if ($m5 == 3) {  # N and its C neighbors are in ring $j
                  $m5 = $j ;
                  last ;
                }
              }
              $nos = 0 ;
              for ($j=0 ; $j < 5 ; $j++) {
                $kk = $nring5[$m5][$j] ;
                if ($elem[$kk] eq "O" || $elem[$kk] eq "S") {$nos++}
              }
              if ($nos == 0) {
                $typ[$i] = "NP" ;   # pyrrole
              }
              else {
                $typ[$i] = "NA" ;   # other aromatic nitrogen
              }
              next ;
            }
          }
        }
        if ($nneib == 1 && $oneib == 1) {
          if ($nbo[$i][$nnei[1]] eq "d" && $nbo[$i][$onei[1]] eq "s") {
            $typ[$i] = "NB" ;  # middle atom in azoxy -N=N-O-
            next ;
          }
        }
        if ($cneib == 1 && $oneib == 1) {
          if ($nbo[$i][$cnei[1]] eq "d" && $nbo[$i][$onei[1]] eq "s") {
            $typ[$i] = "NC" ;  # oxime -C=N-OH 
            next ;
          }
          if ($nbo[$i][$cnei[1]] eq "a" && $nbo[$i][$onei[1]] eq "a") {
            $typ[$i] = "NC" ;  # aromatic oxazole -C=N-O-
            next ;
          }
        }
        if ($cneib == 1 && $nneib == 1) {
          if ($nbo[$i][$cnei[1]] eq "d" && $nbo[$i][$nnei[1]] eq "s") {
            $typ[$i] = "NA" ;
#           $typ[$i] = "N2" ;
            next ;
          }
          if ($nbo[$i][$cnei[1]] eq "s" && $nbo[$i][$nnei[1]] eq "d") {
            $typ[$i] = "NA" ;
#           $typ[$i] = "N2" ;
            next ;
          }
          if ($nbo[$i][$cnei[1]] eq "a" && $nbo[$i][$nnei[1]] eq "a") {
            $typ[$i] = "NA" ;
            next ;
          }
        }
      }
    }
    if ($typ[$i] eq "**") {  # default
#     $typ[$i] = "N?" ; # not defined
      $typ[$i] = "N3" ;
      next ;
    } 
  }

}

# Now assign hydrogen atoms types
for ($i=1 ; $i <= $natoms ; $i++) {
  if ($typ[$i] ne "**") {next}      # type already assigned
  $nn = $nsubst[$i] ;               # number of neighbors
  if ($elem[$i] eq "H" && $nn > 0) {                        # hydrogen
    $kk = $neib[$i][1] ;    # neighbor        
    $neibel = $elem[$kk] ;  # element of neighbor
    if ($neibel ne "N" && $neibel ne "O") {
      $typ[$i] = "H" ;    # H attached to C,Si,B,P,S,halogens
      next ;
    }
    if ($neibel eq "N") {           # connected to nitrogen     
      $typ[$i] = "HN" ;             # default
      if ($nsubst[$kk] > 3) {
        $typ[$i] = "HB" ;           # ammonium
        next ;
      }
      if ($typ[$kk] eq "N3" || $typ[$kk] eq "NC") {
        $typ[$i] = "HN" ;           # nonplanar amine 
        next ;
      }
      if ($typ[$kk] eq "N2" || $typ[$kk] eq "NP") {
        $typ[$i] = "HV" ;           # planar amine 
        next ;
      }
      next ;
    }
    if ($neibel eq "O") {           # connected to oxygen
                            # atom types HO,HX,HV
      $typ[$i] = "HO" ;             # default
      $l1 = $neib[$kk][1] ; # next neighbor to O
      $l2 = $neib[$kk][2] ;
      if ($l1 == $i) {
        $ll = $l2 ;
      }
      else {
        $ll = $l1 ;
      }
      if ($elem[$ll] ne "C") {next}
      if ($nsubst[$ll] == 4) {
        $typ[$i] = "HO" ;  # alcohol oxygen C-C-O-H
        next ;
      }
      if ($nsubst[$ll] != 3) {next}
      for ($j=1 ; $j <= $nsubst[$ll] ; $j++) {
        $mm = $neib[$ll][$j] ;
        if ($elem[$mm] eq "C" && $nbo[$ll][$mm] eq "d") {
          $typ[$i] = "HV" ; # enol oxygen  C=C-O-H 
          next ;
        }
        if ($elem[$mm] eq "O" && $nbo[$ll][$mm] eq "d") {
          $typ[$i] = "HX" ; # carboxylic oxygen  C(=O)OH
          next ;
        }
      }
      next ;
    }
  }
}

# testprint for hin format
if ($debug > 0) {
  print "\n" ;
  for ($j=1 ; $j <= $natoms ; $j++) {
    print "atom $j $elem[$j] $typ[$j] $nsubst[$j] " ;
#   print "atom $j - $elem[$j] $typ[$j] - $charge[$j] 0 0 0 $nsubst[$j] " ;
    if ($nsubst[$j] > 0) {
      for ($k=1 ; $k <= $nsubst[$j] ; $k++) {
        print "$neib[$j][$k] $bond[$j][$k] ";
      }
    }
    print "\n" ;
  }   
  print "\n" ;
}

# $natoms     # number of atoms
# $elem[i]    # element of atom i
# $typ[i]     # atom type of atom i "NA"
# $nsubst[i]  # number of 1st neighbors of atom i (including hydrogens) = number of bonds
# $neib[i][j] # list of 1st neighbors of atom i
# $tn12[i][j] # list of atom types of 1st neighbors, total of $nsubst[i] "NA-CA"
# $nangles[i] # number of angles involving atom i (as first atom)
# $tn13[i][j] # list of atom types of 2nd neighbors, total of $nangles[i] "NA-CA-CA"
# $ntors[i]   # number of torsions involving atom i (as first atom)
# $tn14[i][j] # list of atom types of 3rd neighbors, total of $ntors[i] "NA-CA-CA-BR"

# sequence j-ttj-ttk-ttl
#
for ($j=1 ; $j <= $natoms ; $j++) {
  $nangles[$j] = 0 ;
  $ntors[$j]   = 0 ;
  for ($k=1 ; $k <= $nsubst[$j] ; $k++) {
    $ttj  = $neib[$j][$k] ; # first neighbor
    $ttmp = $typ[$ttj] ; # atom type of first neighbor 
    $tn12[$j][$k] = $typ[$j]."-".$ttmp ; 
#   print "atom $j $typ[$j]  neighbor $ttj $typ[$ttj] tn12 $tn12[$j][$k]\n" ;
##  if ($nsubst[$j] < 2) {next} # only one neighbor
    for ($l=1 ; $l <= $nsubst[$ttj] ; $l++) { # loop over second neighbor
      $ttk = $neib[$ttj][$l] ; # neighbor of second neighbor
      if ($ttk == $j) {next} # k=i path to first atom, skip
      $nangles[$j]++ ;
      $tump = $typ[$ttk] ; # atom type of second neighbor
      $tn13[$j][$nangles[$j]] = $typ[$j]."-".$ttmp."-".$tump ; 
#     print "$j $ttj $ttk $tn13[$j][$nangles[$j]]\n" ;
      for ($m=1 ; $m <= $nsubst[$ttk] ; $m++) { # loop over third neighbor
        $ttl = $neib[$ttk][$m] ; # neighbor of second neighbor
#       if ($ttl == $j) {next}   # l=i path to first atom, skip
#       if ($ttl == $j) {print "l is i\n"}  # seemingly only in 4-membered rings
        if ($ttl == $ttj) {next} # l=j, skip
        $ntors[$j]++ ;
        $tvmp = $typ[$ttl] ; # atom type of second neighbor
        $tn14[$j][$ntors[$j]] = $typ[$j]."-".$ttmp."-".$tump."-".$tvmp ; 
#       print "$j $ttj $ttk $ttl $tn14[$j][$ntors[$j]]\n" ;
      }
    }
  }
}


# neighbor lists 
if ($debug > 0) {
  print "\n" ;
  for ($j=1 ; $j <= $natoms ; $j++) {
    for ($k=1 ; $k <= $nsubst[$j] ; $k++) {
      print "atom $j $tn12[$j][$k]\n" ;
    }
  }
  print "\n" ;
  for ($j=1 ; $j <= $natoms ; $j++) {
    if ($nangles[$j] < 1) {next}
    for ($l=1 ; $l <= $nangles[$j] ; $l++) {
      print "atom $j $tn13[$j][$l]\n" ;
    }
  }
  print "\n" ;
  for ($j=1 ; $j <= $natoms ; $j++) {
    if ($ntors[$j] < 1) {next}
    for ($m=1 ; $m <= $ntors[$j] ; $m++) {
      print "atom $j $tn14[$j][$m]\n" ;
    }
  }
}


} # end of subroutine smi2hin


sub proparray{
#
# setup propability arrays 
# $qf11[i]    # propability of having atom i at this position, total of $natoms
# $qf12[i][j] # propability of having atom j as 1st neighbor, total of $nsubst[i]
# $tn12[i][j] # list of atom types of 1st neighbors, total of $nsubst[i] "NA-CA"
# $qf13[i][j] # propability of having i-j-k as neighbors, total of $nangles[i] 
# $tn13[i][j] # list of atom types of 2nd neighbors, total of $nangles[i] "NA-CA-CA"
# $qf14[i][j] # propability of having i-j-k-l as neighbors, total of $ntors[i] 
# $tn14[i][j] # list of atom types of 3rd neighbors, total of $ntors[i] "NA-CA-CA-BR"
	
for ($j=1 ; $j <= $natoms ; $j++) {
  $ttt = $taty{$typ[$j]} ; # convert string to numerical atom type
  $qf11[$j] = $p11[$ttt]  ; # propability of this atom type
  for ($k=1 ; $k <= $nsubst[$j] ; $k++) {
    @cls = split(/-/,$tn12[$j][$k]); # split at dash sign
    $ttt = $taty{$cls[0]} ; # convert string to numerical atom type
    $uuu = $taty{$cls[1]} ; 
    if ($uuu >= $ttt) {
      $dij = ($ttt-1)*$maxt + $uuu ; # j>=i
    }
    else {
      $dij = ($uuu-1)*$maxt + $ttt ; # i>j
    }
    $qf12[$j][$k] = $p12[$dij]  ; # propability of this bond
#   print "j $j  k $k  cls0 $cls[0] ttt $ttt  cls1 $cls[1] uuu $uuu  dij $dij \n" ;
#   print "j $j k $k $ttt $uuu $cls[0] $cls[1] $p12[$dij]\n" ;
  }
  for ($l=1 ; $l <= $nangles[$j] ; $l++) {
    @cls = split(/-/,$tn13[$j][$l]); # split at dash sign
    $ttt = $taty{$cls[0]} ; # convert string to numerical atom type
    $uuu = $taty{$cls[1]} ;  
    $vvv = $taty{$cls[2]} ;  
    if ($vvv >= $ttt) {
      $dij = ($ttt-1)*$maxt + ($uuu-1)*$maxt2 + $vvv ; # k>=i
    }
    else {
      $dij = ($vvv-1)*$maxt + ($uuu-1)*$maxt2 + $ttt ; # i>k
    }
    $qf13[$j][$l] = $p13[$dij]  ; # propability of this angle
#   print "j $j k $k l $l cls0 $cls[0] ttt $ttt cls1 $cls[1] uuu $uuu cls2 $cls[2] vvv $vvv dij $dij \n" ;
#   print "j $j k $k l $l $ttt $uuu $vvv $cls[0] $cls[1] $cls[2] $p13[$dij]\n" ;
  }
  for ($m=1 ; $m <= $ntors[$j] ; $m++) {
    @cls = split(/-/,$tn14[$j][$m]); # split at dash sign
    $ttt = $taty{$cls[0]} ; # convert string to numerical atom type
    $uuu = $taty{$cls[1]} ; 
    $vvv = $taty{$cls[2]} ; 
    $www = $taty{$cls[3]} ; 
    if ($www >= $ttt) {
      $dij = ($ttt-1)*$maxt + ($uuu-1)*$maxt2 + ($vvv-1)*$maxt3 + $www ; # l>=i, ijkl
    }
    else {
      $dij = ($www-1)*$maxt + ($vvv-1)*$maxt2 + ($uuu-1)*$maxt3 + $ttt ; # i>l, lkji
    }
    $qf14[$j][$m] = $p14[$dij]  ; # propability of this torsion
#   print "j $j k $k l $l m $m $ttt $uuu $vvv $www $cls[0] $cls[1] $cls[2] $cls[3] $p14[$dij]\n" ;
  }
}

} # end of subroutine proparray
