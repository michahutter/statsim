#!/usr/bin/perl -w 

# Compute the distributions/frequencies of atom pair interactions
#
# usage: distrib.pl [-options] list_of_hin_files.txt
#
# options: -<matrix_name>  override "atpairfrq" as matrix name 
#
# <matrix_name>.mtx  matrix files will be generated
$MNAME  = "atpairfrq"; 

$option = $ARGV[0] ;     # the first argument of the command line
$CLIST  = $ARGV[1] ;     # the second argument of the command line


if (substr($option,0,1) ne "-") {
  $CLIST = $option ;
}
else { # evaluate options
  $MNAME = substr($option,1) ; # override atpairfrq as matrix_name
}

if (-e "$CLIST") { print "list of .hin files $CLIST found\n" } 
  else { print "list of .hin files $CLIST not found\n" ; die }


open(INPU,$CLIST) or die "error opening $CLIST : $!\n"; 

$nhin = 0 ; # number of compounds

while(defined($i = <INPU>)) {
  $nhin++ ;
  $compnam[$nhin] = $i ;
}

print "number of .hin files: $nhin\n" ;

close (INPU) or die "error closing $CLIST : $!\n";

$ahin  = 0 ; # actual number of existing .hin files

$ntat   = 0 ; # sum of all atoms (atom types)
$nt12   = 0 ; # sum of all 1-2 relationships (bonds)
$nt123  = 0 ; # sum of all 1-2-3 relationsships (angles)
$nt1234 = 0 ; # sum of all 1-2-3-4 relationships (torsions)

# List of available atom types
# here: those from the MM+ force field plus defaults for each element
#
$t[0]   = "ZZ" ; # undefined atom type, not used

$t[1]   = "C4" ; # Csp3
$t[2]   = "C3" ; # Csp2 alkene
$t[3]   = "CO" ; # Csp2 carbonyl
$t[4]   = "C2" ; # Csp  alkyne and C=C=O
$t[5]   = "CP" ; # Csp3 cyclopropane
$t[6]   = "C+" ; # Carboniumion
$t[7]   = "CZ" ; # Csp2 cyclopropene
$t[8]   = "CA" ; # Csp2 aromatic
$t[9]   = "CB" ; # Csp3 cyclobutane
$t[10]  = "CC" ; # Csp2 cyclobutene
$t[11]  = "CD" ; # Csp2 cyclobutanone C=O
$t[12]  = "CE" ; # Csp2 cyclopropanone C=O

$t[13]  = "N3" ; # Nsp3 
$t[14]  = "N2" ; # Nsp2 amide
$t[15]  = "N1" ; # Nssp
$t[16]  = "NA" ; # N -N= azo, pyridine
$t[17]  = "NH" ; # Nsp3 ammonium 
$t[18]  = "NP" ; # Nsp2 pyrrole
$t[19]  = "NB" ; # N -N=N-O azoxy
$t[20]  = "NZ" ; # N azide central N
$t[21]  = "NO" ; # N NO2 nitro
$t[22]  = "NC" ; # N =N- imine, oxime

$t[23]  = "O2" ; # O C-O-H, C-O-C
$t[24]  = "O1" ; # O =O carbonyl
$t[25]  = "OF" ; # Osp2 furan 
$t[26]  = "OC" ; # O carboxylate 
$t[27]  = "OE" ; # O epoxy
$t[28]  = "ON" ; # O amine oxide 

$t[29]  = "B3" ; # B trigonal  
$t[30]  = "B4" ; # B tetragonal

$t[31]  = "SI" ; # Si silicon

$t[32]  = "P" ;  # P >P- phosphine
$t[33]  = "P5" ; # P (V)

$t[34]  = "S2" ; # S -S- sulfide
$t[35]  = "S+" ; # S+ >S+ Sulfonium
$t[36]  = "SO" ; # S >S=O sulfoxide
$t[37]  = "S4" ; # S >SO2 sulfone
$t[38]  = "SA" ; # Ssp2 thiophene

$t[39]  = "F"  ; # F 
$t[40]  = "CL" ; # Cl
$t[41]  = "BR" ; # Br
$t[42]  = "I"  ; # I

$t[43]  = "H"  ; # hydrogen attached to C,B,Si,P,S
$t[44]  = "HO" ; # alcoholic
$t[45]  = "HN" ; # nonplanar amide
$t[46]  = "HV" ; # planar amide and enol oxygen C=C-OH
$t[47]  = "HB" ; # ammonium, tetravalent nitrogen 
$t[48]  = "HX" ; # carboxyl C(=O)OH 

$maxt = 48 ; # maximum number of atom types (index 0 is used)

$maxth = 38 ; # number without hydrogens and halogens as they cannot 
#               occur as middle atoms in angles and torsions

$maxt2 = $maxt * $maxt ;
$maxt3 = $maxt2 * $maxt ;

# Corresponding valences of the atom types
#

$val[0]  = 1 ; # not used
$val[1]  = 4 ; # Csp3 
$val[2]  = 3 ; # Csp2
$val[3]  = 3 ; # Csp2
$val[4]  = 2 ; # Csp  
$val[5]  = 4 ; # Csp3
$val[6]  = 3 ; # carbo cation
$val[7]  = 3 ; # Csp2
$val[8]  = 4 ; # Csp3
$val[9]  = 3 ; # Csp2
$val[10] = 3 ; # Csp2
$val[11] = 3 ; # Csp2
$val[12] = 3 ; # Csp2
$val[13] = 3 ; # Nsp3 N3
$val[14] = 2 ; # Nsp2 N2
$val[15] = 2 ; # Nsp2 N1
$val[16] = 2 ; # Nsp2 NA
$val[17] = 4 ; # Nsp3 ammonium NH
$val[18] = 3 ; # Nsp2 pyrrole NP
$val[19] = 2 ; # N -N=N-O NB azoxy middle N
$val[20] = 2 ; # N  NZ azide central N
$val[21] = 3 ; # N  NO NO2 nitro
$val[22] = 2 ; # N  NC =N- imine, oxime
$val[23] = 2 ; # O2 
$val[24] = 1 ; # O1 carbonyl oxygen 
$val[25] = 2 ; # OF Osp2 furan 
$val[26] = 1 ; # OC carboxylate 
$val[27] = 2 ; # OE epoxy 
$val[28] = 1 ; # ON amide oxide 
$val[29] = 3 ; # B3 trigonal 
$val[30] = 4 ; # B4 tetragonal 
$val[31] = 4 ; # SI silicon 
$val[32] = 3 ; # P >P- phosphine 
$val[33] = 4 ; # P5 (V)
$val[34] = 2 ; # S2 sulfide
$val[35] = 3 ; # S+ >S+ sulfonium  
$val[36] = 3 ; # SO >S=O sulfoxide
$val[37] = 4 ; # S4 >SO2 sulfon 
$val[38] = 2 ; # SA Ssp2 thiophene 
$val[39] = 1 ; # F 
$val[40] = 1 ; # Cl
$val[41] = 1 ; # Br
$val[42] = 1 ; # I
$val[43] = 1 ; # H attached to C,B,Si,P,S 
$val[44] = 1 ; # HO alcoholic
$val[45] = 1 ; # HN nonplanar amide
$val[46] = 1 ; # HV planar amide and enol oxygen C=C-OH
$val[47] = 1 ; # HB ammonium, tetravalent nitrogen
$val[48] = 1 ; # HX carboxyl C(=O)OH


# check if designated matrix files already exist 
#
for ($i=0 ; $i <= 3 ; $i++) {
  $ssum[$i] = 0 ;
  $tnam = $MNAME.$i.'.txt' ;
# print "$tnam \n" ;
  if (-e "$tnam") {
    print "matrix $tnam already exists ! aborting\n" ;
    die ;
  }
}

# initialize occurence lists and matrices
#

for ($i=1; $i <= $maxt; $i++) {
  $ocp[$i] = 0 ;        # occurrency of atom type p
  for ($j=$i; $j <= $maxt; $j++) { # j>=i
    $dij = ($i-1) * $maxt + $j ; # index ij
    $oc12[$dij] = 0 ; # occurrency of 1-2 relationship (bond)
    $s12[$dij]  = 0 ; # similarity matrix/list
#   print "$i $j  $dij  $t[$i] $t[$j]\n" ; # debug only
  }
}

for ($i=1; $i <= $maxt; $i++) {
  for ($j=1; $j <= $maxth; $j++) { 
    for ($k=$i; $k <= $maxt; $k++) { # k>=i
      $dijk = ($i-1) * $maxt + ($j-1) * $maxt2 + $k ; # index ijk
      $oc123[$dijk] = 0 ; # occurrency of 1-2-3 relationship (angle)
      $s123[$dijk]  = 0 ;   
#     print "$i $j $k  $dijk  $t[$i] $t[$j] $t[$k]\n" ; # debug only
    }
  }
}

for ($i=1; $i <= $maxt; $i++) { # ijkl = lkji, l>=i
  for ($j=1; $j <= $maxth; $j++) { 
    for ($k=1; $k <= $maxth; $k++) {
      for ($l=$i; $l <= $maxt; $l++) { 
        $dijkl = ($i-1) * $maxt + ($j-1) * $maxt2 + ($k-1) * $maxt3 + $l ; # index ijkl
        $oc1234[$dijkl] = 0 ; # occurrency of 1-2-3-4 relationship (torsion)
        $s1234[$dijkl] = 0 ;  
#       print "$i $j $k $l  $dijkl $t[$i] $t[$j] $t[$k] $t[$l]\n" ; # debug only
      }
    }
  }
}


for ($i=1; $i<=$nhin; $i++) { # loop over .hin files

  $name = $compnam[$i] ;
  chomp($name) ;

  if (-e "$name") {
    $ahin++ ;
  } 
  else {
    print "file $name not found, skipping\n" ;
  }

  open(HIN,$name) or die "error opening $name : $!\n"; 

  $nat = 0 ; # number of atoms

  while (defined($jj = <HIN>)) {
    @cols = split(/[\s\t]+/,$jj); # split at spaces or tabs
    if (defined($cols[0])) {
      $f1 = $cols[0] ;
      $c1 = substr($cols[0],0,1) ; # first char of first field 
    }
    else { # empty line
      next ;
    }
    if ($c1 eq ";") {next}    # skip comment line
    if ($f1 ne "atom") {next} # skip anything not declaring an atom

#    atom 1 - C C4 - 0 0 0 0 4 2 s 21 s 22 s 23 s    

    $ntat++ ; # count of all atoms over all .hin files
    $nat++ ;  # count of all atoms in this .hin file
    $el = $cols[3] ;  # element
    $t0 = $cols[4] ;  # assigned atom type

    if ($t0 eq "**") { # unassigned atom types
      if ($el eq "H")  { $t0 = "H" }
      if ($el eq "C")  { $t0 = "C4" }
      if ($el eq "N")  { $t0 = "N3" }
      if ($el eq "O")  { $t0 = "O2" }
      if ($el eq "P")  { $t0 = "P5" }
      if ($el eq "S")  { $t0 = "S2" }
      if ($el eq "B")  { $t0 = "B3" }
      if ($el eq "F")  { $t0 = "F" }
      if ($el eq "Cl") { $t0 = "CL" }
      if ($el eq "Br") { $t0 = "BR" }
      if ($el eq "I")  { $t0 = "I" }
      if ($el eq "Si") { $t0 = "SI" }
    }
    $tt[$nat] = 0 ; 
  
    for ($k=0; $k<=$maxt; $k++) {
      if ($t0 eq $t[$k]) {
        $ocp[$k]++ ; # count of atom types 
        $tt[$nat] = $k ;
      }
    }

    $nn[$nat] = $cols[10] ; # number of neighbours of atom $nat

    if ($nn[$nat] == 0) {
      $nb[$nat][1] = 0 ;
    }
    else {
      for ($k=1 ; $k <= $nn[$nat]; $k++) {
        $nb[$nat][$k] = $cols[9+2*$k] ; # neighbour matrix of atom $nat
      }
   }
 } 

 close (HIN) or die "error closing $name : $!\n";


# determine occupancies and frequencies

for ($j=1 ; $j <= $nat; $j++) {
  if ($tt[$j] != 0) {   
    if ($nn[$j] > 0) {
      for ($kk=1 ; $kk <= $nn[$j] ; $kk++) {
        $q1 = $nb[$j][$kk] ; # first neighbour of atom i (bond)
        if ($tt[$q1] != 0 && $q1 != $j) {
          if ($tt[$j] <= $tt[$q1]) {
            $dij = ($tt[$j]-1) * $maxt + $tt[$q1] ; # j>=i
          }
	  else {
            $dij = ($tt[$q1]-1) * $maxt + $tt[$j] ; # i>j
          }
          $oc12[$dij]++ ; # ij != ji
          $nt12++ ;
	} 
	if ($nn[$q1] != 0 && $q1 != $j) {
	  for ($ll=1 ; $ll <= $nn[$q1] ; $ll++) {
	    $q2 = $nb[$q1][$ll] ; # second neighbour of atom i (angle)
	    if ($tt[$q2] != 0 && $q2 != $j && $q2 != $q1) {
              if ($tt[$j] <= $tt[$q2]) { 
                $dijk = ($tt[$j]-1)*$maxt + ($tt[$q1]-1)*$maxt2 + $tt[$q2] ; # k>=i
              }
	      else {
                $dijk = ($tt[$q2]-1)*$maxt + ($tt[$q1]-1)*$maxt2 + $tt[$j] ; # i>k
              }
              $oc123[$dijk]++ ;
              $nt123++ ;
	    }
	    if ($nn[$q2] != 0 && $q2 != $j) {
	      for ($mm=1 ; $mm <= $nn[$q2] ; $mm++) {
	        $q3 = $nb[$q2][$mm] ; # third neighbour of atom i (torsion)
	        if ($tt[$q3] != 0 && $q3 != $q1) {
                  if ($tt[$j] <= $tt[$q3]) { 
                    $dijkl = ($tt[$j]-1)*$maxt + ($tt[$q1]-1)*$maxt2 + ($tt[$q2]-1)*$maxt3 + $tt[$q3] ; # l>=i ijkl
                  }
                  else {
                    $dijkl = ($tt[$q3]-1)*$maxt + ($tt[$q2]-1)*$maxt2 + ($tt[$q1]-1)*$maxt3 + $tt[$j] ; # i>l lkji
                  }
	          $oc1234[$dijkl]++ ; 
	          $nt1234++ ;
                }
              }		 
            }
          }
        }
      }
    }
  }
}

} # end of loop over all .hin files

print "\n" ;
print "number of compounds $nhin\n" ;
print "number of atoms     $ntat\n" ;
print "number of bonds     $nt12\n" ;
print "number of angles    $nt123\n" ;
print "number of torsions  $nt1234\n" ;
print "\n" ;

# avoid any division or by zero or log(0) later on
if ($nt12 < 1) {$nt12 = 1} # no bonds at all
if ($nt123 < 1) {$nt123 = 1} # no angles  at all
if ($nt1234 < 1) {$nt1234 = 1} # no torsions at all

$trhold     = 1e-6 ; 
$trhold12   = 1e-12; 
$trhold123  = 1e-15 ; 
$trhold1234 = 1e-18 ; 

# determine probabilities p_i and value for null hypothesises

for ($j=1 ; $j<=$maxt; $j++) {
  if ($ocp[$j] < 1) {$ocp[$j] = 1} # avoid multiplication by zero for non-occuring types
  $soc12[$j]   = ($ocp[$j] * $val[$j]) / $nt12 ; 
  $soc123[$j]  = ($ocp[$j] * $val[$j]) / $nt123 ; 
  $soc1234[$j] = ($ocp[$j] * $val[$j]) / $nt1234 ; 
  if ($ocp[$j] < $trhold) {$ocp[$j] = $trhold } # atom type does not occur
  $ocp[$j] = $ocp[$j] / $ntat ; # p_i : probability of atom type i
}


# compute propabilities for matrix entries S_ij = ln(q_ij/(p_i * p_j)
#
# since not all combination of atom type i to atom type j are possible,
# and furthermore i and j are involved in multiple pairwise interactions
# the Null hypothesis p_i * p_j *... is computed as
#
#        product of (count of atomtype i * valence of atomtype i * number of involved atoms)
# p_ij = ------------------------------------------------------------------------------------
#              doubles counts * count of interactions ^ (number of involved atoms)
#
# doubles = 2 for bonds (ji = ij), 1.6875 for angles (kji = ijk), 1.333 for torsions (lkji = ijkl) 
#

for ($i=1; $i <= $maxt; $i++) { 
  for ($j=$i; $j <= $maxt; $j++) { # j>=i
    $dij = ($i-1) * $maxt + $j ; # index ij
    if ($oc12[$dij] > 0) {
      $tmp = $soc12[$i] * $soc12[$j] * 2.0 ;  # p_ij
      if ($tmp < $trhold12) {$tmp = $trhold12} ;
      $s12[$dij] = $tmp ;
    }
  }
}

for ($i=1; $i <= $maxt; $i++) {
  for ($j=1; $j <= $maxth; $j++) { 
    for ($k=$i; $k <= $maxt; $k++) { # k>=i
      $dijk = ($i-1) * $maxt + ($j-1) * $maxt2 + $k ; # index ijk
      if ($oc123[$dijk] > 0) {
#       $tmp = $soc123[$i] * $soc123[$j] * $soc123[$k] * 18.0 ; # p_ijk
        $tmp = $soc123[$i] * $soc123[$j] * $soc123[$k] * 16.0 ; # p_ijk
        if ($tmp < $trhold123) {$tmp = $trhold123} ;
        $s123[$dijk] = $tmp ;
      }
    }
  }
}

for ($i=1; $i <= $maxt; $i++) {
  for ($j=1; $j <= $maxth; $j++) { 
    for ($k=1; $k <= $maxth; $k++) {
      for ($l=$i; $l <= $maxt; $l++) { # l>i
        $dijkl = ($i-1) * $maxt + ($j-1) * $maxt2 + ($k-1) * $maxt3 + $l ; # index ijkl
        if ($oc1234[$dijkl] > 0) {
#         $tmp = $soc1234[$i] * $soc1234[$j] * $soc1234[$k] * $soc1234[$l] * 256.0 ;
          $tmp = $soc1234[$i] * $soc1234[$j] * $soc1234[$k] * $soc1234[$l] * 192.0 ; # p_ijkl
          if ($tmp < $trhold1234) {$tmp = $trhold1234} ;
          $s1234[$dijkl] = $tmp ;
        }
      }
    }
  }
}

# check if sum of probabilities adds up to 1.0

print "sum of probabilities (should be close to 1.0)\n" ;

$ssum = 0 ;
for ($i=1 ; $i<=$maxt; $i++) {
  $ssum = $ssum + $ocp[$i] ; 
}
print "atoms    $ssum\n" ;

$ssum = 0 ;
for ($i=1; $i <= $maxt; $i++) {
  for ($j=$i; $j <= $maxt; $j++) { # j>=i
    $dij = ($i-1) * $maxt + $j ; # index ij
    if ($oc12[$dij] > 0) {
      $ssum = $ssum + $s12[$dij] ;
      if ($s12[$dij] < $trhold12) {$s12[$dij] = $trhold12}
      $s12[$dij] = log($oc12[$dij] / ($nt12 * $s12[$dij])) ;
    }
  }
}
print "bonds    $ssum\n" ;

$ssum = 0 ;
for ($i=1; $i <= $maxt; $i++) {
  for ($j=1; $j <= $maxth; $j++) { 
    for ($k=$i; $k <= $maxt; $k++) { # k>=i
      $dijk = ($i-1) * $maxt + ($j-1) * $maxt2 + $k ; # index ijk
      if ($oc123[$dijk] > 0) {
        $ssum = $ssum + $s123[$dijk] ;
        if ($s123[$dijk] < $trhold123) {$s123[$dijk] = $trhold123}
        $s123[$dijk] = log($oc123[$dijk] / ($nt123 * $s123[$dijk])) ;
      }
    }
  }
}
print "angles   $ssum\n" ;


$ssum = 0 ;
for ($i=1; $i <= $maxt; $i++) {
  for ($j=1; $j <= $maxth; $j++) { 
    for ($k=1; $k <= $maxth; $k++) {
      for ($l=$i; $l <= $maxt; $l++) { # l> i
        $dijkl = ($i-1) * $maxt + ($j-1) * $maxt2 + ($k-1) * $maxt3 + $l ; # index ijkl
        if ($oc1234[$dijkl] > 0) {
          $ssum = $ssum + $s1234[$dijkl] ;
          if ($s1234[$dijkl] < $trhold1234) {$s1234[$dijkl] = $trhold1234}
          $s1234[$dijkl] = log($oc1234[$dijkl] / ($nt1234 * $s1234[$dijkl])) ;
        }
      }
    }
  }
}
print "torsions $ssum\n" ;

# write matrix filesm = $MNAME.$k.'.mtx' ;

for ($n=1 ; $n <= 4 ; $n++) {
  $tnam = $MNAME.$n.'.mtx' ;
  open(DATA,"> $tnam") or die "error writing to $tnam : $!\n"; 

  if ($n == 1) {
    for ($i=1 ; $i<=$maxt; $i++) {
      $tmp = $ocp[$i] * $ntat ;
#     printf DATA ("%-2s%9.6f",$t[$i],$ocp[$i]) ;
      printf DATA ("%-2s%9.6f%7i",$t[$i],$ocp[$i],$tmp) ;
      print DATA "\n" ;
    }
  }

  if ($n == 2) {
    for ($i=1; $i <= $maxt; $i++) {
      for ($j=$i; $j <= $maxt; $j++) { # j>=i
        $dij = ($i-1) * $maxt + $j ; # index ij
	if ($oc12[$dij] > 0) {
#         printf DATA ("%-2s%1s%-2s%8.3f",$t[$i],"-",$t[$j],$s12[$dij]) ;
          printf DATA ("%-2s%1s%-2s%8.3f%7i",$t[$i],"-",$t[$j],$s12[$dij],$oc12[$dij]) ;
          print DATA "\n" ;
        }
      }
    }
  }

  if ($n == 3) {
    for ($i=1; $i <= $maxt; $i++) {
      for ($j=1; $j <= $maxth; $j++) { 
        for ($k=$i; $k <= $maxt; $k++) { # k>=i
          $dijk = ($i-1) * $maxt + ($j-1) * $maxt2 + $k ; # index ijk
          if ($oc123[$dijk] > 0) {
#           printf DATA ("%-2s%1s%-2s%1s%-2s%8.3f",$t[$i],"-",$t[$j],"-",$t[$k],$s123[$dijk]) ;
            printf DATA ("%-2s%1s%-2s%1s%-2s%8.3f%7i",$t[$i],"-",$t[$j],"-",$t[$k],$s123[$dijk],$oc123[$dijk]) ;
            print DATA "\n" ;
          }
        }
      }
    }
  }

  if ($n == 4) {
    for ($i=1; $i <= $maxt; $i++) {
      for ($j=1; $j <= $maxth; $j++) { 
        for ($k=1; $k <= $maxth; $k++) {
          for ($l=$i; $l <= $maxt; $l++) { # l> i
            $dijkl = ($i-1) * $maxt + ($j-1) * $maxt2 + ($k-1) * $maxt3 + $l ; # index ijkl
            if ($oc1234[$dijkl] > 0) {
#             printf DATA ("%-2s%1s%-2s%1s%-2s%1s%-2s%8.3f",$t[$i],"-",$t[$j],"-",$t[$k],"-",$t[$l],$s1234[$dijkl]) ;
              printf DATA ("%-2s%1s%-2s%1s%-2s%1s%-2s%8.3f%7i",$t[$i],"-",$t[$j],"-",$t[$k],"-",$t[$l],$s1234[$dijkl],$oc1234[$dijkl]) ;
              print DATA "\n" ;
            }
          }
        }
      }
    }
  }

  close(DATA) ;
}

print "\nsuccessful completion\n" ;

exit(0) ;
