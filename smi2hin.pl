#!/usr/bin/perl -w 

# convert one or multiple SMILES from a .smi file to multiple .hin files
#
# usage: smi2hin.pl [-options] infile.smi
#
# options: -mm  produce MM2 force field consistent atom types
#               instead of Hyperchem compatible ones (default)
#
# limitations: aromaticity defined by lower case letters does not work (yet)!
#
# SMILES are prepared as follows:
# dative nitro groups to double bonds   [N+]([O-])=O  to N(=O)=O
# add chloride to remaining charged nitrogens  [N+] to N(Cl)
# add chloride to remaining charged phosphorus [P+] to P(Cl)


$option = $ARGV[0] ;     # the first argument of the command line
$INFILE = $ARGV[1] ;     # the second argument of the command line

# check for additional command line arguments 
# to override following default options:
$dative = 1 ;
$neutra = 1 ;  # neutralize acidic oxygens
$aromat = 1 ;  # force aromatic bond information in .hin file
$hyperc = 1 ;  # HYPERCHEM compatible atom type assignment
               # -NO2 groups:  N1(=O1)(=O1) instead of NO(=O1)(=O1)
               # -N3 groups:   NA=N1=NA     instead of NA=NZ=N1 

if (substr($option,0,1) ne "-") {
  $INFILE = $option ;
}
else { # evaluate options
  if ($option eq "-mm") {$hyperc = 0}
}

open(INPU,$INFILE) or die "error opening $INFILE : $!\n"; 

$hinname = " " ;
$dcount = 0 ;     # counter for dummy filenames in case no CID is specified
$tname = "t01234.smi" ;  # temporary file
$hname = "t01234.hin" ;  # temporary file
$natoms = 0           ;  # total number of atoms
$nheavy = 0           ;  # number of non-hydrogen atoms
$nhydr = 0            ;  # number of hydrogen atoms

while(defined($i = <INPU>)) {
  @cols = split(/[\s\t]+/,$i); # split at spaces or tabs

  $smi = " "  ;     # compound SMILES
  $cid = "" ;       # compound CID   will be file name of the .hin file
## print "$cols[0] $cols[1] \n";

  if ($cols[0] ne "") {
    $smi = $cols[0] ;
  }
  if ($cols[1] ne "") {
    $cid = $cols[1] ;
  }
  if (defined($cols[2])) {
    if ($cols[2] ne "") {      # splitted name, e.g. acetic ester
      if ($cols[2] =~ /[A-Z]/) {} # e.g. CHIRAL NATURAL WITHDRAWN
      else {
        $cid = $cid.$cols[2] ;
      }
    }
  }

# print "$cid\n" ; # debug only

# identifier or name of compound
# create uniqe dummy name if no name is given
  if ($cid eq "") {
    $cid = "d".$dcount ;
    $dcount++ ;
  }
  $hinname = $cid.".hin" ;

# print "$hinname " ;

  if ($smi ne "") { 

# substitute [NH] to N
    $smi =~ s/\[NH\]/N/g ;

# substitute some common isotopes
    $smi =~ s/\[2H\]/H/g ;
    $smi =~ s/\[13C\]/C/g ;
    $smi =~ s/\[11C\]/C/g ;
    $smi =~ s/\[10B\]/B/g ;
    $smi =~ s/\[18F\]/F/g ;
    $smi =~ s/\[131I\]/I/g ;

#   print "$smi\n" ; # debug only

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

#print "\n          1         2         3         4         5\n" ;
#print "012345678901234567890123456789012345678901234567890\n" ;
#print "$smi\n\n" ;
#print "Atom\tPos\tSubst\tNeighborlist\n" ;
#for ($i = 0; $i < $pos ; $i++) {
#  if ($el[$i] =~ /[A-Z]/) {
#    print "$el[$i]\t$i\t$nsub[$i]\t$neighbors[$i]\n" ;
#  }
#}
#
#print "\nAtom\tCharge\tBondlist\n" ;
#for ($i = 0; $i < $pos ; $i++) {
#  if ($el[$i] =~ /[A-Z]/) {
#    print "$el[$i]\t$chrg[$i]\t$bondlist[$i]\n" ;
#  }
#}

# now derive atoms sequentially 
# Hier: Verknuepfungsmatrix $neib[][] aufstellen, Wasserstoffe hinzufuegen
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

# testprint for hinfile
#  print "\n" ;
#  for ($j=1 ; $j <= $natoms ; $j++) {
#    print "atom $j - $elem[$j] $typ[$j] - $charge[$j] 0 0 0 $nsubst[$j] " ;
#    if ($nsubst[$j] > 0) {
#      for ($k=1 ; $k <= $nsubst[$j] ; $k++) {
#        print "$neib[$j][$k] $bond[$j][$k] ";
#      }
#    }
#    print "\n" ;
#  }   
#  print "\n" ;


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
    if ($nsubst[$i] == 4) {
      $val = 4 ;
      $charge[$i] = 1 ;
      $typ[$i] = "NH" ; # ammonium
    }
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
          if ($hyperc == 1) {
            $typ[$i] = "N1" ;   # Hyperchem assigns N1 to nitro nitrogens
          }
          else {
            $typ[$i] = "NO" ;   # MM3 assigns NO
          }
          for ($k=1 ; $k <= $oneib ; $k++) {
            $typ[$onei[$k]] = "O1" ;
#   print "atom $onei[$k] $charge[$onei[$k]] $nsubst[$onei[$k]]\n" ;
            if ($nsubst[$onei[$k]] == 1) {
#             if ($bond[$onei[$k]][1] eq "s") {
#               $bond[$onei[$k]][1] = "d" ; #  force N=O
#               $charge[$onei[$k]] = 0 ;
#             }
            }
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
##print "\n" ;

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


# Assign MM+ force field atom types (similar to MM2 and MM3)
#

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
#     $typ[$i] = "C?" ;
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
#         if ($elem[$kk] eq "C" && $nsubst[$kk] == 2) {  corrected 10.01.14
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
#     $typ[$i] = "O?" ;
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
#     $typ[$i] = "S?" ;
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
      if ($hyperc == 1) {
        if ($bond[$i][1] eq "d" && $elem[$neib[$i][1]] eq "N") {
          $typ[$i] = "NA" ; # -N=N=NA assignment of Hyperchem
        }
      }
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
            if ($hyperc == 1) {
              $typ[$i] = "N1" ;  # Hyperchem assigns N1
            }
            else {
              $typ[$i] = "NZ" ;  # MM3 assigns NZ
            }
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
#     $typ[$i] = "N?" ;
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

# Vielleicht hier schon speudo 3D anhand der Atombindungen
# for visual inspection put everything into the xy-plane
# all bond distances are set to 1.4
# bond angles are set to 120 degrees, except small rings
#
#$cpn = 0 ; # angle of compass needle with y-axis
#for ($j=1 ; $j <= $natoms ; $j++) {
#  $x[$j] = 0 ;
#  $y[$j] = 0 ;
#  $z[$j] = 0 ;
#  if ($ring6[$j] == 1) {$cpn = $cpn + 30} 
#  for ($k=1 ; $k <= $nsubst[$j] ; $k++) {
#    $tt = $neib[$j][$k] ;
#    if ($elem[$tt] ne "H") {
#    }
#    else {
#    }
#  }
#}


### testprint for hinfile
#  print "\n" ;
#  for ($j=1 ; $j <= $natoms ; $j++) {
#    print "atom $j - $elem[$j] $typ[$j] - $charge[$j] 0 0 0 $nsubst[$j] " ;
##   print "atom $j - $elem[$j] $typ[$j] - $charge[$j] $x[$j] $y[$j] $z[$j] $nsubst[$j] " ;
#    if ($nsubst[$j] > 0) {
#      for ($k=1 ; $k <= $nsubst[$j] ; $k++) {
#        print "$neib[$j][$k] $bond[$j][$k] ";
#      }
#    }
#    print "\n" ;
#  }   
#  print "\n" ;

# azido groups N=[N+]=[N-]  to N=N=N
# isonitrile groups [N+]#[C-] to N#C
# nitroxides [N+]([O-]) to =N(=O)
# 

# write .hin file

    open(THIN,"> $hname") or die "error opening $hname : $!\n"; 

      print THIN "; $cid\n" ;
      print THIN "; $smi\n" ;
      print THIN "forcefield mm+\n" ;
      print THIN "mol 1\n" ;

      for ($j=1 ; $j <= $natoms ; $j++) {
        print THIN "atom $j - $elem[$j] $typ[$j] - $charge[$j] 0 0 0 $nsubst[$j] " ;
#       print THIN "atom $j - $elem[$j] $typ[$j] - $charge[$j] $x[$j] $y[$j] $z[$j] $nsubst[$j] " ;
        if ($nsubst[$j] > 0) {
          for ($k=1 ; $k <= $nsubst[$j] ; $k++) {
            print THIN "$neib[$j][$k] $bond[$j][$k] ";
          }
        }
        print THIN "\n" ;
      }   

      print THIN "endmol 1\n" ;

    close(THIN) or die "error closing $hname : $!\n";

#   rename temporary .hin file to final name

    rename($hname,$hinname) ;

  }
}

# clean up, remove temporary files

unlink($tname) ;

close (INPU) or die "error closing $INFILE : $!\n";
exit(0) ;
