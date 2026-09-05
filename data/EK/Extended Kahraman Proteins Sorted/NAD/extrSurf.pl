#!/usr/bin/perl -w

#########################################################################
# extrSurf.pl extract surface residues from PDB protein structures	#
# There are several ways to extract surface residues. The following are	#
# implemented in this program.                                 		#
# 1. Obtain the pocket information for the structure and extract the	#
# pocket residues. CastP server is used to calculate pocket residues	#
# 2. If the structure contains a ligand, the surface residues can be	#
# obtained using the information of the ligand. The atoms/residues	#
# around the ligand will be extracted.                  		#
# 3. A file containing the surface accessible area of each atom similar	#
# as pdb file can be read as input and surface residues can be selected.#
# 4. Protein-protein interfaces can be extracted for all interface	#
# atoms.                                                         	#
#                                        				#
# 									#
# Author: Dr. Jinfeng Zhang                				#
# Contact: jinfeng@bioinfo.stat.harvard.edu				#
# Last updated: Jan 24, 2006						#
#########################################################################

# print out usage
$command = "extrSurf.pl";
if (@ARGV==0) {
    print "Usage: ",$command,"
Options: 
\t -e [p,l,s,s2,pdi,ppi,r] extraction type,default is p for pocket, 
l for ligand, 
s for surface residues defined as one atom > 5 A^2, 
s2 for surface residues defined as >25% ASA, 
pdi for protein-DNA interface
ppi for protein-protein interface.
r for adjacent residues for a given residu/residues
\t -f [] protein name, must provide.
\t -n [1-9] number of pocket files; for -e l, it is 1.
\t -ln [] name of ligand, for -e l, must provide.
\t -cn [] chain name of the proteins for pocket or contacting ligand.
\t -t [float] threshold for surface accessible area to be counted as exposed.
\t -lpd [] distance of extracted protein atoms to ligand, default 4.0.\n";

    exit;
}

# initialize parameters

$ExType="p";     # Extraction type, 
                 # "p": extract from pocket; 
                 # "l": extract from ligand;
                 # "s": extract surface residues from *.sf.pdb
                 # "s2": extract surface residues from *.sf.pdb based on a different criterion
                 # "i": extract interface residues from a complex pdb file
$pdbName="";     # pdb name
$numOut=1;       # number of pockets counted
$ligName="ATP";  # name the ligand
$chainName=" ";  # chain identifier
$lig_p_dis=6.0;  # the distance of extracted protein atoms to ligand
$threshold=5.0;  # threshold for surface accessible area to be counted as exposed
$CutOff=25;      # square of the distance cutoff (5.0)
$isResLev=0;     # whether the computation is done at residue level

# get command line arguments

$i=0;
while($i<@ARGV) {
    if($ARGV[$i] eq "-e"){
	$ExType = $ARGV[++$i];
    }
    elsif($ARGV[$i] eq "-f"){
	$pdbName = $ARGV[++$i];
    }
    elsif($ARGV[$i] eq "-n"){
	$numOut = $ARGV[++$i];
    }
    elsif($ARGV[$i] eq "-ln"){
	$ligName = $ARGV[++$i];
    }
    elsif($ARGV[$i] eq "-lpd"){
	$lig_p_dis = $ARGV[++$i];
    }
    elsif($ARGV[$i] eq "-cn"){
	$chainName = $ARGV[++$i];
    }
    elsif($ARGV[$i] eq "-chains") {  # for protein-protein interface
	$Chain1=$ARGV[++$i];
	$Chain2=$ARGV[++$i];
    }
    elsif($ARGV[$i] eq "-dis"){
	$CutOff = $ARGV[++$i];
    }
    elsif($ARGV[$i] eq "-isResLev"){
	$isResLev=1;
    }
    elsif($ARGV[$i] eq "-r"){
	$selectRes{$ARGV[++$i]} = 1;
    }
    $i++;
}

# amino acid and their integer type

$ResLab=0;
$Res{"A"}=$ResLab;$Res{$ResLab++}="A";$Res{"C"}=$ResLab;$Res{$ResLab++}="C";
$Res{"D"}=$ResLab;$Res{$ResLab++}="D";$Res{"E"}=$ResLab;$Res{$ResLab++}="E";
$Res{"F"}=$ResLab;$Res{$ResLab++}="F";$Res{"G"}=$ResLab;$Res{$ResLab++}="G";
$Res{"H"}=$ResLab;$Res{$ResLab++}="H";$Res{"I"}=$ResLab;$Res{$ResLab++}="I";
$Res{"K"}=$ResLab;$Res{$ResLab++}="K";$Res{"L"}=$ResLab;$Res{$ResLab++}="L";
$Res{"M"}=$ResLab;$Res{$ResLab++}="M";$Res{"N"}=$ResLab;$Res{$ResLab++}="N";
$Res{"P"}=$ResLab;$Res{$ResLab++}="P";$Res{"Q"}=$ResLab;$Res{$ResLab++}="Q";
$Res{"R"}=$ResLab;$Res{$ResLab++}="R";$Res{"S"}=$ResLab;$Res{$ResLab++}="S";
$Res{"T"}=$ResLab;$Res{$ResLab++}="T";$Res{"V"}=$ResLab;$Res{$ResLab++}="V";
$Res{"W"}=$ResLab;$Res{$ResLab++}="W";$Res{"Y"}=$ResLab;$Res{$ResLab++}="Y";

#$resSet="ACDEFGHIKLMNPQRSTVWY";

# map three letter aa name to one letter aa name

$one_thr{"ALA"}="A";$one_thr{"CYS"}="C";$one_thr{"ASP"}="D";$one_thr{"GLU"}="E";
$one_thr{"PHE"}="F";$one_thr{"GLY"}="G";$one_thr{"HIS"}="H";$one_thr{"ILE"}="I";
$one_thr{"LYS"}="K";$one_thr{"LEU"}="L";$one_thr{"MET"}="M";$one_thr{"ASN"}="N";
$one_thr{"PRO"}="P";$one_thr{"GLN"}="Q";$one_thr{"ARG"}="R";$one_thr{"SER"}="S";
$one_thr{"THR"}="T";$one_thr{"VAL"}="V";$one_thr{"TRP"}="W";$one_thr{"TYR"}="Y";

$one_thr{"A"}="ALA";$one_thr{"C"}="CYS";$one_thr{"D"}="ASP";$one_thr{"E"}="GLU";
$one_thr{"F"}="PHE";$one_thr{"G"}="GLY";$one_thr{"H"}="HIS";$one_thr{"I"}="ILE";
$one_thr{"K"}="LYS";$one_thr{"L"}="LEU";$one_thr{"M"}="MET";$one_thr{"N"}="ASN";
$one_thr{"P"}="PRO";$one_thr{"Q"}="GLN";$one_thr{"R"}="ARG";$one_thr{"S"}="SER";
$one_thr{"T"}="THR";$one_thr{"V"}="VAL";$one_thr{"W"}="TRP";$one_thr{"Y"}="TYR";

#$dna_sym{"  A"}="a";$dan_sym{"  C"}="c";$dna_sym{"  G"}="g";$dna_sym{"  T"}="t";

# accessible surface area of twenty amino acid when fully exposed
@ASA =   (115,135,150,190,
	  210, 75,195,175,
	  200,170,185,160,
	  145,180,225,115,
	  140,155,255,230);

# number of heavy atoms for each amino acid
@NUM_ATOM =   (5,6,8,9,
	       11,4,10,8,
	       9,8,8,8,
	       7,9,11,6,
	       7,7,14,12);

# get the residue level contact distance from file ContDis.txt
if($isResLev==1){
    if(-e "/cluster/home/jinfeng/SharedData/ContDis.txt"){
	open CONT, "/cluster/home/jinfeng/SharedData/ContDis.txt";
    }
    else {
	print "cannot open file ContDis.txt from /cluster/home/jinfeng/SharedData/ContDis.txt!\n";
	exit;
    }
    $j=0;
    while(defined($line=<CONT>)) {
	chomp($line);
	$line =~ s/^\s+//g;
	@contList=split(/\s+/,$line);
	for($i=0;$i<20;$i++) {
	    $contDis->[$j][$i]=$contList[$i];
	}
	++$j;
    }
}

print "protein name: $pdbName\n";
if($ExType eq "p"){ # extract from pocket
# get the pocket information from *.pocInfo
    $midName=substr($pdbName,1,2);
    if(-e "/cast/$midName/$pdbName.pocInfo"){
	`tail -n $numOut /cast/$midName/$pdbName.pocInfo > pocInfo.tmp`;
	open PocInfo, "pocInfo.tmp";
	$j=0; # number of pocket
# get the coordinates for each pocket from *.poc and output to *.$j.poc
	while(defined($line=<PocInfo>)){
	    @arr=split(/\s+/,$line);
	    $PocInd=$arr[2];
	    open Poc, "/cast/$midName/$pdbName.poc";
	    open OUT, ">$pdbName.$j.poc";
	    while(defined($pocLine=<Poc>)){
		if(length($pocLine) >=70){
		    $tmpNum=substr($pocLine,66,5);
		    $tmpNum=~ s/^\s+|\s+$//g;
		    if($tmpNum == $PocInd){
			printf OUT "$pocLine";
		    }
		}
	    }
	    close(Poc);
	    close(OUT);
	    $j++;
	}
	close(PocInfo);
    }
    else { 
	print "no pocket file found for $pdbName.\n";
	exit;
    }
}
elsif($ExType eq "l"){ # extract from ligand
# get the ligand atoms
    if(-e "$pdbName.pdb"){
	open PDB, "$pdbName.pdb";
    }
    elsif(-e "/home/jinfeng/FuncInf/PDB/$pdbName.pdb"){
	open PDB, "/home/jinfeng/FuncInf/PDB/$pdbName.pdb";
    }
    else {
	print "no pdb file.\n";
	exit;
    }
    $i=0;  # index of ligand atoms
    while(defined($line=<PDB>)){
	chomp($line);
    	$tempName=substr($line,17,3);
	$tempName=~ s/^\s+|\s+$//g;
	if(substr($line,0,6) eq "HETATM" && $tempName eq $ligName && (substr($line,21,1) eq $chainName || $chainName eq " ")){
	    $ligCo->[$i][0]=substr($line,30,8); 
	    $ligCo->[$i][0] =~ s/^\s+|\s+$//g;
	    $ligCo->[$i][1]=substr($line,38,8); 
	    $ligCo->[$i][1] =~ s/^\s+|\s+$//g;
	    $ligCo->[$i][2]=substr($line,46,8); 
	    $ligCo->[$i][2] =~ s/^\s+|\s+$//g;
	    $i++;
	    $chainName=substr($line,21,1);
	}
    }
    $ligAtmCnt=$i;   # number of ligand atoms
    close(PDB);
# get the surface atoms/residues, criteria dis<4A.
    if(-e "$pdbName.sf.pdb"){
	open Surf, "$pdbName.sf.pdb";
    }
    elsif(-e "$pdbName.pdb"){
	open Surf, "$pdbName.pdb";
    }
    else {
	open Surf, "/home/jinfeng/FuncInf/PDB/$pdbName.pdb";
    }
    open OUT, ">$pdbName.lig.$ligName.pdb";
    while(defined($line=<Surf>)){
	chomp($line);
	if(substr($line,0,4) eq "ATOM"){
	    $atmCo[0]=substr($line,30,8); 
	    $atmCo[0] =~ s/^\s+|\s+$//g;
	    $atmCo[1]=substr($line,38,8); 
	    $atmCo[1] =~ s/^\s+|\s+$//g;
	    $atmCo[2]=substr($line,46,8); 
	    $atmCo[2] =~ s/^\s+|\s+$//g;
	    for($i=0;$i<$ligAtmCnt;++$i){
#		print "$atmCo[0] $atmCo[1] $atmCo[2]\n";
		$dis=sqrt(($atmCo[0]-$ligCo->[$i][0])*($atmCo[0]-$ligCo->[$i][0])+($atmCo[1]-$ligCo->[$i][1])*($atmCo[1]-$ligCo->[$i][1])+($atmCo[2]-$ligCo->[$i][2])*($atmCo[2]-$ligCo->[$i][2]));
		if($dis <= $lig_p_dis){
		    printf OUT "%s %8.3f\n",$line, $dis;
		    last;
		}
	    }
	}
    }
    close(Surf);
    close(OUT);
}
elsif($ExType eq "s"){ # extraction type is surface over buried
    $selected=-100;
    $firstNum=-100;
    $lastNum=-100;
    $lastRes="AAA";
    if(-e "$pdbName.sf.pdb"){
	open Surf, "$pdbName.sf.pdb";
    }
    open OUT, ">$pdbName.sfr";
    while(defined($line=<Surf>)){
	if(substr($line,0,6) eq "ATOM  "){
	    $resNum  = substr($line, 22,4); $resNum  =~ s/^\s+|\s+$//g;
	    if($firstNum ==-100){
		$firstNum=$resNum;
	    }
	    if($selected==$resNum){
		next;
	    }
	    if($resNum != $lastNum && $lastNum != $selected){ # $lastNum is buried
		if(defined($one_thr{$resName})){
		    printf OUT "%d %s 0\n",$lastNum-$firstNum+1,$one_thr{$resName};
		}
		else {
		    printf OUT "%d U 0\n",$lastNum-$firstNum+1,$one_thr{$resName};
		}
	    }
	    $atomName= substr($line, 12,4); $atomName=~ s/^\s+|\s+$//g;
	    $resName = substr($line, 17,3); $resName =~ s/^\s+|\s+$//g;
# only residues with side chain atoms exposed on surface are counted 
	    if($atomName eq "N"||$atomName eq "CA"||$atomName eq "C"||$atomName eq "O"){
		$lastNum=$resNum;
		$lastRes=$resName;
		next;
	    }
	    $asa=substr($line,66,7); $asa=~ s/^\s+|\s+$//g;
	    if($asa > $threshold){  # surface residues are labeled as 1
		if(defined($one_thr{$resName})){
		    printf OUT "%d %s 1\n",$resNum-$firstNum+1,$one_thr{$resName};
		}
		else {
		    printf OUT "%d U 1\n",$resNum-$firstNum+1,$one_thr{$resName};
		}
		$selected=$resNum;
	    }
	    $lastNum=$resNum;
	    $lastRes=$resName;
	} # end of if
    } # end of while
    if($resNum != $selected){  
	if(defined($one_thr{$resName})){
	    printf OUT "%d %s 1\n",$resNum-$firstNum+1,$one_thr{$resName};
	}
	else {
	    printf OUT "%d U 0\n",$resNum-$firstNum+1,$one_thr{$resName};
	}
    }
}
elsif($ExType eq "s2"){  # definition of surface residue is that with >25% of relative accessible surface area
    $firstNum=-1000;
    $prevNum=-1000;
    $Threshold=0.25;    # threshold is 25%
    $area=0;
    $numAtom=0;
    if(-e "$pdbName.sf.pdb"){
	open Surf, "$pdbName.sf.pdb";
    }
    else {
	print "cannot open $pdbName.sf.pdb\n";
	exit(0);
    }
    open OUT, ">$pdbName.sfr";
    while(defined($line=<Surf>)){
	if(substr($line,0,6) eq "ATOM  "){
	    $resNum  = substr($line, 22,4); $resNum  =~ s/^\s+|\s+$//g;
	    if($firstNum ==-1000){
		$firstNum=$resNum;
	    }
	    if($prevNum != $resNum){
		if($prevNum!=-1000){  # assign for the last residue

		    $maxA=$ASA[$Res{$one_thr{$resName}}]*$numAtom/$NUM_ATOM[$Res{$one_thr{$resName}}];
		    #printf "%s %6.2f %6.1f %d %4.2f\n",$resName, $area, $maxA, $numAtom,$area/$maxA;
		    if(($area/$maxA)> $Threshold){  # surface residue labeled as 1
			if(defined($one_thr{$resName})){
			    printf OUT "%d %s 1\n",$resNum-$firstNum,$one_thr{$resName};
			}
			else {
			    printf OUT "%d U 1\n",$resNum-$firstNum,$one_thr{$resName};
			}
		    }
		    else {  # core residue labeled as 0
			if(defined($one_thr{$resName})){
			    printf OUT "%d %s 0\n",$resNum-$firstNum,$one_thr{$resName};
			}
			else {
			    printf OUT "%d U 0\n",$resNum-$firstNum,$one_thr{$resName};
			}
		    }
		}
		$resName = substr($line, 17,3); $resName =~ s/^\s+|\s+$//g;
		$area=0;
		$numAtom=1;
		$prevNum=$resNum;
		$asa=substr($line,66,7); $asa=~ s/^\s+|\s+$//g;
		$area+=$asa;
	    }
	    else {
		$asa=substr($line,66,7); $asa=~ s/^\s+|\s+$//g;
		$area+=$asa;
		$numAtom++;
	    }
	} # end of if
    } # end of while
# calculate last residue
    $maxA=$ASA[$Res{$one_thr{$resName}}];
    if(($area/$maxA)>$Threshold){  # surface residue labeled as 1
	if(defined($one_thr{$resName})){
	    printf OUT "%d %s 1\n",$resNum-$firstNum+1,$one_thr{$resName};
	}
	else {
	    printf OUT "%d U 1\n",$resNum-$firstNum+1,$one_thr{$resName};
	}
    }
    else {  # core residue labeled as 0
	if(defined($one_thr{$resName})){
	    printf OUT "%d %s 0\n",$resNum-$firstNum+1,$one_thr{$resName};
	}
	else {
	    printf OUT "%d U 0\n",$resNum-$firstNum+1,$one_thr{$resName};
	}
    }
}
elsif($ExType eq "pdi"){ # extracting type is protein-DNA interface
    if(-e "$pdbName.pdb"){
    	open PDB, "$pdbName.pdb";
    }
    
    
}
elsif($ExType eq "ppi"){ # extracting type is protein-protein interface
    if(-e "$pdbName.pdb"){
	open PDB, "$pdbName.pdb";
    }
    else {
	print "cannot open $pdbName.pdb!\n";
	exit();
    }

    open OUT1, ">$pdbName.intC"; # contact information of interface residues
    open OUT2, ">$pdbName.intR"; # interface residues

# get all atoms
    $lineCnt=0;
    while(defined($tmpline=<PDB>)){
	if(substr($tmpline,0,6) eq "ATOM  "){
	    if(substr($tmpline,16,4) =~ /[BCDE234][A-Z][A-Z][A-Z]/) {next;}
	    if(substr($tmpline,22,1) =~ /[B-Z]/) {next;}
	    $tmpCN=substr($tmpline,21,1);

	    if($isResLev==1){
		$resName=substr($tmpline,17,3);
		$atomName=substr($tmpline, 12,4); $atomName=~ s/^\s+|\s+$//g;
		if($resName eq "GLY" && $atomName ne "CA"){
		    next;
		}
		elsif($resName ne "GLY" && $atomName ne "CB"){
		    next;
		}
	    }
	    if($chainName =~ /$tmpCN/){
		$line[$lineCnt++]=$tmpline;
	    }
	}
    }

# calculate contact atoms at different chains

    for($i=0;$i<$lineCnt;++$i){
	$chainName1=substr($line[$i], 21,1);
	$resName1 = substr($line[$i], 17,3);
	$resNum1  = substr($line[$i], 22,4); $resNum1  =~ s/^\s+|\s+$//g;
	$atomName1= substr($line[$i], 12,4); $atomName1=~ s/^\s+|\s+$//g;
	$p1[0]=substr($line[$i],30,8); 
	$p1[0] =~ s/^\s+|\s+$//g;
	$p1[1]=substr($line[$i],38,8); 
	$p1[1] =~ s/^\s+|\s+$//g;
	$p1[2]=substr($line[$i],46,8); 
	$p1[2] =~ s/^\s+|\s+$//g;
	for($j=$i+1;$j<$lineCnt;++$j){
	    $chainName2=substr($line[$j], 21,1);
	    $resName2 = substr($line[$j], 17,3);
	    if($chainName1 eq $chainName2) {next;}
	    $p2[0]=substr($line[$j],30,8); 
	    $p2[0] =~ s/^\s+|\s+$//g;
	    $p2[1]=substr($line[$j],38,8); 
	    $p2[1] =~ s/^\s+|\s+$//g;
	    $p2[2]=substr($line[$j],46,8); 
	    $p2[2] =~ s/^\s+|\s+$//g;
	    $dis=sqrt(($p1[0]-$p2[0])*($p1[0]-$p2[0]) + ($p1[1]-$p2[1])*($p1[1]-$p2[1]) + ($p1[2]-$p2[2])*($p1[2]-$p2[2]));
	    if($isResLev == 1){
		#print "$chainName1 $chainName2 $resName1 $resName2 $Res{$one_thr{$resName1}} $Res{$one_thr{$resName2}} $contDis->[$Res{$one_thr{$resName1}}][$Res{$one_thr{$resName2}}] $dis\n";
		$CutOff=1.4*$contDis->[$Res{$one_thr{$resName1}}][$Res{$one_thr{$resName2}}];
	    }
	    if($dis<=$CutOff){
		$resNum2  = substr($line[$j], 22,4); $resNum2  =~ s/^\s+|\s+$//g;
		$atomName2= substr($line[$j], 12,4); $atomName2=~ s/^\s+|\s+$//g;
		if($isResLev==1){
		    printf OUT1 "$chainName1 %5d $resName1 $chainName2 %5d $resName2\n", $resNum1, $resNum2;
		    $resRecord{"$chainName1" . "$resNum1"} = 1;
		    $resRecord{"$chainName2" . "$resNum2"} = 1;
		}
		else {
		    printf OUT "$chainName1 %5d $resName1 %4s $chainName2 %5d $resName2 %4s\n",$resNum1,$atomName1,$resNum2,$atomName2;
} 
	    }
	}
    }
    if($isResLev==1){
	@outRecord = %resRecord;
	$j=0;
	for($i=0;$i<@outRecord;$i=$i+2){
	    $not_sorted[$j++]=$outRecord[$i];
	    #print "$outRecord[$i]\n";
	}
	@sorted = sort { lc($a) cmp lc($b) } @not_sorted;
	for($i=0;$i<@sorted;$i++){ 
	    printf OUT2 "$sorted[$i]\n";
	}
    }
}
elsif($ExType eq "r") { # extract adjacent residues of a given residue

    if(-e "$pdbName.pdb"){
	open PDB, "$pdbName.pdb";
    }
    elsif(-e "/home/jinfeng/SharedData/PDB/$pdbName.pdb"){
	open PDB, "/home/jinfeng/SharedData/PDB/$pdbName.pdb";
    }
    else {
	print "no pdb file.\n";
	exit;
    }
    open OUT, ">$pdbName.selRes";
    $lineCnt=0;  # index of selected residue atoms
    $numSelAtom=0; # number of atoms in selected residues
    while(defined($line=<PDB>)){
	if(substr($line,0,6) eq "ATOM  "){
	    if(substr($line, 16,4) =~ /[BCDE234][A-Z][A-Z][A-Z]/) {next;}
	    if(substr($line,22,1) =~ /[B-Z]/) {next;}
	    $line[$lineCnt++]=$line;
	    $resNum  = substr($line, 22,4); $resNum  =~ s/^\s+|\s+$//g;
	    $atomName= substr($line, 12,4); $atomName=~ s/^\s+|\s+$//g;
	    if(defined($selectRes{$resNum})&&$atomName eq "SG"){
		$selResCo->[$numSelAtom][0]=substr($line,30,8); 
		$selResCo->[$numSelAtom][0] =~ s/^\s+|\s+$//g;
		$selResCo->[$numSelAtom][1]=substr($line,38,8); 
		$selResCo->[$numSelAtom][1] =~ s/^\s+|\s+$//g;
		$selResCo->[$numSelAtom][2]=substr($line,46,8); 
		$selResCo->[$numSelAtom][2] =~ s/^\s+|\s+$//g;
		$numSelAtom++;
	    }
	}
    }

    close(PDB);

    for($i=0;$i<$lineCnt;++$i){
	$resNum  = substr($line[$i], 22,4); $resNum  =~ s/^\s+|\s+$//g;
	$resName = substr($line[$i], 17,3);
	$atomName= substr($line[$i], 12,4); $atomName=~ s/^\s+|\s+$//g;
	#if(defined($extrRes{$resNum})){next;}
	if(defined($selectRes{$resNum})){
	    $extrRes{$resNum} = $resName;
	    #printf OUT "$resNum $resName\n";
	}
	else {
	    #if($atomName eq "CA" || $atomName eq "C" || $atomName eq "O" || $atomName eq "N" || $atomName eq "CB"){next;}
	    $atmCo[0]=substr($line[$i],30,8); 
	    $atmCo[0] =~ s/^\s+|\s+$//g;
	    $atmCo[1]=substr($line[$i],38,8); 
	    $atmCo[1] =~ s/^\s+|\s+$//g;
	    $atmCo[2]=substr($line[$i],46,8); 
	    $atmCo[2] =~ s/^\s+|\s+$//g;
	    #print "atomco: $atmCo[0] $atmCo[1] $atmCo[2]\n";
	    for($j=0;$j<$numSelAtom;++$j){
		#print "$selResCo->[$j][0] $selResCo->[$j][1] $selResCo->[$j][2]\n"; 
		$dis=sqrt(($atmCo[0]-$selResCo->[$j][0])*($atmCo[0]-$selResCo->[$j][0])+($atmCo[1]-$selResCo->[$j][1])*($atmCo[1]-$selResCo->[$j][1])+($atmCo[2]-$selResCo->[$j][2])*($atmCo[2]-$selResCo->[$j][2]));
		if($dis <= $lig_p_dis){
		    $extrRes{$resNum} = $resName;
		    #printf OUT "$line";
		    printf OUT "$j $resNum $resName $atomName $dis\n"; 
		    #last;
		}
	    }
	}
    }
    close(OUT);
}
else {
    print "Extraction type given is wrong, choose from p, l, s, ppi, or pdi.\n";
}

exit;

# backup codes
for($i=0;$i<4;++$i){
    $chainName=substr($line,21,1);
    $resName = substr($line, 17,3);
    $resNum  = substr($line, 22,4); $resNum  =~ s/^\s+|\s+$//g;
    if($chainName eq $Chain1){
	$ch1Co->[$ch1Cnt][0]=substr($line,30,8); 
	$ch1Co->[$ch1Cnt][0] =~ s/^\s+|\s+$//g;
	$ch1Co->[$ch1Cnt][1]=substr($line,38,8); 
	$ch1Co->[$ch1Cnt][1] =~ s/^\s+|\s+$//g;
	$ch1Co->[$ch1Cnt][2]=substr($line,46,8); 
	$ch1Co->[$ch1Cnt][2] =~ s/^\s+|\s+$//g;
	$ch1Res[$ch1Cnt]=$resName;
	$ch1Num[$ch1Cnt]=$resNum;
	$ch1Cnt++;
    }
    if($chainName eq $Chain2){
	$ch2Co->[$ch2Cnt][0]=substr($line,30,8); 
	$ch2Co->[$ch2Cnt][0] =~ s/^\s+|\s+$//g;
	$ch2Co->[$ch2Cnt][1]=substr($line,38,8); 
	$ch2Co->[$ch2Cnt][1] =~ s/^\s+|\s+$//g;
	$ch2Co->[$ch2Cnt][2]=substr($line,46,8); 
	$ch2Co->[$ch2Cnt][2] =~ s/^\s+|\s+$//g;
	$ch2Res[$ch2Cnt]=$resName;
	$ch2Num[$ch2Cnt]=$resNum;
	$ch2Cnt++;
    }
}
# calculate interface residues based on atom-atom distances
$ch1IntCnt=0;
$ch2IntCnt=0;
for($i=0;$i<$ch1Cnt;++$i){
    for($j=0;$j<$ch2Cnt;++$j){
	    $dis=($ch1Co->[$i][0]-$ch2Co->[$j][0])*($ch1Co->[$i][0]-$ch2Co->[$j][0])
		+ ($ch1Co->[$i][1]-$ch2Co->[$j][1])*($ch1Co->[$i][1]-$ch2Co->[$j][1])
		+ ($ch1Co->[$i][2]-$ch2Co->[$j][2])*($ch1Co->[$i][2]-$ch2Co->[$j][2]);
	    if($dis<$CutOff){
		if($ch1IntCnt==0){
		    $ch1IF[$ch1IntCnt]=$ch1Num[$i];
		    $ch1IFRes[$ch1IntCnt]=$ch1Res[$i];
		    $ch1IntCnt++;
		}
		else {
		    if($ch1Num[$i]!=$ch1IF[$ch1IntCnt-1]){
			$ch1IF[$ch1IntCnt]=$ch1Num[$i];
			$ch1IFRes[$ch1IntCnt]=$ch1Res[$i];
			$ch1IntCnt++;			
		    }
		}
		if($ch2IntCnt==0){
		    $ch2IF[$ch2IntCnt]=$ch2Num[$j];
		    $ch2IFRes[$ch2IntCnt]=$ch2Res[$j];
		    $ch2IntCnt++;
		}
		else {
		    if($ch2Num[$j]!=$ch2IF[$ch2IntCnt-1]){
			$ch2IF[$ch2IntCnt]=$ch2Num[$j];
			$ch2IFRes[$ch2IntCnt]=$ch2Res[$j];
			$ch2IntCnt++;
		    }
		}
	    }    
	}
    }
# output results
    for($i=0;$i<$ch1IntCnt;++$i){
	printf OUT "$Chain1 $ch1IF[$i] $ch1IFRes[$i]\n";
    }
    for($i=0;$i<$ch2IntCnt;++$i){
	printf OUT "$Chain2 $ch2IF[$i] $ch2IFRes[$i]\n";
    }
