#!/bin/bash
unset
##########################################################################################################################
EvoProD=`pwd`											
AADB=$EvoProD/data/CS-SCORE_REFERENCE.dat							
ASDB=$EvoProD/data/CSS-SCORE_REFERENCE.dat							
##########################################################################################################################
JOB=$1
################################# CHECK IF THE INPUT IS SEQUENCE ONLY OR SEQ+SS  #########################################
cd $EvoProD/$JOB
SQ=$(grep "^>$JOB\:SEQ" $JOB.fasta -A1 | grep -v "^>" | tr -d '[[:space:]]')
SS=$(grep "^>$JOB\:SEC" $JOB.fasta -A1 | grep -v "^>" | tr -d '[[:space:]]')
if [ ! -z $SQ ] && [ -z $SS ];
 then
   #echo "Secondary Structure in not available for $JOB Sequence, Please perform Secondary Structure Prediction..."
   #echo "If the SS is provided, please provide it in the format mentioned in README"
   exit
#else
   #echo "Sequence and Secondary Structure are available for $JOB.fasta"
fi
###########################################################################################################################
######################         CONVERTING THE INPUT PROTEIN SS from 8-CLASS TO 3-CLASS ####################################
SQ=$(grep "^>$JOB\:SEQ" $JOB.fasta -A1 | grep -v "^>" | tr -d '[[:space:]]' | sed 's/./& /g')
SS=$(grep "^>$JOB\:SEC" $JOB.fasta -A1 | grep -v "^>" | tr -d '[[:space:]]' | sed 's:[G,I]:H:g' | sed 's:[b,B]:E:g' | sed 's:[T,S]:C:g')
TMP=$(echo $SS | sed 's:CHC:CCC:g' | sed 's:CEC:CCC:g' | sed 's:CEH:CCH:g' | sed 's:CHE:CCE:g' | sed 's:EHC:ECC:g' | sed 's:EHE:ECE:g' | sed 's:HEC:HCC:g' | sed 's:HEH:HCH:g')
SS=$(echo $TMP | sed 's/./& /g')

declare -a SOQ=($SQ)
LEN1=$(echo $SQ | tr -d '[[:space:]]'| wc -c)
declare -a SOC=($SS)
LEN2=$(echo $SS | tr -d '[[:space:]]'| wc -c)
if [ "$LEN1" -ne "$LEN2" ]; 
 then
     echo "Sequence and Coresponding Secondary Structure are not matching in length"
     exit
fi
###########################################################################################################################
#################    RUNING COMPETENCIES SCORE COMPUTATION FROM PRE-COMPILED DATABASES    #################################
if [ -f $JOB.SQSS ]; then rm $JOB.SQSS; fi
MAX=$(expr $LEN1 - 2)
for (( X=1; X <= "$MAX"; X++ ))
    do
     N1=$(expr $X - 1); C1=$(expr $X + 1);
     NTR=$(echo ${SOQ[@]:$N1:1} | tr -d '[[:space:]]')					# N-terminal Residue
     NTS=$(echo ${SOC[@]:$N1:1} | tr -d '[[:space:]]')					# N-terminal Secondary Structure

     CTR=$(echo ${SOQ[@]:$C1:1} | tr -d '[[:space:]]')					# C-terminal Residue
     CTS=$(echo ${SOC[@]:$C1:1} | tr -d '[[:space:]]')					# C-terminal Secondary Structure

     MDR=$(echo ${SOQ[@]:$X:1} | tr -d '[[:space:]]')					# Middle Residue
     MDS=$(echo ${SOC[@]:$X:1} | tr -d '[[:space:]]')                           	# Middle Secondary Structure

     TRI=$(echo ${SOQ[@]:$N1:3} | tr -d '[[:space:]]')					# Tripeptide Sequence
     SSS=$(echo ${SOC[@]:$N1:3} | tr -d '[[:space:]]')					# Tripeptide SS

     CTRI=$(grep -m1 "^$TRI=" $AADB | awk -F'=' '{printf "%0.2f\n", $2}')   	        # Competency Score of tripeptide pre-computed 
     CSSS=$(grep -m1 "^$TRI\-$SSS\=" $ASDB | awk -F'=' '{printf "%0.2f\n", $2}')        # Competency Score of tripeptide-SecStr pre-computed
     echo -e "$NTR$MDR$CTR-$NTS$MDS$CTS\t\t$CTRI\t\t$CSSS" >> $JOB.SQSS
   done
   CSSQPR=$(awk '{sum += $2; sumsq += ($2)^2} END {printf "%0.3f %0.3f \n", sum/NR, sqrt((sumsq-sum^2/NR)/NR)}' $JOB.SQSS)
   CSSSPR=$(awk '{sum += $3; sumsq += ($3)^2} END {printf "%0.3f %0.3f \n", sum/NR, sqrt((sumsq-sum^2/NR)/NR)}' $JOB.SQSS)
   PRBCS=$(awk '{if ($2 < 32.15) print}' $JOB.SQSS | wc -l  | awk -v LN=$LEN1 '{printf "%0.2F\n", (100*$1/LN)}')
   PRBCSS=$(awk '{if ($3 < 15.50) print}' $JOB.SQSS | wc -l | awk -v LN=$LEN1 '{printf "%0.2F\n", (100*$1/LN)}')
   #echo -e "$JOB\t\t$CSSQPR\t$PRBCS\t\t$CSSSPR\t$PRBCSS"  > $JOB.DAT 
   echo -e "Name of Input Protein:\t$1.fasta" > $JOB.DAT
   echo -e "Average CS-Score for $1.fasta :\t $CSSQPR" >> $JOB.DAT
   echo -e "Percentage Residues Scoring Below CS-Score Threshold :\t$PRBCS" >> $JOB.DAT
   echo -e "Average CSS-Score for $1.fasta :\t $CSSSPR" >> $JOB.DAT
   echo -e "Percentage Residues Scoring Below CSS-Score Threshold :\t$PRBCSS" >> $JOB.DAT
