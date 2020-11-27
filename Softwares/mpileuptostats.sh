## FROM SAMTOOLS MPILEUP TO STATS ##
## Eva Maria Novoa, 2019 ##

if [[ $1 == "-h" ]]; then
echo "Usage: $0 <bam>"
echo "Expects mpileup in the same folder, in the form of <file.bam>.mpileup"
fi

# INPUT
i=$1

# SCRIPT

###############################
## PART 1: QUALITIES
###############################

echo "# Convert qualities to ordinal numbers"
python scripts/1.base_quality_chr2ord.keep_all.py $i.mpileup > $i.mpileup.digit

echo "# Get mean and median quality scores (produces output with header)"
python scripts/2.get_mean_median_qualities.py $i.mpileup.digit > $i.mean_median_quals

###############################
## PART 2: COVERAGE AND REF_NUC
###############################

echo "# Get general file info (coverage, ref_nuc)"
echo -ne "chr\tpos\tref_nuc\tcoverage\n" > $i.info
cut -f 1-4 $i.mpileup >> $i.info

###############################
## PART 3: RT-STOP
###############################

echo "# Get number of rt-stop, insertion and deletions (produces file with header)"
echo "0" > $i.mpileup.digit.rtstop.tmp
cat $i.mpileup.digit | awk {'print $5'} | awk -F "^" '{print NF-1}' >> $i.mpileup.digit.rtstop.tmp 
cat $i.mpileup.digit | awk {'print $5'} | awk -F "-" '{print NF-1}'  > $i.mpileup.digit.del
cat $i.mpileup.digit | awk {'print $5'} | awk -F "+" '{print NF-1}'  > $i.mpileup.digit.ins
sed '$ d' $i.mpileup.digit.rtstop.tmp > $i.mpileup.digit.rtstop # remove last line as we pushed everything 1
echo -ne "rtstop\tins\tdel\n" > $i.3moreinfo
paste $i.mpileup.digit.rtstop $i.mpileup.digit.del $i.mpileup.digit.ins >> $i.3moreinfo


###############################
## PART 6: CLEAN UP
###############################

# Merge all files
paste $i.info $i.mean_median_quals $i.3moreinfo > $i.STATS      

# Remove temporary files
rm $i.info $i.mean_median_quals $i.3moreinfo $i.mpileup.digit.rtstop $i.mpileup.digit.del $i.mpileup.digit.ins $i.mpileup.digit.rtstop.tmp $i.mpileup.digit sample1.txt

# Script completed!
echo "Script completed! Have a nice day! :-)"