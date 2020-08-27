#!/usr/bin/perl
use strict; use warnings;
use POSIX;
#use List::Util 'shuffle';

my $datetime = localtime();  
print "Checkpoint 1 - Current date and time according to the system : $datetime\n";  

#Getting arguments from command line
my ($Human_fastq1, $Human_fastq2, $Mouse_fastq1, $Mouse_fastq2, $chimera_genome, $barcode_whitelist, $gtf_human, $gtf_mouse, $mouse_index, $human_index, $threshold, $root_dir) = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4], $ARGV[5], $ARGV[6], $ARGV[7], $ARGV[8], $ARGV[9], $ARGV[10], $ARGV[11]);
my $system_cmd;
$system_cmd = system("mkdir $root_dir");

#Paths to the genome index directories, whitelist barcodes and human-mouse annotation gtf files
my $genome_dir = $chimera_genome; #"/home/msnaveed/sra_local_repo/chimera_index/v3";
my $white_list = $barcode_whitelist; #"/home/msnaveed/sra_local_repo/10x_genomics/*.txt";
my $base_output_dir = "$root_dir/output_dir";
my $species_output_dir = "";

my $human_gtf = $gtf_human; #"/home/msnaveed/sra_local_repo/chimera_genome/human_genome/v3/*.gtf";
my $mouse_gtf = $gtf_mouse; #"/home/msnaveed/sra_local_repo/chimera_genome/altered_mouse_genome/v3/*.gtf";





$system_cmd = system("mkdir $base_output_dir");

# Run STAR on each human and mouse dataset
for my $specie (qw(human mouse)) {
	
	$species_output_dir = "$base_output_dir/$specie/";
	$system_cmd = system("mkdir $species_output_dir");
	if($specie eq "mouse"){
		$system_cmd = system("./runSTARAnalysis.sh $genome_dir $white_list $Mouse_fastq1 $Mouse_fastq2 $species_output_dir > $base_output_dir/Log_STARAnalysis_Mouse.txt");
		
	}
	elsif ($specie eq "human") {
		$system_cmd = system("./runSTARAnalysis.sh $genome_dir $white_list $Human_fastq1 $Human_fastq2 $species_output_dir > $base_output_dir/Log_STARAnalysis_Human.txt");
	}

}

# Classify the barcodes for each STAR result set
$system_cmd = system("perl classify_barcodes_pipeline.pl $mouse_gtf $human_gtf $base_output_dir nonchimeric $threshold");


# Putting the human and mouse barcodes extracted in the previous steps into individual BAM files and then doing BAM to fastq conversion
for my $specie (qw(human mouse)) {
	
	$species_output_dir = "$base_output_dir/$specie/analysis";
	my $currentLabelledFile = "$species_output_dir/$specie" . "_barcodes_classification.csv";

	open my $BARCODE, "<$currentLabelledFile" or die "can't open $currentLabelledFile\n";

	# Reading the classified barcodes in the labelled file into a dictionary
	my %barcodes;
	my $header = <$BARCODE>;
	my @header_ar = split(",", $header);
	chomp(@header_ar);

	while (<$BARCODE>) {
		chomp;
		my @line = split(",", $_);
		my $barcode = $line[0];
		my $species = $line[1];
		$barcodes{$species}{$barcode}+=1;

	}

	if(exists $barcodes{"mouse"}){


		my $currentSamFile = "$species_output_dir/Aligned.sortedByCoord.out.bam";

		my $prem_currentSamFile = "$species_output_dir/MouseFilteredAligned.out.sam";

		#Writing the header of the sam file to the prem_currentSamFile
		$system_cmd = system("samtools view -H -@ 32 $currentSamFile > $prem_currentSamFile");
		
		#Opening the prem_currentSamFile for writing
		open my $OUT, ">>$prem_currentSamFile" or die "can't open $prem_currentSamFile\n";

		#Opening the currentSamFile for reading
		open my $IN, "samtools view -@ 32 $currentSamFile |" or die "can't open $currentSamFile\n";

		while (<$IN>) {

			my ($barcode) = $_ =~ /CB:Z:([ACGT]*)/;
			
			if ($barcode && (exists $barcodes{"mouse"}{$barcode})) {
				print $OUT $_;
			}
		
		}
		close $IN;
		close $OUT;

		$system_cmd = system("samtools view -b -@ 32 $prem_currentSamFile > $prem_currentSamFile.bam");
		$system_cmd = system("rm $prem_currentSamFile");


		my $sys1 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Mouse_fastq1 $species_output_dir/mouse_R2.fastq";
		my $sys2 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Mouse_fastq2 $species_output_dir/mouse_R1.fastq";
		
		foreach my $i (1, 2) {
			my $pid = fork();
			if ($pid==0) { # child
				if($i==1){
					exec("$sys1");
				}
				if($i==2){
					exec("$sys2");	
				}
				die "Exec $i failed: $!\n";
			} elsif (!defined $pid) {
				warn "Fork $i failed: $!\n";
			}
		}

		1 while wait() >= 0;

		#system("$sys1");
		#system("$sys2");

		

	}

	if(exists $barcodes{"human"}){

		my $currentSamFile = "$species_output_dir/Aligned.sortedByCoord.out.bam";

		my $prem_currentSamFile = "$species_output_dir/HumanFilteredAligned.out.sam";

		#Writing the header of the sam file to the filteredSamFile
		$system_cmd = system("samtools view -H -@ 32 $currentSamFile > $prem_currentSamFile");
		
		#Opening the filteredSamFile for writing
		open my $OUT, ">>$prem_currentSamFile" or die "can't open $prem_currentSamFile\n";

		#Opening the currentSamFile for reading
		open my $IN, "samtools view -@ 32 $currentSamFile |" or die "can't open $currentSamFile\n";

		while (<$IN>) {

			my ($barcode) = $_ =~ /CB:Z:([ACGT]*)/;
			
			if ($barcode && (exists $barcodes{"human"}{$barcode})) {
				print $OUT $_;
			}
		
		}
		close $IN;
		close $OUT;

		$system_cmd = system("samtools view -b -@ 32 $prem_currentSamFile > $prem_currentSamFile.bam");
		$system_cmd = system("rm $prem_currentSamFile");

		my $sys1 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Human_fastq1 $species_output_dir/human_R2.fastq";
		my $sys2 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $Human_fastq2 $species_output_dir/human_R1.fastq";
		
		#Running the conversion of BAM to fastq in parallel because we need to individualy convert for R1 and R2. Converting in parallel would save us a lot of time
		#Logic: Creating two child processes and letting these complete
		foreach my $i (1, 2) {
			my $pid = fork();
			if ($pid==0) { # child
				if($i==1){
					exec("$sys1"); # Exec shall exit the child process
				}
				if($i==2){
					exec("$sys2"); # Exec shall exit the child process
				}
				die "Exec $i failed: $!\n";
			} elsif (!defined $pid) {
				warn "Fork $i failed: $!\n";
			}
		}

		1 while wait() >= 0;

		#system("$sys1");
		#system("$sys2");

	}
	
}

$datetime = localtime();  
print "Checkpoint 2 - Current date and time according to the system : $datetime\n";  

#Path to the diretory containing the final results
my $results_dir = "$root_dir/results";
my $exp_dir = "";
my $mouse_genome_dir = $mouse_index; #"/home/msnaveed/sra_local_repo/chimera_genome/mouse_index/v3";
my $human_genome_dir = $human_index; #"/home/msnaveed/sra_local_repo/chimera_genome/human_index/v3";

system("mkdir $results_dir");

# Running STAR on each human and mouse set produced from individual raw human and mouse dataset input i.e. 
# 3 - Mouse Fasqs input -> 3h(human) and 3m(mouse) dataset
# 4 - Human Fasqs input -> 4h(human) and 4m(mouse) dataset
for my $exp (qw(3h 3m 4h 4m)) {
	system("mkdir $results_dir/$exp");
	$exp_dir = "$results_dir/$exp/";
	if($exp eq "3h"){
		if (-e "$base_output_dir/mouse/analysis/human_R1.fastq.gz") { #Checking if there were any human reads founds inside the dataset 3
			
			$system_cmd = system("./runFinalSTAR.sh $human_genome_dir $white_list '$base_output_dir/mouse/analysis/human_R2.fastq.gz' '$base_output_dir/mouse/analysis/human_R1.fastq.gz' $exp_dir > $results_dir/3hLog_STARAnalysis.txt");
		}
		else{
			print ("3h fastq doesn't exist\n");
		}
	}
	if($exp eq "3m"){ 
		if (-e "$base_output_dir/mouse/analysis/mouse_R1.fastq.gz") { #Checking if there were any mouse reads founds inside the dataset 3
			$system_cmd = system("./runFinalSTAR.sh $mouse_genome_dir $white_list '$base_output_dir/mouse/analysis/mouse_R2.fastq.gz' '$base_output_dir/mouse/analysis/mouse_R1.fastq.gz' $exp_dir > $results_dir/3mLog_STARAnalysis.txt");
		}
		else{
			print ("3m fastq doesn't exist\n");
		}
	}
	if($exp eq "4m"){
		if (-e "$base_output_dir/human/analysis/mouse_R1.fastq.gz") { #Checking if there were any mouse reads founds inside the dataset 4
			$system_cmd = system("./runFinalSTAR.sh $mouse_genome_dir $white_list '$base_output_dir/human/analysis/mouse_R2.fastq.gz' '$base_output_dir/human/analysis/mouse_R1.fastq.gz' $exp_dir > $results_dir/4mLog_STARAnalysis.txt");
		}
		else{
			print ("4m fastq doesn't exist\n");
		}
	}
	if($exp eq "4h"){ 
		if (-e "$base_output_dir/human/analysis/human_R1.fastq.gz") { #Checking if there were any human reads founds inside the dataset 4
			$system_cmd = system("./runFinalSTAR.sh $human_genome_dir $white_list '$base_output_dir/human/analysis/human_R2.fastq.gz' '$base_output_dir/human/analysis/human_R1.fastq.gz' $exp_dir > $results_dir/4hLog_STARAnalysis.txt");
		}
		else{
			print ("4h fastq doesn't exist\n");
		}
	}
}


$datetime = localtime();  
print "Checkpoint 3 - Current date and time according to the system : $datetime\n";  

my %barcodes;

#Reading the labelled csv file to a dictionary %barcodes
for my $specie (qw(human mouse)) {
	
	$species_output_dir = "$base_output_dir/$specie/analysis";
	my $currentLabelledFile = "$species_output_dir/$specie" . "_barcodes_classification.csv";

	open my $BARCODE, "<$currentLabelledFile" or die "can't open $currentLabelledFile\n";

	
	my $header = <$BARCODE>;
	my @header_ar = split(",", $header);
	chomp(@header_ar);

	while (<$BARCODE>) {
		chomp;
		my @line = split(",", $_);
		my $barcode = $line[0];
		my $species = $line[1];
		$barcodes{$specie}{$species}{$barcode}+=1;
	}

	close $BARCODE;
}


$datetime = localtime();   
print "Checkpoint 4 - Current date and time according to the system : $datetime\n";  

my $exp5files = "$base_output_dir/exp5files";
system("mkdir $exp5files");


#Opening up individual human and mouse BAM files that we formed by running the STARSolo on the intiail raw input fastqs and then using the threshold to classify the barcodes belonging to either of the species
open my $MOUSE, "samtools view -@ 32 $base_output_dir/mouse/analysis/MouseFilteredAligned.out.sam.bam |" or die "can't open $base_output_dir/mouse/analysis/MouseFilteredAligned.out.sam.bam\n";
open my $HUMAN, "samtools view -@ 32 $base_output_dir/human/analysis/HumanFilteredAligned.out.sam.bam |" or die "can't open $base_output_dir/human/analysis/HumanFilteredAligned.out.sam.bam\n";


#Putting all the reads of mouse barcodes that are present in the csv files into a dictionary %mouse_barcodes
my %mouse_barcodes;
while (<$MOUSE>) {
	
	my ($barcode) = $_ =~ /CB:Z:([ACGT]*)/;
	
	if ($barcode && exists($barcodes{"mouse"}{"mouse"}{$barcode})) {
		$mouse_barcodes{$barcode}{$_}=1;
	}
    
}

close $MOUSE;



my %only_mb; # The barcode dictionary for mouse barcodes only

#Creating a directory to store each mouse barcode as a file and the reads associated with that barcode within that file
system("mkdir mb");
my $mb = "mb";
for my $barcode (keys %mouse_barcodes) {
	open my $OUT, ">>$mb/$barcode" or die "can't open $mb/$barcode\n";
	for my $line (keys %{$mouse_barcodes{$barcode}}){
		print $OUT $line;
	}
	$only_mb{$barcode}=1;
	close $OUT;
}

#No longer need %mouse_barcodes dict as we already have %only_mb
%mouse_barcodes = ();
undef %mouse_barcodes;


$datetime = localtime();   
print "Checkpoint 5 - Current date and time according to the system : $datetime\n";  


#Putting all the reads of human barcodes that are present in the csv files into a dictionary %human_barcodes
my %human_barcodes;
while (<$HUMAN>) {
	
	my ($barcode) = $_ =~ /CB:Z:([ACGT]*)/;
	
	if ($barcode && exists($barcodes{"human"}{"human"}{$barcode})) {
		$human_barcodes{$barcode}{$_}=1;
	}
}

close $HUMAN;


my %only_hb; # The barcode dictionary for human barcodes only

#Creating a directory to store each human barcode as a file and the reads associated with that barcode within that file
system("mkdir hb");
my $hb = "hb";

for my $barcode (keys %human_barcodes) {
	open my $OUT, ">>$hb/$barcode" or die "can't open $hb/$barcode\n";
	for my $line (keys %{$human_barcodes{$barcode}}){
		print $OUT $line;
	}
	$only_hb{$barcode}=1;
	close $OUT;
	close $OUT;
}

#No longer need %human_barcodes dict as we already have %only_hb
%human_barcodes = ();
undef %human_barcodes;

$datetime = localtime();  
print "Checkpoint 6 - Current date and time according to the system : $datetime\n";   



my $tot = 0;
my $human_reads = 0;
my $mouse_reads = 0;
my $count = 0;


#Producing chimeric-exp5 dataset by mixing proportions of human and mouse reads
# exp5.1: 50% human - 50% mouse
# exp5.2: 30% human - 70% mouse
# exp5.3: 10% human - 90% mouse
for my $exp (qw(5.1 5.2 5.3)) {
	if($exp eq "5.1"){
		$human_reads = keys %{$barcodes{"human"}{"human"}}; #Finding ho many human reads we have in total from the original raw input
		$human_reads = ceil($human_reads*0.5); #Rounding up to the nearest integer
		$mouse_reads = keys %{$barcodes{"mouse"}{"mouse"}}; #Finding ho many mouse reads we have in total from the original raw input
		$mouse_reads = ceil($mouse_reads*0.5); #Rounding up to the nearest integer
		print("Number of Human Barcodes in $exp: $human_reads\n");
		print("Number of Mouse Barcodes in $exp: $mouse_reads\n");
	}
	if($exp eq "5.2"){
		$human_reads = keys %{$barcodes{"human"}{"human"}};
		$human_reads = ceil($human_reads*0.3);
		$mouse_reads = keys %{$barcodes{"mouse"}{"mouse"}};
		$mouse_reads = ceil($mouse_reads*0.7);
		print("Number of Human Barcodes in $exp: $human_reads\n");
		print("Number of Mouse Barcodes in $exp: $mouse_reads\n");
	}
	if($exp eq "5.3"){
		$human_reads = keys %{$barcodes{"human"}{"human"}};
		$human_reads = ceil($human_reads*0.1);
		$mouse_reads = keys %{$barcodes{"mouse"}{"mouse"}};
		$mouse_reads = ceil($mouse_reads*0.9);
		print("Number of Human Barcodes in $exp: $human_reads\n");
		print("Number of Mouse Barcodes in $exp: $mouse_reads\n");
	}

	# Getting BAM and BAM-fastq data ready for both human and mouse 
	for my $specie (qw(mouse human)){
		$species_output_dir = "$base_output_dir/$specie/analysis";
		my $currentSamFile  = "";
		if($specie eq "mouse"){
			$currentSamFile = "$species_output_dir/MouseFilteredAligned.out.sam.bam"; #Name of the Mouse bam file to be produced
		}
		if($specie eq "human"){
			$currentSamFile = "$species_output_dir/HumanFilteredAligned.out.sam.bam"; #Name of the Human bam file to be produced
		}
		

		my $prem_currentSamFile = "$exp5files/$exp.$specie" . ".sam";
		#Writing the header of the sam file to the filteredSamFile
		#system("samtools view -H -@ 32 $currentSamFile > $prem_currentSamFile");
		my $output = qx/samtools view -H -@ 32 $currentSamFile > $prem_currentSamFile/;

		#Opening the filteredSamFile for writing
		#open my $OUT, ">>$prem_currentSamFile" or die "can't open $prem_currentSamFile\n";

		#system("samtools view -H -@ 32 $currentSamFile > $prem_currentSamFile");

		
		if($specie eq "mouse"){
			$tot = $mouse_reads;
			$count = 0;
			for my $barcode (keys %only_mb) {
				system("cat $mb/$barcode >> $prem_currentSamFile"); #Copy the all the barcode reads present in the file to the sam file
				#print $OUT $line;
				$count+=1;
				if($count==$tot){ # Do it until the proportion of reads is reached
					last;
				}
			}
		}

		print("Passed Mouse\n");

		if($specie eq "human"){
			$tot = $human_reads;
			$count = 0;
			for my $barcode (keys %only_hb) {
				system("cat $hb/$barcode >> $prem_currentSamFile");
				#print $OUT $line;
				$count+=1;
				if($count==$tot){
					last;
				}
			}
		}
		#close $OUT;
		
		print("Passed Human\n");
		

		system("samtools view -b -@ 32 $prem_currentSamFile > $prem_currentSamFile.bam"); #Convert SAM to BAM
		system("rm $prem_currentSamFile"); # Remove the SAM file as we already have the BAM version

		my $curr_fastq1 = "";
		my $curr_fastq2 = "";
		if($specie eq "mouse"){
			$curr_fastq1 = $Mouse_fastq1;
			$curr_fastq2 = $Mouse_fastq2;
		}
		if($specie eq "human"){
			$curr_fastq1 = $Human_fastq1;
			$curr_fastq2 = $Human_fastq2;
		}

		# Bam to fastq conversion commands for R1 and R2 separately
		my $sys1 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $curr_fastq1 $exp5files/$exp.$specie._R2.fastq";
		my $sys2 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $curr_fastq2 $exp5files/$exp.$specie._R1.fastq";
		
		foreach my $i (1, 2) {
			my $pid = fork();
			if ($pid==0) { # child
				if($i==1){
					exec("$sys1");
				}
				if($i==2){
					exec("$sys2");	
				}
				die "Exec $i failed: $!\n";
			} elsif (!defined $pid) {
				warn "Fork $i failed: $!\n";
			}
		}

		1 while wait() >= 0;

		#system("$sys1");
		#system("$sys2");
	}

	print("Got out of the loop\n");
	system("cat $exp5files/$exp.mouse._R1.fastq.gz $exp5files/$exp.human._R1.fastq.gz >> $exp5files/$exp._R1.fastq.gz"); #Merging the individual human and mouse R1 dataset produced to create exp5 data
	system("cat $exp5files/$exp.mouse._R2.fastq.gz $exp5files/$exp.human._R2.fastq.gz >> $exp5files/$exp._R2.fastq.gz"); #Merging the individual human and mouse R2 dataset produced to create exp5 data
	

}

#Removing the directories containing the each barcode as a file
system("rm -rf mb hb");


$datetime = localtime();   
print "Checkpoint 7 - Current date and time according to the system : $datetime\n";  
print("Done creating first round of exp5 fastq files-chimeric\n");

#Running initial STAR on chimeric - exp 5s to filter out human and mouse reads from each experiment set
for my $specie (qw(5.1 5.2 5.3)) {
	
	$species_output_dir = "$base_output_dir/$specie/";
	$system_cmd = system("mkdir $species_output_dir");

	$system_cmd = system("./runSTARAnalysis.sh $genome_dir $white_list $exp5files/$specie._R2.fastq.gz $exp5files/$specie._R1.fastq.gz $species_output_dir > $base_output_dir/Log_STARAnalysis_$specie.txt");
		
}

#Producing a spare matrix of barcode-gene and a labelled file for each STAR Experiment where each barcode is classified on the threshold level
$system_cmd = system("perl classify_barcodes_pipeline.pl $mouse_gtf $human_gtf $base_output_dir chimera $threshold");


#Creating individual human and mmouse sets out of the initial STARsolo results on each exp5 dataset
for my $specie (qw(5.1 5.2 5.3)) {
	
	$species_output_dir = "$base_output_dir/$specie/analysis";
	my $currentLabelledFile = "$species_output_dir/$specie" . "_barcodes_classification.csv";

	open my $BARCODE, "<$currentLabelledFile" or die "can't open $currentLabelledFile\n";

	#Reading barcode labelled csv file
	my %barcodes;
	my $header = <$BARCODE>;
	my @header_ar = split(",", $header);
	chomp(@header_ar);

	while (<$BARCODE>) {
		chomp;
		my @line = split(",", $_);
		my $barcode = $line[0];
		my $species = $line[1];
		$barcodes{$species}{$barcode}+=1;

	}

	#If there are mouse barcodes present in the exp5 dataset, then produce a bam file and then do BAM-to-fastq conversion
	if(exists $barcodes{"mouse"}){


		my $currentSamFile = "$species_output_dir/Aligned.sortedByCoord.out.bam";

		my $prem_currentSamFile = "$species_output_dir/MouseFilteredAligned.out.sam";

		#Writing the header of the sam file to the filteredSamFile
		$system_cmd = system("samtools view -H -@ 32 $currentSamFile > $prem_currentSamFile");
		
		#Opening the filteredSamFile for writing
		open my $OUT, ">>$prem_currentSamFile" or die "can't open $prem_currentSamFile\n";

		#Opening the currentSamFile for reading
		open my $IN, "samtools view -@ 32 $currentSamFile |" or die "can't open $currentSamFile\n";

		while (<$IN>) {

			my ($barcode) = $_ =~ /CB:Z:([ACGT]*)/;
			
			if ($barcode && (exists $barcodes{"mouse"}{$barcode})) {
				print $OUT $_;
			}
		
		}
		close $IN;
		close $OUT;

		$system_cmd = system("samtools view -b -@ 32 $prem_currentSamFile > $prem_currentSamFile.bam");
		$system_cmd = system("rm $prem_currentSamFile");


		
		my $sys1 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $exp5files/$specie._R2.fastq.gz $species_output_dir/$specie.mouse_R2.fastq";
		my $sys2 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $exp5files/$specie._R1.fastq.gz $species_output_dir/$specie.mouse_R1.fastq";
		
		#Child processes to achieve parallelism
		foreach my $i (1, 2) {
			my $pid = fork();
			if ($pid==0) { # child
				if($i==1){
					exec("$sys1");
				}
				if($i==2){
					exec("$sys2");	
				}
				die "Exec $i failed: $!\n";
			} elsif (!defined $pid) {
				warn "Fork $i failed: $!\n";
			}
		}

		1 while wait() >= 0;

		#system("$sys1");
		#system("$sys2");

	}

	#If there are human barcodes present in the exp5 dataset, then produce a bam file and then do BAM-to-fastq conversion
	if(exists $barcodes{"human"}){

		my $currentSamFile = "$species_output_dir/Aligned.sortedByCoord.out.bam";

		my $prem_currentSamFile = "$species_output_dir/HumanFilteredAligned.out.sam";

		#Writing the header of the sam file to the prem_currentSamFile
		$system_cmd = system("samtools view -H -@ 32 $currentSamFile > $prem_currentSamFile");
		
		#Opening the prem_currentSamFile for writing
		open my $OUT, ">>$prem_currentSamFile" or die "can't open $prem_currentSamFile\n";

		#Opening the currentSamFile for reading
		open my $IN, "samtools view -@ 32 $currentSamFile |" or die "can't open $currentSamFile\n";

		while (<$IN>) {

			my ($barcode) = $_ =~ /CB:Z:([ACGT]*)/;
			
			if ($barcode && (exists $barcodes{"human"}{$barcode})) {
				print $OUT $_;
			}
		
		}
		close $IN;
		close $OUT;

		$system_cmd = system("samtools view -b -@ 32 $prem_currentSamFile > $prem_currentSamFile.bam");
		$system_cmd = system("rm $prem_currentSamFile");

		
		my $sys1 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $exp5files/$specie._R2.fastq.gz $species_output_dir/$specie.human_R2.fastq";
		my $sys2 = "perl convertBamtToFastq.pl $prem_currentSamFile.bam $exp5files/$specie._R1.fastq.gz $species_output_dir/$specie.human_R1.fastq";
		
		
		foreach my $i (1, 2) {
			my $pid = fork();
			if ($pid==0) { # child
				if($i==1){
					exec("$sys1");
				}
				if($i==2){
					exec("$sys2");	
				}
				die "Exec $i failed: $!\n";
			} elsif (!defined $pid) {
				warn "Fork $i failed: $!\n";
			}
		}

		1 while wait() >= 0;

		#system("$sys1");
		#system("$sys2");
		

	}
	
	
}


$datetime = localtime();  
print "Checkpoint 8 - Current date and time according to the system : $datetime\n";  
my $exp = "";

#Final STARsolo run on individual human and mouse set produced from each exp5 set
for my $set (qw(5.1h 5.1m 5.2h 5.2m 5.3h 5.3m)) {
	system("mkdir $results_dir/$set");
	$exp_dir = "$results_dir/$set/";


	if(($set eq "5.1h") || ($set eq "5.2h") || ($set eq "5.3h")){ #If it is a human set, neeed to align against human only genome
		$exp = $set; #Can't change the readable-only variable of forloop so equating to a changeable dynamic variable $exp
		chop($exp); #Want to scrape the last letter out
		if (-e "$base_output_dir/$exp/analysis/$exp.human_R1.fastq.gz") {
			
			$system_cmd = system("./runFinalSTAR.sh $human_genome_dir $white_list '$base_output_dir/$exp/analysis/$exp.human_R2.fastq.gz' '$base_output_dir/$exp/analysis/$exp.human_R1.fastq.gz' $exp_dir > $results_dir/$exp.Log_STARAnalysis.txt");
		}
		else{ #Just in case human reads were not present in a particular exp5 set
			print ("$set fastq doesn't exist");
		}
	}
	
	if(($set eq "5.1m") || ($set eq "5.2m") || ($set eq "5.3m")){ #If it is a mouse set, neeed to align against mouse only genome
		$exp = $set; #Can't change the readable-only variable of forloop so equating to a changeable dynamic variable $exp
		chop($exp); #Want to scrape the last letter out
		if (-e "$base_output_dir/$exp/analysis/$exp.mouse_R1.fastq.gz") {
			
			$system_cmd = system("./runFinalSTAR.sh $mouse_genome_dir $white_list '$base_output_dir/$exp/analysis/$exp.mouse_R2.fastq.gz' '$base_output_dir/$exp/analysis/$exp.mouse_R1.fastq.gz' $exp_dir > $results_dir/$exp.Log_STARAnalysis.txt");
		}
		else{ #Just in case mousereads were not present in a particular exp5 set
			print ("$set fastq doesn't exist");
		}
	}
}

$datetime = localtime();    
print "Checkpoint 9(final) - Current date and time according to the system : $datetime\n";  