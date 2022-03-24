#!/usr/bin/perl
##############################################################################
#                           SCASA WRAP
#             Single Cell Alternative Splicing Analyser
#             Version V1.0.1
#             Author: Lu Pan, Trung Nghia Vu, Yudi Pawitan
#             Last Update: 2022-03-24
##############################################################################
use strict;
use warnings;
use File::Basename;
use FindBin;
#use File::Find::Rule;
use lib "$FindBin::Bin";
use Cwd 'abs_path';
use Getopt::Std;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use Term::ANSIColor qw(:constants);

print BRIGHT_YELLOW, "\n##############################################################\n";
print BRIGHT_YELLOW, "#\tSCASA V1.0.1\n";
print BRIGHT_YELLOW, "#\tSINGLE CELL TRANSCRIPT QUANTIFICATION TOOL\n";
print BRIGHT_YELLOW, "#\tVersion Date: 2022-03-24\n";
print BRIGHT_YELLOW, "#\tFOR ANY ISSUES, CONTACT: LU.PAN\@KI.SE\n";
print BRIGHT_YELLOW, "#\thttps://github.com/eudoraleer/scasa/\n";
print BRIGHT_YELLOW, "##############################################################\n";

# conda install -c bioconda salmon
# required: coreutils, salmon, kallisto, bustools

print GREEN, "\n", RESET;
# print GREEN, ON_CYAN, BLINK, "garish!\n";

my $scasa_dir = `which scasa`;
# my $scasa_dir = "/Users/luadmpan/Dropbox/KI/Studies/Annual_Follow_Up/Half_Time/Manuscript/SCASA/scasa";
chomp $scasa_dir;
$scasa_dir=~/(.*)\/scasa$/i;
$scasa_dir = $1;
my $scasa_xmatrix_script1 = "$scasa_dir/SCRIPTS/Xmatrix_Gen/scasa_simulate_xmatrix_step1_v1.0.0.R";
my $scasa_xmatrix_script2 = "$scasa_dir/SCRIPTS/Xmatrix_Gen/scasa_simulate_xmatrix_step2_v1.0.0.R";
my $scasa_xmatrix_script31 = "$scasa_dir/SCRIPTS/Xmatrix_Gen/scasa_simulate_xmatrix_step3_1_v1.0.0.R";
my $scasa_xmatrix_script32 = "$scasa_dir/SCRIPTS/Xmatrix_Gen/scasa_simulate_xmatrix_step3_2_v1.0.0.R";
my $scasa_xmatrix_script4 = "$scasa_dir/SCRIPTS/Xmatrix_Gen/scasa_simulate_xmatrix_step4_v1.0.0.R";

my $scasa_quant_script11 = "$scasa_dir/SCRIPTS/scasa_quant_step1_1_v1.0.0.R";
my $scasa_quant_script12 = "$scasa_dir/SCRIPTS/scasa_quant_step1_2_v1.0.0.R";
my $scasa_quant_script2 = "$scasa_dir/SCRIPTS/scasa_quant_step2_v1.0.0.R";
my $scasa_quant_script3 = "$scasa_dir/SCRIPTS/scasa_quant_step3_v1.0.0.R";

my $scasa_Rsource = "$scasa_dir/SCRIPTS/Scasa_Rsource.R";
my $reformat = "reformat.sh";

my $scasa_ref = "$scasa_dir/REFERENCE/";
my $trans2gen = "$scasa_dir/REFERENCE/txp2gene_hg38.tsv";
my $xmatrix_alevin = "$scasa_dir/REFERENCE/Xmatrix_hg38_alevin.RData";
my $xmatrix_bustools = "$scasa_dir/REFERENCE/Xmatrix_hg38_bustools.RData";
my $isoform_info_alevin = "$scasa_dir/REFERENCE/isoform_groups_UCSC_hg38_alevin.RData";
my $isoform_info_bustools = "$scasa_dir/REFERENCE/isoform_groups_UCSC_hg38_bustools.RData";

my $ref = "";
my $index_path = "";

# Input parameters
my $help = 0;
my $mapper = "salmon_alevin";
my $project_name = "My_Project";
my $input_dir = "./";
my $root_out = "./";
my $to_index = "YES";
my $to_align = "YES";
my $to_quant = "YES";
my $createxmatrix = "NO";
my $xmatrix = "alevin";
my $samplesheet = "NULL";
my @fastq_files = ();
my $tech = "10xv3";
my $white_list = "";
my $expectCells = "none";
my $num_threads = 4;
my $postalign_dir = "";
my $extra_para = "";

GetOptions('help|h' => \$help,
           'mapper|m=s' => \$mapper,
           'project|p=s' => \$project_name,
           'in|i=s' => \$input_dir,
           'out|o=s'=> \$root_out,
           'createxmatrix|b=s'=> \$createxmatrix,
           'xmatrix|x=s' => \$xmatrix,
           'index|e=s' => \$to_index,
           'index_dir|d=s' => \$index_path,
           'align|a=s'=> \$to_align,
           'quant|q=s'=> \$to_quant,
           'postalign_dir|g=s' => \$postalign_dir,
           'samplesheet|s=s' => \$samplesheet,
           'fastq|f=s' => \@fastq_files,
           'ref|r=s' => \$ref,
           'tech|t=s' => \$tech,
           'whitelist|w=s' => \$white_list,
           'cellthreshold|c=s' => \$expectCells,
           'nthreads|n=i' => \$num_threads) or pod2usage(1);

pod2usage(1) if $help;

if(($ref eq "") | ($ref eq "NULL")){
    die RED, "ERROR: Please provide a reference.\n";
}
my $output_dir = $root_out;
my $time = strftime("%Y%m%d%H%M%S", localtime(time));
if($project_name !~ /.*\_202[0-9]+$/){
    $project_name = "$project_name\_$time";
    $output_dir = $root_out."/SCASA_".$project_name."/";
}elsif($project_name =~ /^SCASA_.*\_202[0-9]+$/i){
    $output_dir = $root_out."/$project_name/";
}
$input_dir =~s/\/\//\//g;
$to_index = uc($to_index);
$createxmatrix = uc($createxmatrix);
$to_align = uc($to_align);
$to_quant = uc($to_quant);
$tech = lc($tech);
my $presets_out = $output_dir."0PRESETS/";
my $align_out = $output_dir."1ALIGN/";

if($createxmatrix eq "NO"){
    if(lc($xmatrix) eq "alevin"){
        $xmatrix = $xmatrix_alevin;
    }elsif (lc($xmatrix) eq "bustools"){
        $xmatrix = $xmatrix_bustools;
    }
}elsif(($createxmatrix eq "YES") & ((lc($xmatrix) eq "alevin")|(lc($xmatrix) eq "bustools"))){
    die RED, "ERROR: Please make a decision if you want to create Xmatrix or use prebuilt Xmatrix.\n";
}

@fastq_files = split(/,/,join(",",@fastq_files));
if(($to_align eq "NO") & ($postalign_dir ne "")){
  $align_out = $postalign_dir;
}
my $quant_out = $output_dir."2QUANT/";
my $log_dir = $output_dir."LOG/";

if((lc($mapper) eq "salmon_alevin") | lc($mapper) eq ""){
    $mapper = "salmon alevin";
}elsif (lc($mapper) eq "kallisto_bus"){
    $mapper = "kallisto bus";
}else{
    pod2usage(1);
    die RED, "\nERROR: Check your mapper argument before proceeding the pipeline, else go to Scasa on Github (https://github.com/eudoraleer/scasa/wiki) for more information.\n";
}

my @alignout_names = ();
my $align_script = "";
my $current_name = "";
my $pid_path = "";
my $current_out = "";
my $start_run = "";
my $end_run = "";
my $run_time = "";

chdir $input_dir;

# my $pipe_start = time();

if (-d $root_out){
  print RED, "Directory $root_out already exists. Writing into existing directory..\n";
}else{
  system("mkdir $root_out");
}

if (-d $output_dir){
  print RED, "Directory $output_dir already exists. Writing into existing directory..\n";
}else{
  system("mkdir $output_dir");
}

if (-d $project_name){
  print RED, "Directory $output_dir already exists. Writing into existing directory..\n";
}else{
  system("mkdir $output_dir");
}

if($to_align eq "YES"){
if (-d $align_out){
  print RED, "Directory $align_out already exists. Writing into existing directory..\n";
}else{
  system("mkdir $align_out");
}
}

if (-d $quant_out){
  print RED, "Directory $quant_out already exists. Writing into existing directory..\n";
}else{
  system("mkdir $quant_out");
}

if (-d $log_dir){
	print RED, "Directory $log_dir already exists. Writing into existing directory..\n";
}else{
	system("mkdir $log_dir");
}

if (-d $presets_out){
    print RED, "Directory $presets_out already exists. Writing into existing directory..\n";
}else{
    system("mkdir $presets_out");
}

if(($index_path eq "") | ($index_path eq "NULL")){
    $index_path = "$presets_out/REF_INDEX/";
    if (-d $index_path){
      print RED, "Directory $index_path already exists. Writing into existing directory..\n";
    }else{
      system("mkdir $index_path");
    }
}

if($createxmatrix eq "YES"){
    my $xmatrix_fasta = "$presets_out/XMATRIX_SIM_CDNA.fa";
    my $xmatrix_fasta_index = "$presets_out/XMATRIX_SIM_CDNA.fa.idx";
    system("Rscript $scasa_xmatrix_script1 $presets_out $ref $xmatrix_fasta");
    
    my $xmatrix_dir = "$presets_out/SIMULATE_XMATRIX/";
    my $polyester_out =`ls $xmatrix_dir/*_1.fasta`;
    chomp $polyester_out;
    $polyester_out =~/(.*)_1.fasta/;
    my $name = $1;
        system("awk '\{print \> \(NR \% 2 \? \"$name.odd.txt\" \: \"$name.even.txt\"\) \}' $polyester_out");
        system("awk -F ' ' '{print \$1}' $name.odd.txt | awk -F '/' '{print \$2}' > $name.odd.clean.txt");
    if(($white_list eq "NULL") | ($white_list eq "")){
        die RED, "ERROR: Please provide a whitelist in order to proceed with Xmatrix generation.\n";
    }
        system("Rscript $scasa_xmatrix_script2 $white_list $xmatrix_dir $name.odd.clean.txt");
        system("paste -d '' $xmatrix_dir/sample_01.white.list.txt $xmatrix_dir/sample_01.umi.txt > $xmatrix_dir/sample_01.barcodes.txt");
        system("paste -d \\\\n $name.odd.txt $xmatrix_dir/sample_01.barcodes.txt > $name.temp.fasta");
        system("perl -i -pe 's/mate1:(\\d+)-(\\d+)/\"mate1:\$1-\".(\$1 + 28 -1)/ge' $name.temp.fasta");
    system("perl -i -pe 's/mate1:(\\d+)-(\\d+)/\"mate1:\$1-\".(\$1 + 28 -1)/ge' $name\_2.fasta");
    system("mv $name.temp.fasta $polyester_out");
    system("mv $name\_1.fasta $xmatrix_dir/simulation_sampling_01.fasta");
    system("mv $name\_2.fasta $xmatrix_dir/simulation_sampling_02.fasta");
    system("$reformat in1=$xmatrix_dir/simulation_sampling_01.fasta in2=$xmatrix_dir/simulation_sampling_02.fasta out1=$xmatrix_dir/simulation_sampling_01.fastq out2=$xmatrix_dir/simulation_sampling_02.fastq overwrite=true");
    system("rm $xmatrix_dir/simulation_sampling_01.fasta $xmatrix_dir/simulation_sampling_02.fasta $name.odd.clean.txt");
    system("rm $name.even.txt $name.odd.txt");
    
    if($mapper eq "salmon alevin"){
        system("salmon index -t $xmatrix_fasta -i $xmatrix_fasta_index");
        if($tech eq "10xv3"){
            system("$mapper -l ISR -1 $xmatrix_dir/simulation_sampling_01.fastq -2 $xmatrix_dir/simulation_sampling_02.fastq --chromiumV3 -i $xmatrix_fasta_index -p $num_threads -o $xmatrix_dir/simulation_sampling_alevin_output --tgMap $trans2gen --whitelist $white_list --dumpBfh --dumpFeatures --dumpMtx");
        }elsif(($tech eq "10xv1") | ($tech eq "10xv2")){
            system("$mapper -l ISR -1 $xmatrix_dir/simulation_sampling_01.fastq -2 $xmatrix_dir/simulation_sampling_02.fastq --chromium -i $xmatrix_fasta_index -p $num_threads -o $xmatrix_dir/simulation_sampling_alevin_output --tgMap $trans2gen --whitelist $white_list --dumpBfh --dumpFeatures --dumpMtx");
        }
        
        system("Rscript $scasa_xmatrix_script31 $xmatrix_dir/simulation_sampling_alevin_output/alevin/bfh.txt $xmatrix_dir/TxCellID_mapping.RData $num_threads $xmatrix_dir");
        
    }elsif($mapper eq "kallisto bus"){
        system("kallisto index -i $xmatrix_fasta_index $xmatrix_fasta");
        system("$mapper -i $xmatrix_fasta_index -o $xmatrix_dir/simulation_sampling_bus_output/ -x $tech -t $num_threads $xmatrix_dir/simulation_sampling_01.fastq $xmatrix_dir/simulation_sampling_02.fastq");
        system("bustools text $xmatrix_dir/simulation_sampling_bus_output/output.bus -o $xmatrix_dir/simulation_sampling_bus_output/output.bus.txt");
        system("awk '{print \$1\$3 }' $xmatrix_dir/simulation_sampling_bus_output/output.bus.txt > $xmatrix_dir/simulation_sampling_bus_output/output.barcode_eq.bus.txt");
        system("Rscript $scasa_xmatrix_script32 $xmatrix_dir/simulation_sampling_bus_output/ $xmatrix_dir/TxCellID_mapping.RData $xmatrix_dir");
    }
    
    system("Rscript $scasa_xmatrix_script4 in=$xmatrix_dir/Xmatrix_eqClass.txt H=0.01 out=$xmatrix_dir/Xmatrix.RData Rsource=$scasa_Rsource ncore=$num_threads");
    $xmatrix = "$xmatrix_dir/Xmatrix.RData";
}

if($to_align eq "YES"){
    print GREEN, "\nPreparing for alignment..\n";
    
    if($mapper eq "salmon alevin"){
    if(($white_list ne "") & ($white_list ne "NULL")){
        $extra_para = "--whitelist $white_list $extra_para";
    }
    
    if($tech eq "10xv3"){
        $extra_para = "--chromiumV3 $extra_para";
    }elsif(($tech eq "10xv1") | ($tech eq "10xv2")){
        $extra_para = "--chromium $extra_para";
    }
    
    if($expectCells ne "none"){
        $extra_para = "--expectCells $expectCells $extra_para";
    }
    }
    
    if($to_index eq "YES"){
        print "Indexing reference..\n";
        if (-d $index_path){
          print RED, "Directory $index_path already exists. Writing into existing directory..\n";
        }else{
          system("mkdir $index_path");
        }
        if($mapper eq "kallisto bus"){
            system("kallisto index -i $index_path -k 31 $ref");
        }elsif($mapper eq "salmon alevin"){
            system("salmon index -t $ref -i $index_path");
        }else{
            die RED, "ERROR: Other alignment tools are not supported.\n";
            pod2usage(1);
        }
        print "Finnished indexing reference..\n";
    }
    
    print "Begins pseudo-alignment..\n";
    my $fastq_len = scalar @fastq_files;
    my @input_fastqs = ();
    if($samplesheet ne "NULL"){
    open SAMPLESHEET, "<$samplesheet" or die "Unable to open $samplesheet\n";
    while(my $a = <SAMPLESHEET>){
    chomp $a;
    my $current_fastqs =~s/,|\t|\s+/ /g;
    push(@input_fastqs, $current_fastqs);
    }
    close SAMPLESHEET;
    }elsif ($samplesheet eq "NULL"){
        if($fastq_len == 0){
        @fastq_files = `find $input_dir \\\( -name "*.fastq" -or -name "*.fq" -or -name "*.fq.gz" -or -name "*.fastq.gz" \\\)`;
        }
        chomp @fastq_files;
        my @r1_fastqs =  grep(/.*S*_L*_R1_0.*|.*_R1\..*|.*\.R1\..*|.*_R1_.*/, @fastq_files);
        @r1_fastqs =  sort @r1_fastqs;
        my @r2_fastqs =  grep(/.*S*_L*_R2_0.*|.*_R2\..*|.*\.R2\..*|.*_R2_.*/, @fastq_files);
        @r2_fastqs =  sort @r2_fastqs;
        
        if(scalar @r1_fastqs != scalar @r2_fastqs){
            print "FASTQ files found in your input directory did not come in pairs! Please check and resubmit.\n";
            exit;
        }else{
        for(my $i = 0; $i <= $#r1_fastqs; $i++){
            push(@input_fastqs, "$r1_fastqs[$i] $r2_fastqs[$i]");
        }
        }
    }
    
    $pid_path = $log_dir."scasa.align.$time.pid";
    my $align_all = $log_dir."scasa.$project_name.alignment.all.sh";
    open RUNALL, ">$align_all" or die "Unable to open $align_all\n";
    print RUNALL "#!/usr/bin/env bash\n";

    foreach my $i(@input_fastqs){

    $i =~/(.*)R1.*\.f.*q.*\s+.*/i;
    $current_name = $1;
    $current_name =~s/.*\/(.*)/$1/g;
    $current_name =~s/(.*)_$/$1/g;
    $current_name =~s/(.*)-$/$1/g;
    $current_name =~s/(.*)\.$/$1/g;
    $current_name = $1;
    my $temp_name = $current_name;
    $align_script = $log_dir."$project_name.$current_name.alignment.sh";

    push(@alignout_names, $current_name);
    $current_name = "$align_out/$current_name\_alignout/";
    
    if (-d $current_name){
          print RED, "Directory $current_name already exists. Writing into existing directory..\n";
        }else{
          system("mkdir $current_name");
    }
        
    print RUNALL "sh $align_script & echo \$\! >> $pid_path\n";
    open ALIGN, ">$align_script" or die "Unable to open $align_all\n";
    print ALIGN "#!/usr/bin/env bash\n";
    print ALIGN "cd $input_dir\n";
    print ALIGN "echo \"Begin alignment for sample $current_name..\"\n";
        if($mapper eq "kallisto bus"){
    print ALIGN "$mapper -i $index_path -o $current_name -x $tech -t $num_threads $i 2> $log_dir/$project_name.align.$temp_name.$time.o 1> $log_dir/$project_name.align.$temp_name.$time.e\n";
            print ALIGN "echo \"Finnish alignment for sample $temp_name!\"\n";
            
            print ALIGN "echo \"Preparing for correction step for sample $temp_name..\"\n";

            if($white_list eq ""){
            print ALIGN "bustools whitelist -o $log_dir/$temp_name\_whitelist.txt $current_name/output.sort.bus 2> $log_dir/$project_name.wl.$temp_name.$time.o 1> $log_dir/$project_name.wl.$temp_name.$time.e\n";
                print ALIGN "echo \"Finnish building whitelist for sample $temp_name!\"\n";
                $white_list = "$log_dir/$temp_name\_whitelist.txt";
            }

            print ALIGN "bustools correct -o $current_name/output.correct.bus -w $white_list $current_name/output.bus 2>$log_dir/$project_name.correct.$temp_name.$time.e 1> $log_dir/$project_name.correct.$temp_name.$time.o\n";
            print ALIGN "echo \"Finnish barcode correction for sample $temp_name!\"\n";
   print ALIGN "echo \"Begin bus sorting for sample $temp_name..\"\n";
            print ALIGN "bustools sort -t $num_threads -o $current_name/output.correct.sort.bus $current_name/output.correct.bus 2> $log_dir/$project_name.sort.$temp_name.$time.o 1> $log_dir/$project_name.sort.$temp_name.$time.e\n";
            print ALIGN "echo \"Finnish bus output sorting for sample $temp_name!\"\n";
            print ALIGN "bustools text $current_name/output.correct.sort.bus -o $current_name/output.correct.sort.bus.txt 2>$log_dir/$project_name.text.$temp_name.$time.o 1> $log_dir/$project_name.text.$temp_name.$time.e\n";
            print ALIGN "echo \"Finnish texting for sample $temp_name!\"\n";
            print ALIGN "awk -F '\t' '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$1\$2}' $current_name/output.correct.sort.bus.txt > $current_name/output.correct.sort.bus.final.txt 2>$log_dir/$project_name.text.$temp_name.$time.o 1> $log_dir/$project_name.text.$temp_name.$time.e\n";
            close ALIGN;
        }elsif($mapper eq "salmon alevin"){
            my @current_fastqs = split(/ /,$i);
            print ALIGN "$mapper -l ISR -1 $current_fastqs[0] -2 $current_fastqs[1] -i $index_path -p $num_threads -o $current_name --tgMap $trans2gen --dumpBfh --dumpFeatures --dumpMtx $extra_para 2> $log_dir/$project_name.align.$temp_name.$time.o 1> $log_dir/$project_name.align.$temp_name.$time.e\n";
            print ALIGN "echo \"Finnish alignment for sample $temp_name!\"\n";
    }

        close ALIGN;
}
    
close RUNALL;
$start_run = time();
system("nohup sh $align_all > $log_dir/$project_name.batch_align.o &");
sleep 30;
&checkPs("$pid_path");
$end_run = time();
$run_time = $end_run - $start_run;
print "Congratulations! Pseudo-alignment has completed in $run_time seconds!\n";

}

################### Quantification ########################
my $searchname = "";
if($to_quant eq "YES"){
	print "Scasa quantification has started..\n";
	$pid_path = $log_dir."scasa.quant.$time.pid";
	my $quant_all = $log_dir."scasa.$project_name.quant.all.sh";
	open RUNALL, ">$quant_all" or die "Unable to open $quant_all\n";
	print RUNALL "#!/usr/bin/env bash\n";
	my $quant_script = "";

    if($to_align eq "NO"){
        if($mapper eq "salmon alevin"){
            $searchname = "bfh.txt";
        }elsif ($mapper eq "kallisto bus"){
            $searchname = "output.bus";
        }else{
            die RED, "ERROR: Name of alignment tool is not provided or not valid for scasa.\n";
        }
        if($postalign_dir eq ""){
        my @current = `find $align_out -name "$searchname"`;
        chomp @current;
        s/(.*)_alignout\/.*/$1/ for @current;
        s/.*\/(.*)/$1/ for @current;
        @alignout_names = @current;
        }else{
            if($mapper eq "salmon alevin"){
                my @current = `find $postalign_dir -name "bfh.txt"`;
                chomp @current;
                s/(.*)\/alevin\/$searchname/$1/ for @current;
                s/.*\/(.*)/$1/ for @current;
                @alignout_names = @current;
            }elsif ($mapper eq "kallisto bus"){
                my @current = `find $postalign_dir -name "$searchname"`;
                chomp @current;
                s/(.*)\/$searchname/$1/ for @current;
                s/.*\/(.*)/$1/ for @current;
                @alignout_names = @current;
            }
        }
    }
	foreach my $i (@alignout_names){
	$current_out = "$quant_out/$i\_quant/";
        if (-d $current_out){
              print RED, "Directory $current_out already exists. Writing into existing directory..\n";
            }else{
              system("mkdir $current_out");
        }
	$quant_script = $log_dir."scasa.$project_name.$i.quant.sh";
	print RUNALL "sh $quant_script & echo \$\! >> $pid_path\n";
	open QUANT, ">$quant_script" or die "Unable to open $quant_all\n";
	print QUANT "#!/usr/bin/env bash\n";
	print QUANT "echo \"Begin Scasa quantification for sample $i..\"\n";
        if($mapper eq "salmon alevin"){
            $current_name = "$align_out/$i\_alignout/alevin/bfh.txt";
            $current_name =~ s/_alignout_alignout/_alignout/g;
	print QUANT "Rscript $scasa_quant_script11 $current_name $current_out $num_threads 1>$log_dir/$project_name.batch_quant.step1.o\n";
        }elsif ($mapper eq "kallisto bus"){
            $current_name = "$align_out/$i\_alignout/output.correct.sort.bus.txt";
            $current_name =~ s/_alignout_alignout/_alignout/g;
            # Current lib.size threshold is > 200 per cell
            print QUANT "Rscript $scasa_quant_script12 $current_name $current_out 200 1>$log_dir/$project_name.batch_quant.step1.o\n";
        }
    print QUANT "Rscript $scasa_quant_script2 in=$current_out design.matrix=$xmatrix Rsource=$scasa_Rsource core=$num_threads 1>$log_dir/$project_name.batch_quant.step2.o\n";
        if($mapper eq "salmon alevin"){
            print QUANT "Rscript $scasa_quant_script3 workdir=$current_out txgroup=$isoform_info_alevin 1>$log_dir/$project_name.batch_quant.step3.o\n";
        }elsif($mapper eq "kallisto bus"){
            print QUANT "Rscript $scasa_quant_script3 workdir=$current_out txgroup=$isoform_info_bustools 1>$log_dir/$project_name.batch_quant.step3.o\n";
        }
	close QUANT;
	}

	close RUNALL;

	$start_run = time();
	system("sh $quant_all");
	sleep 30;
	&checkPs("$pid_path");
	$end_run = time();
	$run_time = $end_run - $start_run;
    print "Congratulations! Scasa single cell RNA-Seq transcript quantification has completed in $run_time seconds!\n";
}


print "All done!\n";

sub checkPs{
	my $pid_path = shift;
	my $stop = 1;
	open PID, "<$pid_path" or die "Unable to open $pid_path\n";
	my @pid_list = <PID>;
	chomp @pid_list;
	my $k = scalar @pid_list;
	while ($stop){
		foreach my $i (@pid_list){
			my @pid_status = `ps -p $i`;
			chomp @pid_status;
			my $grp = grep(/$i/, @pid_status);
			if (grep(/$i/, @pid_status)){
				sleep 120;
			} else {
				$k = $k -1;
				@pid_list =  grep {!/$i/} @pid_list;
			}
		}
		if ($k == 0){
			$stop = 0;
		}
	}
	close PID;
}

__END__

=head1 NAME

=head1 SYNOPSIS

scasa [options] [argument]

For detailed information on the flags, visit: https://github.com/eudoraleer/scasa/wiki

=head1 OPTIONS

=over 20

=item B<--help,-h>

Print scasa options and their arguments and descriptions.

=item B<--mapper,-m>

Alignment tool <STRING, optional argument, default is set to salmon_alevin. Current scasa only supports these arguments to the option: salmon_alevin, kallisto_bus>

=item B<--project,-p>

Project name <STRING, optional, default = My_Project. No space is allowed, please use '-' or '_' or "." symbols to replace space. If you want to rerun existing folder, provide the timestamp suffix project folder name (not project directory)>

=item B<--in,-i>

Provide an input directory containing input FASTQ files<STRING, optional, no space in directory path is allowed, default is set to current directory>

=item B<--out,-o>

Output directory <STRING, optional, no space is allowed. Output will be generated under output directory with given project name with timestamp suffix, default is set to current directory>

=item B<--createxmatrix,-b>

YES/NO to generate a X-matrix <STRING, optional, provide YES or NO to the option. Default is set to NO>

=item B<--xmatrix,-x>

State Xmatrix reference file if --createxmatrix is set to NO <STRING, optional, two preset options: alevin or bustools. Give argument alevin to use Xmatrix for alevin alignment output or give bustools as argument to use Xmatrix for bustools alignment output data. no space in directory path is allowed. If option is not set, default is set to use scasa prebuilt Xmatrix for alevin alignment output data>

=item B<--ref,-r>

Provide a directory to reference fasta file <STRING, required, provide a fasta reference file, currently scasa only supports hg38. Users could download fasta reference via scasa Github. No default>

=item B<--index,-e>

YES/NO to run indexing for reference fasta file <STRING, optional, provide YES or NO to the option. Default is set to YES>

=item B<--index_dir,-d>

Provide a directory to reference fasta index <STRING, optional, this option is needed if --index option is set to NO. No space in directory path is allowed, no default>

=item B<--align,-a>

YES/NO to run pseudo-alignment step <STRING, optional, if set to YES, please state in the --mapper option which pseudo-alignment to use. Currently, scasa only supports two alignment tools: salmon alevin, kallisto bus. Default is set to YES>

=item B<--quant,-q>

YES/NO to run scasa quantification step to produce transcript counts <STRING, optional, default is set to YES>

=item B<--postalign_dir,-g>

Provide a directory to post-alignment files <STRING, optional, this option is only needed if --align is set to NO. No space in directory path is allowed, no default. Currently we only supports output from alevin or bustools>

=item B<--samplesheet,-s>

Provide a directory to a samplesheet file containing input FASTQ file names, one row one pair of fastq files, separated by a comma. See scasa Github for a sample of samplesheet file <STRING, optional, this option is for users with many fastq files to run, if this option is not used, please supply fastq names to option --fastq. If both --samplesheet and --fastq options are not provided, scasa will look for fastq files in the input directory supplied by user via option --in. No default for this option>

=item B<--fastq,-f>

Provide fastq file names, separate them by commas, make sure that you have R1 and R2 together the same prefix name within each pair of fastq files <STRING, optional, this option is for users with few fastq files to run, provide argument to either --samplesheet or --fastq but not both. If both --samplesheet and --fastq options are not provided, scasa will look for fastq files in the input directory supplied by user via option --in>

=item B<--tech,-t>

Provide the technlogy used for sequencing <STRING, optional, arguments to this option: 10xv1, 10xv2, 10xv3. Currently, scasa only supports sequencing output from 10X 3'  Chromium V1, V2 and V3 chemistries. Default is set to 10xv3>

=item B<--whitelist,-w>

Provide a directory to white list file for barcode correction. Visit scasa on Github for more information <STRING, optional, However, if --xmatrix is set to YES, this option is required for Xmatrix generation. No space in directory path is allowed, no default>

=item B<--cellthreshold,-c>

Provide a threshold for number of expected cells to be produced <NUMERIC, optional, this option is valid for alevin alignment step only, default is set to no expected cells.>

=item B<--nthreads,-n>

Number of threads to be used for running scasa <NUMERIC, optional, default is set to 4>

=back

=head1 DESCRIPTION

=cut
