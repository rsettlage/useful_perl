#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Honngseok Tae
# Date: 2010
######################################

use strict;

if($#ARGV < 1){ print STDERR "Needs parameters\n  $0 <job ID> <command> <options for the command....>\n\n"; exit; }

my $jobID = shift;
my $fn = "qsub.$jobID.sh";
my $out;
open($out, "> $fn") || die "Could not create a script file : $!\n";

my $i = 0;
while($i == 0){
	$i = int(rand()*10%7+1); ### to avoid node 8
}

my $node = "node00$i.cm.cluster";

print $out "#!/bin/bash\n";
print $out "#PBS -N $jobID\n";
print $out "#PBS -l walltime=20:00:00 \n"; #maximum : 150:00:00
#print $out "#PBS -l pmem=8000mb\n";
#print $out "#PBS -l nodes=1:ppn=5\n"; ## any 1 node and 5 processes
#print $out "#PBS -l nodes=$node:ppn=5\n"; ## specific node and 5 processes
#print $out "#PBS -l nodes=$node\n"; ## specific node and 1 process
#print $out "#PBS -m ae\n#PBS -M htae\@vbi.vt.edu\n"; ## send an email when this job ends

print $out "echo start $jobID `date`\n";
my $query = "";
foreach my $str (@ARGV) {
	if($str =~ / /) { $query .= "\'$str\' ";}
	else { $query .= $str . " "; }
}
print $out "echo \"### $query\"\n$query\n";
print $out "echo end $jobID `date`\n";
close($out);

system("qsub -d $ENV{'PWD'} $fn");

unlink($fn);
