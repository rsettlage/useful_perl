#!/usr/bin/perl -w
use strict;


my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
my $year = 1900 + $yearOffset;
my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year"; 

open (OUT, ">syncdone.txt") || die "Couldn't open syncdone\n";
print OUT "Last synced at $theTime\n";
