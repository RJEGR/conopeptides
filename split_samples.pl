#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

# to run perl split_samples.pl
# Function to split file by n number of lines
sub split_file_by_lines {
    my ($file_path, $lines_per_chunk) = @_;
    
    open my $fh, '<', $file_path or die "Could not open file '$file_path': $!";
    
    my $file_count = 1;
    my $line_count = 0;
    my $output_file;
    my $output_fh;

    while (my $line = <$fh>) {
        if ($line_count % $lines_per_chunk == 0) {
            close $output_fh if defined $output_fh;
            my $output_file_name = sprintf("%s_part_%02d.txt", basename($file_path, ".txt"), $file_count);
            open $output_fh, '>', $output_file_name or die "Could not open file '$output_file_name': $!";
            $file_count++;
        }
        print $output_fh $line;
        $line_count++;
    }
    
    close $output_fh if defined $output_fh;
    close $fh;
}

# Read input_file and lines_per_chunk from standard input
print "Enter the input file path: ";
my $input_file = <STDIN>;
chomp($input_file);

print "Enter the number of lines per chunk: ";
my $lines_per_chunk = <STDIN>;
chomp($lines_per_chunk);

# Validate that lines_per_chunk is a positive integer
unless ($lines_per_chunk =~ /^\d+$/ && $lines_per_chunk > 0) {
    die "Invalid input for lines per chunk. Please enter a positive integer.\n";
}

# Split the input file by the specified number of lines
split_file_by_lines($input_file, $lines_per_chunk);

print "File split completed.\n";