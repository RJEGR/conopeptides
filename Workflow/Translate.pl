#!/usr/bin/perl

use strict;
use warnings;
use LWP::UserAgent;

# Create a user agent object
my $browser = LWP::UserAgent->new;

# Open the input FASTA file
my $input_file = 'input.fasta';  # Specify your input FASTA file here
open(my $fh, '<', $input_file) or die "Could not open file '$input_file' $!";

# Prepare to read sequences
my $header;
my $sequence = '';

while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ /^>(.*)/) {
        # If we encounter a header line, process the previous sequence (if any)
        if ($sequence) {
            # Send the nucleotide sequence for translation
            my $response = $browser->post(
                'https://web.expasy.org/cgi-bin/translate/dna2aa.cgi',
                [
                    'dna_sequence'    => $sequence,
                    'output_format'   => 'fasta'
                ]
            );

            if ($response->is_success) {
                print ">$header\n";  # Print header
                print $response->content;  # Print translated content
            } else {
                warn "Failed to translate sequence: " . $response->status_line;
            }
        }
        # Start a new sequence
        $header = substr($line, 1);  # Remove '>' from header
        $sequence = '';  # Reset sequence for new entry
    } else {
        # Append to the sequence
        $sequence .= uc($line);  # Convert to uppercase for consistency
    }
}

# Process the last sequence in the file (if any)
if ($sequence) {
    my $response = $browser->post(
        'https://web.expasy.org/cgi-bin/translate/dna2aa.cgi',
        [
            'dna_sequence'    => $sequence,
            'output_format'   => 'fasta'
        ]
    );

    if ($response->is_success) {
        print ">$header\n";  # Print header for last sequence
        print $response->content;  # Print translated content for last sequence
    } else {
        warn "Failed to translate sequence: " . $response->status_line;
    }
}

# Close file handle
close($fh);