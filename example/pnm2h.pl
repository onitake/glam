#!/usr/bin/perl

use strict;
use warnings;
use IO::File;

my $in = (defined($ARGV[0]) ? IO::File->new($ARGV[0], 'r') : \*STDIN) || die("Can't open input");
my $out = (defined($ARGV[1]) ? IO::File->new($ARGV[1], 'w') : \*STDOUT) || die("Can't open output");

sub readheader {
	my ($fd) = @_;
	my %header;
	while ((!defined($header{type}) || !defined($header{width}) || !defined($header{height}) || !defined($header{colors})) && defined(my $line = <$fd>)) {
		if ($line !~ /^#/) {
			if (!defined($header{type})) {
				$line =~ /P([0-9]+)/;
				if (defined($1)) {
					$header{type} = $1;
				} else {
					die("Invalid map type");
				}
			} elsif (!defined($header{width})) {
				$line =~ /([0-9]+)\s+([0-9]+)/;
				if (defined($1) && defined($2)) {
					$header{width} = $1;
					$header{height} = $2;
				} else {
					die("Invalid resolution line");
				}
			} elsif (!defined($header{colors})) {
				$line =~ /([0-9]+)/;
				if (defined($1)) {
					$header{colors} = $1;
				} else {
					die("Invalid color line");
				}
			}
		}
	}
	return %header;
}

sub readpixels {
	my ($fd, $width, $height) = @_;
	my @pixels;
	while (!eof($fd)) {
		read($fd, my $rgb, 3);
		my @pixel = unpack('C3', $rgb);
		push(@pixels, \@pixel);
	}
	return @pixels;
}

my %header = readheader($in);

if ($header{type} != 6) {
	die("Unsupported map type $header{type}");
}
if ($header{colors} != 255) {
	die("Unsupported number of bits per pixel");
}

my @data = readpixels($in, $header{width}, $header{height});

print($out "unsigned int image_width = $header{width};\n");
print($out "unsigned int image_height = $header{height};\n");
print($out "unsigned char image_samples[] = { ");
for my $rgb (@data) {
	printf($out "0x%02x,0x%02x,0x%02x,", @$rgb);
}
print($out "};\n");
