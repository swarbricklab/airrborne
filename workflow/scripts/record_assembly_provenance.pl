#!/usr/bin/env perl
# Write the producer record only after this job has completed a new assembly.
# The host supplies its source/configuration context; the running image supplies
# its installed TRUST4 version, and the completed files supply their identities.
#
# Usage: record_assembly_provenance.pl --context JSON --reference FILE
#   --output FILE --airr FILE --cdr3 FILE --final FILE
#
# The run_trust4 rule runs it in the TRUST4 image as the job's last command.
# The record (schema_version 2.0) holds the library_id, the source, image and
# settings from the context, the read pairs, the measured TRUST4 version, and
# the name, SHA-256 and size of the reference and the three retained outputs.
# It is written to a temporary file and renamed into place.

use strict;
use warnings;
use Cwd qw(abs_path);
use Digest::SHA;
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempfile);
use Getopt::Long qw(GetOptions);
use JSON::PP;

my ($context_text, $reference, $output, $airr, $cdr3, $final);
GetOptions(
    'context=s' => \$context_text,
    'reference=s' => \$reference,
    'output=s' => \$output,
    'airr=s' => \$airr,
    'cdr3=s' => \$cdr3,
    'final=s' => \$final,
) or die "Invalid assembly provenance arguments\n";
die "Required: --context JSON --reference FILE --output FILE --airr FILE --cdr3 FILE --final FILE\n"
    if @ARGV || grep { !defined($_) || !length($_) } ($context_text, $reference, $output, $airr, $cdr3, $final);

my $json = JSON::PP->new->utf8->canonical->pretty;
my $context = $json->decode($context_text);
die "Assembly context must be a JSON object\n" unless ref($context) eq 'HASH';
my %allowed = map { $_ => 1 } qw(library_id source image settings read_pairs);
die "Assembly context contains undeclared fields\n" if grep { !$allowed{$_} } keys %{$context};
die "Assembly context needs a string library_id\n"
    unless is_string($context->{library_id}) && length($context->{library_id});
for my $field (qw(source image settings)) {
    die "Assembly context needs a $field object\n" unless ref($context->{$field}) eq 'HASH';
}
die "Assembly context needs a source revision and dirty state\n"
    unless is_string($context->{source}{revision}) && length($context->{source}{revision})
        && JSON::PP::is_bool($context->{source}{dirty});
die "Assembly context needs an immutable actual image identity\n"
    unless ($context->{image}{format} // '') =~ /^(?:oci|docker|sif)$/
        && ($context->{image}{digest} // '') =~ /^sha256:[0-9a-f]{64}$/
        && ($context->{image}{configured_reference} // '') =~ /\@sha256:[0-9a-f]{64}$/;
die "Assembly context needs explicit paired input filenames\n"
    unless ref($context->{read_pairs}) eq 'ARRAY' && @{$context->{read_pairs}};
for my $pair (@{$context->{read_pairs}}) {
    die "Every input pair needs only r1 and r2 filenames\n"
        unless ref($pair) eq 'HASH' && keys(%{$pair}) == 2
            && is_string($pair->{r1}) && length($pair->{r1})
            && is_string($pair->{r2}) && length($pair->{r2});
}

my %paths = (airr => $airr, cdr3 => $cdr3, final_out => $final);
my $output_absolute = abs_path($output) // File::Spec->rel2abs($output);
for my $path ($reference, values %paths) {
    my $absolute = abs_path($path);
    die "Assembly provenance input is missing: $path\n" unless defined($absolute) && -f $path;
    die "The producer sidecar must not replace an assembly input or output\n" if $absolute eq $output_absolute;
}

my $record = {
    schema_version => '2.0',
    stage => 'assembly',
    library_id => $context->{library_id},
    producer => {
        source => $context->{source},
        image => $context->{image},
        settings => $context->{settings},
        software => [{
            name => 'TRUST4',
            version => installed_trust4_version(),
            status => 'measured',
            method => 'installed run-trust4 version declaration',
        }],
    },
    inputs => {
        read_pairs => $context->{read_pairs},
        reference => file_identity($reference),
    },
    outputs => { map { $_ => file_identity($paths{$_}) } sort keys %paths },
};

# A partial write cannot become a producer record. The temporary file shares
# the destination directory, so rename commits the complete JSON in one step.
my $rendered = $json->encode($record);
my $directory = dirname($output);
make_path($directory) unless -d $directory;
my ($handle, $temporary) = tempfile('.assembly-provenance-XXXXXX', DIR => $directory, UNLINK => 0);
my $success = eval {
    binmode($handle, ':raw') or die "Cannot set raw output mode: $!\n";
    print {$handle} $rendered or die "Cannot write producer record: $!\n";
    close($handle) or die "Cannot close producer record: $!\n";
    chmod(0666 & ~umask(), $temporary) or die "Cannot set producer record permissions: $!\n";
    rename($temporary, $output) or die "Cannot publish producer record: $!\n";
    1;
};
if (!$success) {
    my $error = $@;
    unlink($temporary);
    die $error;
}

sub is_string {
    my ($value) = @_;
    return 0 if !defined($value) || ref($value);
    # JSON::PP preserves the distinction between "0460" and the number 460.
    return JSON::PP->new->allow_nonref->encode($value) =~ /^"/;
}

sub file_identity {
    my ($path) = @_;
    open(my $handle, '<:raw', $path) or die "Cannot read $path: $!\n";
    my $digest = Digest::SHA->new(256)->addfile($handle)->hexdigest;
    close($handle) or die "Cannot close $path: $!\n";
    return {name => basename($path), sha256 => $digest, size_bytes => -s $path};
}

sub installed_trust4_version {
    my ($launcher) = grep { -f $_ && -x $_ }
        map { File::Spec->catfile($_, 'run-trust4') } File::Spec->path();
    die "Cannot find installed run-trust4\n" unless defined($launcher);
    open(my $handle, '<', $launcher) or die "Cannot read installed run-trust4: $!\n";
    my @versions;
    while (my $line = <$handle>) {
        push(@versions, $1) if $line =~ /^\s*my\s+\$version\s*=\s*["']([^"']+)["']\s*;/;
    }
    close($handle) or die "Cannot close installed run-trust4: $!\n";
    die "Cannot identify one installed TRUST4 version declaration\n" unless @versions == 1;
    return $versions[0];
}
