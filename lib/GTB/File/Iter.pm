# $Id: Iter.pm 6569 2013-09-17 21:00:15Z pchines $
package GTB::File::Iter;
# Possible CPAN name: Iterator::File::Delimited?

use strict;
use Carp qw(carp croak);
use GTB::File qw(Open);

our $VERSION = '0.05';

sub new {
    my ($pkg, %param) = @_;
    if (!$param{fh} && !$param{file}) {
        croak "GTB::File::Iter::new: either fh or file is required";
    }
    my $self = bless {}, ref $pkg || $pkg;
    for my $f (qw(file fh fs rs cols has_header retain_skipped skip)) {
        $self->{$f} = delete $param{$f};
    }
    if (keys %param) {
        carp "GTB::File::Iter::new: unknown parameters "
            . join (" ,", map { "'$_'" } sort keys %param);
    }
    if (!$self->{fh}) {
        $self->{fh} = Open($self->{file});
    }
    $self->{fs} ||= "\t";
    $self->{rs} ||= "\n";
    if (!defined $self->{has_header}) {
        $self->{has_header} = $self->{cols} ? 0 : 1;
    }
    if ($self->{has_header}) {
        my $line = $self->_line();
        $self->{cols} ||= [split $self->{fs}, $line];
    }
    return $self;
}

sub skipped_lines {
    my ($self) = @_;
    return delete($self->{_skipped}) || '';
}

sub _line {
    my ($self) = @_;
    return if $self->{_is_done};
    my $line;
    my $fh = $self->{fh};
    local $/ = $self->{rs};
    while ($line = <$fh>) {
        if ($self->{skip}) {
            my $skip = 0;
            if (UNIVERSAL::isa($self->{skip}, 'CODE')) {
                $skip = $self->{skip}->($line);
            }
            else {
                $skip = ($line =~ /$self->{skip}/);
            }
            if ($skip) {
                if ($self->{retain_skipped}) {
                    $self->{_skipped} .= $line;
                }
                next;
            }
        }
        last;
    }
    if (!defined $line) {
        $self->{_is_done} = 1;
    }
    else {
        chomp $line;
    }
    return $line;
}

sub columns {
    my ($self) = @_;
    if (!$self->{cols}) {
        if (!$self->{_line}) {
            # This is an error because reading the line to find out would
            # potentially skip the line; this could be remedied by either
            # storing the read line in a buffer, or by reopening the file.
            # The workaround is for the caller to call some other method
            # (next_line, current_line, etc.) before calling this method.
            croak "columns: called before number of columns is known. "
                . "Please call\na retrieval method (e.g. next_line, "
                . "current_line, etc. before calling columns()";
        }
        my $ra = $self->current_array();
        my @cols = ('A'..'Z');
        my $prefix = 'A';
        while (@$ra > @cols) {
            push @cols, ("${prefix}A".."${prefix}Z");
            ++$prefix;
        }
        splice @cols, @$ra;
        $self->{cols} = \@cols;
    }
    return wantarray ? @{ $self->{cols} } : $self->{cols};
}

sub current_array {
    my ($self) = @_;
    if (!$self->{_array}) {
        if (!defined $self->{_line}) {
            return if !$self->next_line();
        }
        $self->{_array} = [ split $self->{fs}, $self->{_line} ];
    }
    return wantarray ? @{ $self->{_array} } : $self->{_array};
}

sub current_hash {
    my ($self) = @_;
    if (!$self->{_hash}) {
        my $ra_data = $self->current_array();
        return if !$ra_data;
        my $ra_cols = $self->columns();
        # TODO: test for same size array?
        my %h;
        @h{@$ra_cols} = @$ra_data;
        $self->{_hash} = \%h;
    }
    return wantarray ? %{ $self->{_hash} } : $self->{_hash};
}

sub current_line {
    my ($self) = @_;
    if (!$self->{_line}) {
        return if !$self->next_line();
    }
    return $self->{_line} . $self->{rs};
}

sub next_array {
    my ($self) = @_;
    $self->next_line();
    return $self->current_array();
}

sub next_hash {
    my ($self) = @_;
    $self->next_line();
    return $self->current_hash();
}

sub next_line {
    my ($self) = @_;
    $self->{_array} = $self->{_hash} = undef;
    $self->{_line} = $self->_line();
    return defined $self->{_line} ? $self->{_line} . $self->{rs}
                                  : undef;
}

sub is_done {
    my ($self) = @_;
    return $self->{_is_done} ? 1 : 0;
}

sub reopen {
    my ($self, $file) = @_;
    if ($file) {
        $self->{file} = $file;
    }
    else {
        $file = $self->{file};
    }
    if ($file) {
        $self->{fh} = Open($file);
    }
    else {
        seek $self->{fh}, 0, 0
            or die "Failed to seek to beginning of stream, $!\n";
    }
    delete $self->{skipped};
    delete $self->{_is_done};
    $self->{_hash} = $self->{_array} = $self->{_line} = undef;
    if ($self->{has_header}) {
        my $line = $self->_line();
        $self->{cols} = [split $self->{fs}, $line];
    }
}

1;
__END__

=head1 NAME

GTB::File::Iter - read records from file

=head1 SYNOPSIS

    # simplest case: tab-delimited file with header line:

    use GTB::File::Iter;
    my $iter = GTB::File::Iter->new(file => '/path/to/my.txt');
    while (my $rh = $iter->next_hash()) {
        # work with hashref
    }

=head1 DESCRIPTION

Iterator for delimited files.

=head1 METHODS

=head2 new

Create new GTB::File::Iter object.

  Usage: my $iter = GTB::File::Iter->new(%params);
   Args: file   => path to filename (require file or fh)
         fh     => already opened filehandle
         fs     => field separator (default is tab)
         rs     => record separator (default is newline)
         cols   => arrayref of column names
     has_header => whether file has a header line, 1/0
 retain_skipped => whether object should remember "skipped"
                    lines for later use, 1/0; if set, these
                    lines are accessible via skipped_lines()
         skip   => regexp or subroutine ref that returns 1/0
                    will ignore lines where value is 1
 Return: iterator object

If 'cols' are provided, the default is to assume there is no header line,
otherwise, the first line of the file that is not skipped is assumed to be
a header with column names.  To change this behavior, you may explicitly set
'has_header'.  If has_header==1 and cols are provided, the first line is
ignored.  If has_header==0 and no cols are provided, column names are named
like Excel columns, i.e. A-Z, AA-AZ, BA-BZ, etc.

Note: For sorted files where you wish to read multiple files in a
synchronized manner, or to advance to a particular record, use the
GTB::File::SortedIter instead; you will need to define a 'key' that
corresponds to the sort order of the file.  See L<GTB::File::SortedIter>.

=head2 columns

  Usage: my @cols = $iter->columns();
 Return: array of column names

=head2 current_array

Return array or arrayref for current_line.

=head2 current_hash

Return hash or hashref for current_line.

=head2 current_line

Returns current line from the file, without advancing.  This is the same as
the value returned in the preceding call to next_line.  If no prior call to
next_line was made, current_line behaves like next_line, returning the first
line of the file.

=head2 next_array

Return array or array reference with fields from the next line.  Fields are
split on the "field separator".

=head2 next_hash

Return hash or hash reference with fields from the next line in the file.
Fields are split on the "field separator", and assigned values based on the
list of columns (either provided, or read from a header line).

=head2 next_line

Returns next line from file, undef when all lines have been read.

=head2 is_done

Return true if done reading this file (at EOF).

=head2 reopen

Rewind file to beginning, and reinitialize.

=head2 skipped_lines

If 'retain_skipped' option is set, object will remember skipped lines, and
this method will return a single string containing all lines skipped (so
far).  Assuming all of the skipped lines are near the beginning of
the file, it is often best to call this method after a call to columns() or
one of the data access methods.  Each call to skipped_lines will clear the
buffer, so that the next call to skipped lines only includes lines skipped
since the previous call.

=head1 AUTHORS

 Peter Chines <pchines@mail.nih.gov>

=head1 SEE ALSO

L<GTB::File::SortedIter>, L<GTB::File>

=head1 LEGAL

This software is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S. Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation. 

=cut
