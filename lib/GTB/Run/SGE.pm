# $Id: SGE.pm 6587 2013-10-29 20:58:48Z pchines $
package GTB::Run::SGE;

use strict;
use Carp qw(carp croak confess);
use File::Temp qw(tempfile);
use File::Spec;
use GTB::Run qw(get_output as_number);
use GTB::File qw(Open);
# requires XML::Twig, if qstat() method is used.

our $VERSION = '0.12';
our $DEFAULT_MEMORY_LIMIT = 1024 * 1024 * 1024; # 1Gb
our %QstatFieldMap = (
        JB_job_number   => 'job_number',
        JB_name         => 'job_name',
        JB_owner        => 'user',
        JB_project      => 'project',
        );
our $QACCT_FILE_GLOB = $ENV{SGE_ROOT} ?
    "$ENV{SGE_ROOT}/$ENV{SGE_CELL}/common/accounting*" : "/dev/null";
our @QacctFields = qw( queue_name host_name group user job_name job_number
     account priority submission_t start_t end_t
     failed exit_status wallclock utime stime
     ru_maxrss ru_ixrss ru_ismrss ru_idrss ru_isrss ru_minflt ru_majflt
     ru_nswap ru_inblock ru_oublock ru_msgsnd ru_msgrcv ru_nsignals
     ru_nvcsw ru_nivcsw
     project department granted_pe slots task_id
     cpu mem io category iow pe_taskid max_vmem
     arid ar_submission_time );
our %QacctFieldMap = (
        qname           => 'queue_name',
        hostname        => 'host_name',
        owner           => 'user',
        jobname         => 'job_name',
        jobnumber       => 'job_number',
        taskid          => 'task_id',
        ru_wallclock    => 'wallclock',
        maxvmem         => 'max_vmem',
        );

my $EMPTY = q{};
my $SPACE = q{ };
my $COMMA = q{,};
our $RE_NAME = qr/^([\w.-]+)$/;     # single name/id, no wildcards
our $QSTAT_MAX_WAIT = 2;            # seconds to wait (retry once per sec)
our $QACCT_MAX_WAIT = 12;           # seconds total waiting for qacct jobs

sub new {
    my ($pkg, %param) = @_;
    my $self = bless { %param }, ref $pkg || $pkg;
    return $self;
}

# Handle parameters for OO or non-OO calls, default params
# should be in base class?
sub _params {
    my $rhp = shift;
    my $self = __PACKAGE__;
    if (ref $_[0] eq __PACKAGE__ || $_[0] eq __PACKAGE__) {
        $self = shift;
    }
    my %p;
    if (!$rhp->{defaults}) {
        # no defaults
    }
    elsif ($rhp->{defaults} eq 'self') {
        if (ref $self) {
            %p = %$self;
        }
    }
    elsif (ref $rhp->{defaults} eq 'HASH') {
        %p = %{ $rhp->{defaults} };
    }
    else {
        confess("defaults should be 'self' or reference to hash");
    }

    if (@_ == 1 && $rhp->{single}) {
        $p{$rhp->{single}} = shift;
    }
    elsif (@_ % 2 == 0) {
        my %n = @_;
        for my $k (keys %n) {
            $p{$k} = $n{$k};
        }
    }
    else {
        my $from_sub = (caller(1))[3];
        confess("$from_sub: odd number of params (expected hash)");
    }
    # TODO: add support for $rhp->{required} (array ref of required params)
    # TODO: add support for $rhp->{strict} (warn about unspecified params,
    #   i.e. not required, defaults, or this list, if a list is given)
    return ($self, %p);
}

sub submit_job {
    my ($self, %p) = _params({single => 'cmd', defaults => 'self'}, @_);
    my $cmd  = delete $p{cmd} or confess "submit_job: 'cmd' required";
    # automatically generate a script if redirection is employed
    # or if umask or preserve_umask is set
    if (defined $p{'umask'} || $p{preserve_umask} || $p{strict_mem}
            || $cmd =~ /[>|<]/ ) {
        return $self->submit_shell_script(script => $cmd, %p);
    }
    _set_mem_free(\%p);
    my $sge_cmd = _make_sge_cmd(\%p, $cmd);
    if ($p{dryrun}) {
        warn "Would have run: $sge_cmd\n";
        return 0;
    }
    my ($jid,$err,$rv) = get_output($sge_cmd);
    if ($err || $rv) {
        if ($p{sge_cmd} && $p{sge_cmd} eq 'qrsh' && $cmd =~ /^\S*\btrue$/) {
            if ($rv == 129 * 256) {
                warn "Encountered spurious SGE error 129 "
                    . "(HUP after job completion)\n";
            }
            elsif ($rv) {
                warn "New (probably harmless) SGE error code while "
                    . "running qrsh hold: $rv\n";
            }
            elsif ($err !~ /X11 forwarding/ and	$err !~ /Permanently added.*\(RSA\) to the list of known hosts/ and
		   $err !~ /Connection to.* closed\./ ) {
                warn "New (but harmless) warning while running qrsh hold: "
                    . "$err\n";
            }
        }
        else {
            die "Error running $sge_cmd\nError is: $jid $err $rv\n";
        }
    }
    if ($jid && $jid =~ /^(\d+)(?:\.\d+-\d+:\d+)?$/) {
        $jid = $1;
        if (ref $self) {
            push @{ $self->{_jobs} }, $jid;
        }
        return $jid;
    }
    elsif ($jid) {
        warn "Job ID '$jid' is not recognized";
    }
    return 0;
}

sub _make_sge_cmd {
    my ($rh_p, $cmd) = @_;
    my $job  = $rh_p->{sge_cmd}  || 'qsub';
    my %opts = %{ $rh_p->{sge_opts} || {} };
    if ($job eq 'qsub') {
        $job .= " -terse";
    }
    elsif ($job eq 'qrsh') {
        $job .= " -now no";
    }
    else {
        confess "submit_job: sge_cmd '$job' not recognized; "
            . "should be 'qsub' or 'qrsh'";
    }
    my @cmd = split " ", $cmd, 2;  # split cmd from args on any whitespace
    if ($cmd[0] !~ m{/}) {
        my ($fp,$err,$rv) = get_output("which $cmd[0]");
        if ($err) {
            croak "submit_job: can't find $cmd[0].\n$err";
        }
        if ($fp) {
            chomp $fp;
            $cmd[0] = $fp;
            $cmd = join $SPACE, @cmd;
        }
    }
    if (!$opts{'-b'}) {
        my $ft = `file $cmd[0]`;
        if ($ft && $ft !~ /ERROR|script|text/) {
            $opts{'-b'} = 'y';
        }
    }
    if ($rh_p->{hold_jid}) {
        my $hold;
        if (ref $rh_p->{hold_jid}) {
            $hold = join($COMMA, grep { $_ } @{ $rh_p->{hold_jid} });
        }
        else {
            $hold = $rh_p->{hold_jid};
        }
        if ($hold) {
            $job .= " -hold_jid $hold";
        }
    }
    if (!$opts{'-N'}) {
        if ($rh_p->{job_name}) {
            $opts{'-N'} = $rh_p->{job_name};
        }
        elsif ($cmd[0] =~ m{(?:^|/)([^/]+)$}) {
            $opts{'-N'} = $1;
        }
    }
    if ($rh_p->{job_prefix}) {
        $opts{'-N'} = $rh_p->{job_prefix} . $opts{'-N'};
    }
    if ($opts{'-N'} =~ m{^\d}) {
        $opts{'-N'} = "_$opts{'-N'}";
    }
    my $default_opts = $SPACE . ($ENV{SGE_OPT} || $EMPTY);
    $job .= join($SPACE, $default_opts, map { ($_, $opts{$_}) } keys %opts);
    return "$job $cmd";
}

sub _set_mem_free {
    my ($rh_p) = @_;
    if ($rh_p->{memory}) {
        my $mem_opts = "mem_free=$rh_p->{memory}";
        if ($rh_p->{strict_mem}) {
            $mem_opts .= ",h_vmem=$rh_p->{memory}";
        }
        my $pmem = as_number($rh_p->{memory});
        if ($rh_p->{sge_opts}{'-l'}) {
            if ($rh_p->{sge_opts}{'-l'} =~ /\b(?:mf|mem_free)=([\w.]+)/) {
                my $lmem = as_number($1);
                if ($pmem > $lmem) {
                    $rh_p->{sge_opts}{'-l'} =~ s/,?h_vmem=(?:[\w.]+)//;
                    $rh_p->{sge_opts}{'-l'} =~ s/^,//;
                    $rh_p->{sge_opts}{'-l'} =~
                        s/\b(?:mf|mem_free)=[\w.]+/$mem_opts/;
                }
            }
            else {
                $rh_p->{sge_opts}{'-l'} .= ",$mem_opts";
            }
        }
        else {
            $rh_p->{sge_opts}{'-l'} = $mem_opts;
        }
        delete $rh_p->{memory}; # no longer needed, avoid resetting again
    }
}

sub submit_shell_script {
    my ($self, %p) = _params({single => 'script', defaults => 'self'}, @_);
    if (!$p{script}) {
        croak "submit_shell_script: 'script' required";
    }
    if (!$p{job_name}) {
        if ($p{script} =~ m{^\s*\S*?([\w.-]+)}) {
            $p{job_name} = $1;
        }
    }
    my @template;
    if ($p{job_prefix}) {
        @template = ("$p{job_prefix}XXXXX");
    }
    # need to specify directory, since default TMPDIR may not be the
    # same place on a node
    push @template, DIR => File::Spec->rel2abs($p{tempdir} || '.');
    $p{script} =~ s/\s+$//;
    if (defined $p{'umask'} || $p{preserve_umask}) {
        $p{preserve_umask} = 0;
        my $umask = defined $p{'umask'} ? $p{'umask'} : umask;
        $p{'umask'} = undef;
        # handle non-executable shell scripts
        my @cmd = split " ", $p{script};
        if (-f $cmd[0] && !-x $cmd[0]) {
            my $fh = Open($cmd[0]);
            my $shebang = <$fh>;
            if ($shebang !~ s/^#!//) {
                $shebang = "source";
            }
            $p{script} = "$shebang " . $p{script};
        }
        $p{script} = sprintf("umask %o\n", $umask) . $p{script};
    }
    _set_mem_free(\%p);
    if ($p{strict_mem}) {
        # TODO: remove this, when make strict_mem the default
        # this is belt-and-suspenders territory; ideally the h_vmem SGE 
        # resource should take care of this, and in a nicer way.
        my $kmem;
        if ($p{sge_opts}{'-l'} =~ /\b(?:mf|mem_free)=([\w.]+)/) {
            $kmem = int(as_number($1) / 1024)+1;
        }
        else {
            $kmem = $DEFAULT_MEMORY_LIMIT / 1024;
        }
        $p{script} = "ulimit -v $kmem\n" . $p{script};
        $p{strict_mem} = 0; # avoid endless loop
    }
    if ($p{dryrun}) {
        my $temp = "$template[2]/$template[0]";
        my $job = _make_sge_cmd(\%p,$temp);
        warn "Dry-run: would have written temporary script $temp\n"
            . "...with contents:\n$p{script} && rm $temp\n"
            . "...and executed via:\n$job\n";
        return 0;
    }
    my ($tfh, $temp) = tempfile(@template);
    print $tfh "#!/bin/bash\n";
    print $tfh "$p{script} && rm $temp\n";
    close $tfh or die "Can't write tempfile $temp, $!\n";
    system "chmod +x $temp" and die "Error making $temp executable, $!\n";
    delete $p{script};
    return $self->submit_job(%p, cmd => $temp);
}

sub wait {
    my ($self, @par) = @_;
    my $ra_jobs;
    if (@par == 0) {    # wait for all jobs to complete
        $ra_jobs = [ $self->job_ids_submitted() ];
    }
    elsif (@par == 1) { # wait for single job
        if (ref $par[0]) {
            croak "wait: when called with single argument, must be single "
                . "job ID";
        }
        $ra_jobs = [$par[0]];
    }
    elsif (@par > 1) {  # wait for specified jobs
        my %p = @par;
        if (!$p{jobs} || ref($p{jobs}) ne 'ARRAY') {
            croak "wait: when called with multiple arguments, must be hash "
                . "with key 'jobs' pointing to reference to array of job IDs";
        }
        $ra_jobs = [@{ $p{jobs} }];
    }
    if (@$ra_jobs) {
        submit_job(
            cmd      => '/bin/true',
            hold_jid => $ra_jobs,
            sge_cmd  => 'qrsh',
            sge_opts => { -b => 'y'},
            preserve_umask => 0,
            );
    }
    return 1;
}

sub job_ids_submitted {
    my ($self) = @_;
    if (ref $self) {
        return @{ $self->{_jobs}||[] };
    }
    return;
}

sub clear_job_list {
    my ($self) = @_;
    if (ref $self) {
        $self->{_jobs} = [];
    }
}

sub qstat {
    my ($self, %param) = @_;
    if (!defined $param{max_wait}) {
        $param{max_wait} = $QSTAT_MAX_WAIT;
    }
    my @qstat;
    my $cmd = 'qstat -xml -ext';
    if ($param{user} && $param{user} =~ $RE_NAME) {
        $cmd .= " -u $1";
    }
    elsif ($param{user}) {
        croak "qstat: '$param{user}' is not a valid user name";
    }
  {   # bare block works like a loop for redo
    my ($output, $err, $rc) = get_output($cmd);
    if ($rc) {
        warn "qstat failed with error code $rc, $err\n";
    }
    else {
        require XML::Twig;
        my $t = XML::Twig->new( twig_handlers => {
                'job_list' => sub { __parse_qstat_twig($_, \@qstat) },
                })->parse($output);
        if ($param{job} && $param{job} =~ /^\d+\n?$/) {
            @qstat = grep { $_->{job_number} == $param{job} } @qstat;
        }
        elsif (ref $param{job}) {
            my %sjob;
            @sjob{ @{ $param{job} } } = ();
            @qstat = grep { exists $sjob{$_->{job_number}}
                                || exists $sjob{$_->{job_name}} } @qstat;
        }
        elsif ($param{job}) {
            @qstat = grep { $_->{job_name} eq $param{job} } @qstat;
        }
        if ($param{submitted}) {
            my %sjob;
            @sjob{ $self->job_ids_submitted() } = ();
            @qstat = grep { exists $sjob{ $_->{job_number} } } @qstat;
        }
        if ($param{state}) {
            @qstat = grep { $_->{state} eq $param{state} } @qstat;
        }
        if (defined $param{project}) {
            @qstat = grep { $_->{project} eq $param{project} } @qstat;
        }
    }
    if (!@qstat && $param{max_wait}) {
        --$param{max_wait};
        sleep 1;
        redo;
    }
  }
    return scalar(@qstat) ? \@qstat : undef;
}

sub __parse_qstat_twig {
    my ($twig, $ra) = @_;
    my %q;
    for my $child ($twig->children()) {
        my $tag = $child->name;
        if ($QstatFieldMap{$tag}) {
            $tag = $QstatFieldMap{$tag};
        }
        my $val = $child->text_only;
        $q{$tag} = $val;
    }
    push @$ra, \%q;
}

sub qacct {
    my ($self, %param) = @_;
    if (!defined $param{max_wait}) {
        $param{max_wait} = $QACCT_MAX_WAIT;
    }
    my (@qacct, %jids, %jnames);
    if ($param{job}) {
        if (ref $param{job}) {
            my @errs;
            for my $job (@{ $param{job} }) {
                if ($job =~ /^\d+$/) {
                    $jids{$job} = undef;
                }
                elsif ($job =~ $RE_NAME) {
                    $jnames{$1} = undef;
                }
                else {
                    push @errs, $job;
                }
            }
            if (@errs) {
                croak "qacct: " . join($COMMA, map { "'$_'" } @errs)
                    . " do not look like valid job ids";
            }
        }
        elsif ($param{job} =~ /^\d+$/) {
            $jids{ $param{job} } = undef;
        }
        elsif ($param{job} =~ $RE_NAME) {
            $jnames{$1} = undef;
        }
        else {
            croak "qacct: '$param{job}' does not look like a valid job "
                . "name or number";
        }
    }
    if ($param{submitted}) {
        @jids{ $self->job_ids_submitted() } = undef;
    }
    my $min_jid_requested = 0;
    if (!keys %jnames) {
        if (keys %jids) {
            ($min_jid_requested) = sort { $a <=> $b } keys %jids;
        }
        else {
            croak "qacct: job name/number is required (for now)";
        }
    }
    my @qacct_files = sort glob $QACCT_FILE_GLOB;
    my $t0 = time;
    for my $ac_file (@qacct_files) {
        my $max_jid = 0;
        push @qacct, __parse_qacct_file(
                file       => $ac_file,
                max        => \$max_jid,
                job_number => \%jids,
                job_name   => \%jnames,
                );
        if ($min_jid_requested > $max_jid) {
            last;   # no need to look in older files
        }
    }
    my @missed_jids   = grep { !$jids{$_} } keys %jids;
    my @missed_jnames = grep { !$jnames{$_} } keys %jnames;
    if (@missed_jids || @missed_jnames || !@qacct) {
        %jids = ();
        %jnames = ();
        @jids{@missed_jids} = undef;
        @jnames{@missed_jnames} = undef;
        my $wait = $param{max_wait} - time + $t0;
        if ($wait > 0) {
            sleep $wait;
        }
        push @qacct, __parse_qacct_file(
                file       => $qacct_files[0],
                job_number => \%jids,
                job_name   => \%jnames,
                );
    }
    if (defined $param{failed}) {
        my @keep;
        for my $job (@qacct) {
            my $f = $job->{failed} || $job->{exit_status};
            if (($param{failed} && $f) || (!$param{failed} && !$f)) {
                push @keep, $job;
            }
        }
        @qacct = @keep;
    }
    return scalar(@qacct) ? \@qacct : undef;
}

sub __parse_qacct_file {
    my (%p) = @_;
    my @matched;
    my $fh = Open($p{file});
    while (<$fh>) {
        next if /^#/;
        chomp;
        my %j;
        @j{ @QacctFields } = split /:/;
        my $match = 0;
        for my $f (qw(job_number job_name)) {
            if ($p{$f} && $j{$f} && exists $p{$f}{$j{$f}}) {
                $p{$f}{$j{$f}} = $j{job_number};
                $match = 1;
            }
        }
        if (exists $p{max} && ${$p{max}} < $j{job_number}) {
            ${$p{max}} = $j{job_number};
        }
        if ($match) {
            push @matched, \%j;
        }
    }
    return @matched;
}

1;
__END__

=head1 NAME

GTB::Run::SGE - object for submitting SGE jobs

=head1 SYNOPSIS

    use GTB::Run::SGE;
    my $sge = GTB::Run::SGE->new();
    my $jid1 = $sge->submit_job($cmd1);
    my $jid2 = $sge->submit_job(cmd => $cmd2, hold_jid => $jid1);
    my $jid3 = $sge->submit_shell_script( script => $many_cmds,
            hold_jid => $jid1 );
    $sge->wait();

=head1 DESCRIPTION

This is an object for submitting SGE jobs.  It is intended to be far simpler
than the alternative ways of submitting and waiting for the completion of SGE
jobs.  It allows for submitting only batch jobs (no interactive sessions),
and allows you to specify dependencies among those jobs.

This module simplifies using SGE by supplying full paths for commands where
the path is not included, preserving your umask across SGE calls, and
automatically deciding whether the C<-b> (binary executable) option needs to
be passed to SGE commands.  It also simplifies access to qstat (active SGE
jobs) and qacct (completed SGE jobs, including failures) information.

New in version 0.12: the SGE_OPT environment variable may be used to pass
certain options through to the qsub or qrsh command.  These options will be
overridden by explicit option settings at the object or method level.

=head1 METHODS

=head2 new

Creates a new GTB::Run::SGE object.
    
    $sge = GTB::Run::SGE->new();

You may set default values by passing a hash to this constructor method.  Any
of the parameters used by submit_job() may have default values supplied via
this method, however the following are of particular interest:

=over 8

=item dryrun

If "dryrun" is given a true value, this SGE object will report the
commands that would have been submitted to STDERR, but will not actually
submit them.

=item job_prefix

This option adds a prefix to every job name that is submitted via this
object.  This can be useful for multi-step processes.

=item tempdir

tempdir, is a location that is reachable by the same path on the head
(submitting) node and all of the compute nodes.  Note that this temporary
directory is only used to write temporary script files to be submitted to
SGE, and does not affect the temporary directories used by the running
programs themselves.  By default, this module uses the current directory,
which should work as long as you have write permissions, and are not
in a local scratch or ramdisk directory.

=item preserve_umask

If preserve_umask is set to a true value, your current umask will be
propagated to the SGE job, rather than being set to an arbitrary value
(0022, in the case of the trek cluster).  The upside of setting this value
to true is that in group-shared directories, other members of your group will
be more easily able to work with the files and directories created by your
job.  The downside is that this requires creating a temporary script for
every command submitted.

=item umask

This option allows you to specify an arbitrary umask to apply to all files
and directories created by your SGE job.  Setting this to a value like 2 will
ensure that all files and directories are group-writable (unless they are
otherwise restricted).  Note that the value of this option must be a number,
so if you want to specify permissions using the usual octal representation,
you should prefix the octal number with a zero, e.g. '027' to restrict group
members from writing and deny all access to non-group members.

=item strict_mem

When set to a true value, each process will be stictly limited to the amount
of memory requested using the 'memory' parameter or the SGE "-l mem_free=XX"
resource request.  Jobs that exceed this amount of virtual memory will be
killed.  Note that while this is currently an option that you must turn on,
it will likely become the default in the near future.

WARNING: when a job is killed, this results in an "Out of Memory" error
message and a job exit_code of 1, but it is the user's responsibility to
check for these conditions.  Dependent jobs will proceed as if the job
succeeded. (see TODO's)

=back

=head2 submit_job

Submits an SGE job.  Must be a full command that can be run at the
commandline, and not include "qsub" or its parameters.  Certain parameters
can be set by passing hash.

  Usage: my $jid = $sge->submit_job($cmd);
     or: my $jid = $sge->submit_job(%job_params);
   Args: cmd        => commandline; required
         hold_jid   => job(s) this jobs depends on; may be single JID
                       or reference to array of JIDs
         sge_cmd    => SGE command (qsub, qrsh); default is qsub
         sge_opts   => hash ref of SGE options to be used; initial dash
                       is required on keys, e.g. { '-o' => '/dev/null' }
         job_name   => name of job in SGE queue (will not override -N
                       option in sge_opts)
         tempdir    => temporary directory that points to same place on
                       head and nodes (default is current directory)
     preserve_umask => flag that tells whether to use your current umask
                       value when running the job remotely via SGE;
                       default is true, but may be changed when creating
                       SGE object [see new()].
         umask      => specific umask value to use; must be numeric,
                       not a string (see perldoc -f umask)
         memory     => synonym for sge_opts { -l => mem_free=X }
         strict_mem => enforce memory limit: kill jobs that use more than
                       requested amount
 Return: job identifer

Note: be careful about quotes and escaping variables, especially those used
in awk or perl; while submit_job() does not expand these variables, the shell
does.  In most cases, it probably would be better to submit such commands
using submit_shell_script(), even if they consist of only a single command.
Also note that GTB::Run::SGE will "promote" commands that use any form of
redirection symbols to submit_shell_script(), even if these are quoted.
Finally, beginning with version 0.09 of this module, unless preserve_umask
is set to false, every command will be submitted using submit_shell_script(),
and will thus create a temporary file.

Starting with version 0.12, the contents of the SGE_OPT environment variable
will be used as defaults for all SGE qsub or qrsh commands.

=head2 submit_shell_script

Creates a temporary shell script file and submits it as an SGE job.  This is
useful when the command that you want to run includes pipelines or shell
redirects.  The temporary script is deleted as the last step in the process.

  Usage: my $jid = $sge->submit_shell_script($script);
     or: my $jid = $sge->submit_job(%job_params);
   Args: script     => shell script text; required
         hold_jid   => job(s) this jobs depends on; may be single JID
                       or reference to array of JIDs
         sge_cmd    => SGE command (qsub, qrsh); default is qsub
         sge_opts   => hash ref of SGE options to be used; initial dash
                       is required on keys, e.g. { '-o' => '/dev/null' }
         job_name   => name of job in SGE queue (will not override -N
                       option in sge_opts)
         tempdir    => temporary directory that points to same place on
                       head and nodes (default is current directory)
     preserve_umask => flag that tells whether to use your current umask
                       value when running the job remotely via SGE
         umask      => specific umask value to use; must be numeric,
                       not a string (see perldoc -f umask)
         memory     => synonym for sge_opts { -l => mem_free=X }
         strict_mem => enforce memory limit: kill jobs that use more than
                       requested amount
 Return: job identifer

(Due to implementation details, if the last line in your script is a comment,
the temporary script will not be deleted.  This is a "feature" that should
not be relied upon.)

Starting with version 0.12, the contents of the SGE_OPT environment variable
will be used as defaults for all SGE qsub or qrsh commands.

=head2 wait

Waits for all jobs (or specific jobs, if job IDs are provided) to complete
before returning.

  Usage: $sge->wait();
     or: $sge->wait($jid);
     or: $sge->wait(jobs => \@jids);
 Return: true

=head2 job_ids_submitted

Returns a list of job ID for all SGE jobs submitted through this object
instance.  Note that jobs may or may not be running.

  Usage: my @jids = $sge->job_ids_submitted();

=head2 clear_job_list

Resets list of job IDs, so that job_ids_submitted() returns an empty list,
until more jobs are submitted.

=head2 qstat

Get status of jobs in SGE queue.  Returns a reference to an array of hashes,
each of which represents a job, and has (at least) the following keys:

   job_number
   job_name
   user
   project
   queue_name
   slots
   state

 Usage: $ra = $sge->qstat();                # all queued jobs
        $ra = $sge->qstat(submitted => 1);  # submitted by this object
        $ra = $sge->qstat(user => 'username');    # user's jobs
        $ra = $sge->qstat(job  => 'job_name');    # job id or name
        $ra = $sge->qstat(job  => \@jids);        # multiple jobs
        $ra = $sge->qstat(project => 'prj_name'); # jobs for project
        $ra = $sge->qstat(state   => 'r');        # with given state code
    Args: may add max_wait => $secs to retry each second for specified
          time period; only retries if no jobs would be returned.  Can
          be useful when looking for a specific job, since it takes a
          few seconds for each job to be entered into SGE database.
          Defaults to two seconds; use max_wait => 0 to disable.  Whole,
          positive integers only.

 Returns: arrayref of hashes.  undef if no jobs match criteria.

Criteria may be combined.  Currently there is no way to specify multiple
values for any criterion besides job IDs, e.g. two different project names.
It is recommended that for complex criteria like these, you parse the results
yourself.

=head2 qacct

Get accounting information about an SGE job that has finished.  Returns a
reference to an array of hash refs, where each hash includes (at least)
the following keys:

    job_number
    task_id
    job_name
    cpu
    end_time
    exit_status
    failed
    host_name
    max_vmem
    project
    queue_name
    qsub_time
    slots
    start_time
    wallclock

 Usage: $ra = $sge->qacct(job => $job_number);  # no wildcards
        $ra = $sge->qacct(job => $job_name);    # no wildcards
        $ra = $sge->qacct(job => \@jids);       # multiple jobs
        $ra = $sge->qacct(submitted => 1);      # all submitted jobs
    Args: optional arguments that may be added to above:
          max_wait => $secs
            to retry each second for specified time period; is frequently
            important, since it takes a few seconds for each job to be
            entered into the SGE database.  Defaults to 10 seconds; use
            max_wait => 0 to disable.  Whole, positive integers only.
            Does not wait for all jobs; only that at least one is found.
          failed => 1
            to only return jobs that failed, i.e. have non-zero failed
            flag or exit_status
 Returns: reference to array of hash references;
          undef if no jobs match criteria.

=head1 TODO

Consider making wait process (optionally?) confirm not only that the
precursor processes finished, but that they finished successfully.
Currently, this can be done by using qacct() after wait().

=head1 AUTHORS

 Peter Chines <pchines@mail.nih.gov>

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
