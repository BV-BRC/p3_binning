
#
# Module to encapsulate metagenome binning code.
#
# Name badness. app_params used in two different ways down there. Fix that.
#

package Bio::P3::MetagenomeBinning::MetagenomeBinning;

use P3DataAPI;
use gjoseqlib;
use strict;
use Data::Dumper;
use POSIX;
use Cwd;
use base 'Class::Accessor';
use JSON::XS;
use Module::Metadata;
use IPC::Run;
use File::Basename;
use File::Copy qw(copy);
use File::Path qw(make_path remove_tree);
use Template;
use Text::CSV_XS qw(csv);
use Bio::KBase::AppService::ClientExt;
use Bio::KBase::AppService::AppConfig qw(data_api_url db_host db_user db_pass db_name
					 binning_spades_threads binning_spades_ram
					 bebop_binning_user bebop_binning_key
					 seedtk binning_genome_annotation_clientgroup
					 binning_data_api_url
					);
use DBI;
use File::Slurp;

use Bio::P3::MetagenomeBinning::BebopBinning;

push @INC, seedtk . "/modules/RASTtk/lib";
push @INC, seedtk . "/modules/p3_seedtk/lib";
require BinningReports;
require GEO;

__PACKAGE__->mk_accessors(qw(app app_def params token task_id
			     work_dir assembly_dir stage_dir
			     output_base output_folder 
			     assembly_params spades
			     contigs app_params bebop
			    ));

sub new
{
    my($class) = @_;

    my $self = {
	assembly_params => [],
	app_params => [],
    };

    if (bebop_binning_user && bebop_binning_key)
    {
	my $bebop = Bio::P3::MetagenomeBinning::BebopBinning->new(user => bebop_binning_user,
							       key => bebop_binning_key);
	$self->{bebop} = $bebop;
    }

    return bless $self, $class;
}

#
# Preflight. The CGA app itself has fairly small requirements; it spends most of its
# time waiting on other applications.
#
# We don't mark as a control task, however, because it does have some signficant
# cpu use.
#
sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    my $cpu = 8;

    # if (!$params->{force_local_assembly} && bebop_binning_user && bebop_binning_key)
    # {
    # 	$cpu = 1;
    # }

    #
    # Require the checkv database if we are doing viral binning.
    # Also bump the cpu request up.
    #
    if ($params->{perform_viral_binning})
    {
	$cpu = 8;
	my $db = $ENV{CHECKVDB};
	if (!defined($db))
	{
	    die "Checkv database environment CHECKVDB not defined\n";
	}
	elsif (! -s $db)
	{
	    die "Checkv database not found at $db\n";
	}
    }


    my $pf = {
	cpu => $cpu,
	memory => "128G",
	runtime => 0,
	storage => 0,
	is_control_task => 0,
    };
    return $pf;
}

sub process
{
    my($self, $app, $app_def, $raw_params, $params) = @_;

    $self->app($app);
    $self->app_def($app_def);
    $self->params($params);
    $self->token($app->token);
    $self->task_id($app->task_id);

    print "Process metagenome binning run ", Dumper($app_def, $raw_params, $params);

    my $cwd = getcwd();
    my $assembly_dir = "$cwd/assembly";
    my $work_dir = "$cwd/work";
    my $stage_dir = "$cwd/stage";
    
    -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
    -d $assembly_dir or mkdir $assembly_dir or die "Cannot mkdir $assembly_dir: $!";
    -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";

    $self->work_dir($work_dir);
    $self->assembly_dir($assembly_dir);
    $self->stage_dir($stage_dir);

    my $output_base = $self->params->{output_file};
    my $output_folder = $self->app->result_folder();

    $self->output_base($output_base);
    $self->output_folder($output_folder);

    #
    # Check for exactly one of our input types;
    my @input = grep { $_ } @$params{qw(paired_end_libs single_end_libs contigs srr_ids)};
    if (@input == 0)
    {
	die "No input data specified";
    }
    elsif (@input > 1)
    {
	die "Only one input data item may be specified";
    }

    if (my $val = $params->{paired_end_libs})
    {
	#
	# Check for new form
	#
	if (@$val == 2 && !ref($val->[0]) && !ref($val->[1]))
	{
	    $val = [{ read1 => $val->[0], read2 => $val->[1] }];
	}

	my $size = $self->compute_paired_end_lib_size($val);
	printf STDERR "Paired end library size is %.2f GB\n", $size/1e9;

	if (!$params->{force_local_assembly} &&
	    $self->bebop &&
	    $size > 10_000_000_000 &&
	    @$val == 1 &&
	    ($params->{assembler} eq 'auto' || $params->{assembler} eq 'metaspades')
	   )
	{
	    $self->bebop->assemble_paired_end_libs($output_folder, $val->[0], $self->app->task_id);
	    my $local_contigs = "$assembly_dir/contigs.fasta";
	    $self->contigs($local_contigs);
	    print "Copy from " . $self->output_folder . "/contigs.fasta to $local_contigs\n";
	    $app->workspace->download_file($self->output_folder . "/contigs.fasta",
					   $local_contigs,
					   1, $self->token);
	    if (! -f $local_contigs)
	    {
		die "Local data not found\n";
	    }
	}
	else
	{
	    $self->stage_paired_end_libs($val);
	    $self->assemble();
	}
    }
    elsif (my $val = $params->{single_end_libs})
    {
	$self->stage_single_end_libs($val);
	$self->assemble();
    }
    elsif (my $val = $params->{srr_ids})
    {
	# current bebop code doesn't grok SRA
	
	if (0 && $self->bebop)
	{
	    $self->bebop->assemble_srr_ids($val);
	}
	else
	{
	    $self->stage_srr_ids($val);
	    $self->assemble();
	}
    }
    else
    {
	$self->stage_contigs($params->{contigs});
    }

    my $all_bins = [];
    my @good_results;
    if ($params->{perform_bacterial_annotation} || $params->{perform_bacterial_binning})
    {
	$self->compute_coverage();
	$self->compute_bins();
	my $all_bins = $self->extract_fasta();

	if (@$all_bins == 0)
	{
	    $self->write_empty_bin_report();
	}
	elsif ($params->{perform_bacterial_annotation})
	{
	    @good_results = $self->annotate_bins($all_bins);
	}
    }
    else
    {
	IPC::Run::run(["bins_coverage", $self->contigs(), $self->work_dir]);
	symlink("contigs.fasta", $self->work_dir . "/unbinned.fasta") || die "symlink failed $!";
    }

    if ($params->{perform_viral_annotation} || $params->{perform_viral_binning})
    {
	my $bins = $self->bin_viruses();

	if ($params->{perform_viral_annotation})
	{
	    my @good_viruses = $self->annotate_viruses($bins);
	    $self->write_viral_summary_report(\@good_viruses, $all_bins);
	}
    }
    
    $self->write_summary_report(\@good_results, $all_bins, $self->app->workspace, $self->token);
}

#
# Stage the paired end library data as given in parameters. Ensure we
# have a single pair of contigs (spades metagenome assembler only
# handles one pair of paired-end read sets).
#
# 'paired_end_libs' => [
#		                                        '/olson@patricbrc.org/Binning/Data/SRR2188006_1.fastq.gz',
#		                                        '/olson@patricbrc.org/Binning/Data/SRR2188006_2.fastq.gz'
#		                                      ],
    

sub stage_paired_end_libs
{
    my($self, $libs) = @_;

    my @reads;

    #
    # Check for new form
    #

    my @pairs_1;
    my @pairs_2;

    if (@$libs == 2 && !ref($libs->[0]) && !ref($libs->[1]))
    {
	@reads = @$libs;
    }
    else
    {
	if (@$libs == 0)
	{
	    die "MetagenomeBinning:: stage_paired_end_libs - no libs provided";
	}
	elsif (@$libs > 1)
	{
	    # die "MetagenomeBinning:: stage_paired_end_libs - only one lib may be provided";
	    $self->{megahit_mode} = 1;
	}

	for my $lib (@$libs)
	{
	    @reads = @$lib{qw(read1 read2)};
	    my $staged = $self->app->stage_in(\@reads, $self->stage_dir, 1);
	    push(@pairs_1, $staged->{$lib->{read1}});
	    push(@pairs_2, $staged->{$lib->{read2}});
	}
    }
    my $p1 = join(",", @pairs_1);
    my $p2 = join(",", @pairs_2);
    push(@{$self->assembly_params},
	 ($p1 ? ("-1", $p1) : ()),
	 ($p2 ? ("-2", $p2) : ()),
     );
}

sub compute_paired_end_lib_size
{
    my($self, $libs) = @_;

    my $total = 0;
    
    if (@$libs == 0)
    {
	return $total;
    }
    for my $lib (@$libs)
    {
	my @reads = @$lib{qw(read1 read2)};
	for my $r (@reads)
	{
	    my $stat = eval { $self->app->workspace->stat($r); };
	    $total += $stat if ref($stat);
	}
    }
    return $total;
}

sub stage_single_end_libs
{
    my($self, $libs) = @_;

    my @reads;

    if (@$libs == 0)
    {
	die "MetagenomeBinning:: stage_paired_end_libs - no libs provided";
    }

    $self->{megahit_mode} = 1;

    my @staged;
    for my $lib (@$libs)
    {
	my $reads = $lib->{read};
	my $staged = $self->app->stage_in([$reads], $self->stage_dir, 1);
	push(@staged, $staged->{$lib->{read}});
    }

    my $p1 = join(",", @staged);

    push(@{$self->assembly_params},
	 "-r", $p1);
}

sub stage_srr_ids
{
    my($self, $srr_ids) = @_;

    my $stage = $self->stage_dir;

    my @pairs_1;
    my @pairs_2;
    my @unpaired;
    for my $id (@$srr_ids)
    {
	my $dir = "$stage/$id";
	remove_tree($dir);
	make_path($dir);
	my $ok = IPC::Run::run(["p3-sra", "--id", $id, "--out", $dir]);
	if (!$ok)
	{
	    die "Error staging SRA sample $id\n";
	}
	my @files = glob("$dir/*fastq");
	if (@files == 1)
	{
	    push(@unpaired, $files[0]);
	}
	elsif (@files == 2)
	{
	    push(@pairs_1, $files[0]);
	    push(@pairs_2, $files[1]);
	}
	else
	{
	    die "Unxpected file list from loading $id: @files\n";
	}
    }
    if (@unpaired || @pairs_1 > 1)
    {
	if ($self->params->{assembler} eq 'metaspades')
	{
	    die "This job cannot be processed using metaspades. It includes unpaired reads or more than one paired-read library\n";
	}
	$self->{megahit_mode} = 1;
    }
    my $p1 = join(",", @pairs_1);
    my $p2 = join(",", @pairs_2);
    my $up = join(",", @unpaired);
    push(@{$self->assembly_params},
	 ($p1 ? ("-1", $p1) : ()),
	 ($p2 ? ("-2", $p2) : ()),
	 ($up ? ("-r", $up) : ()),
     );
}

#
# Stage the assembled contigs.
#
sub stage_contigs
{
    my($self, $contigs) = @_;

    my $staged = $self->app->stage_in([$contigs], $self->stage_dir, 1);

    my $file = $staged->{$contigs};
    if (my($unzipped) = $file =~ /(.*)\.gz$/)
    {
	print STDERR "Unzipping $file => $unzipped\n";
	my $rc = system("gunzip", $file);
	if ($rc != 0)
	{
	    die "Error unzipping $file: $rc\n";
	}
	elsif (-s $unzipped)
	{
	    $self->contigs($unzipped);
	}
	else
	{
	    die "Zero-length file $unzipped resulting from unzipping $file\n";
	}
    }
    else
    {
	$self->contigs($staged->{$contigs});
    }
}

#
# Invoke the assembler. We've built a list of assembler parameters during
# the stage-in process. Complete the set of parameters for our
# current configuration and run the assembly.
#
sub assemble
{
    my ($self) = @_;

    if ($self->{megahit_mode} || $self->params->{assembler} eq 'megahit')
    {
	return $self->assemble_with_megahit();
    }
    else
    {
	return $self->assemble_with_spades();
    }
}

sub assemble_with_spades
{
    my($self) = @_;

    my $params = $self->assembly_params;
    push(@$params,
	 "--meta",
	 "--only-assemble",
	 "-o", $self->assembly_dir);

    #
    # Prior code looked at binning_spades_threads and binning_spades_ram
    # to set these parameters; use the scheduler-allocated value instead.
    #
    if (my $cpu = $ENV{P3_ALLOCATED_CPU})
    {
	push(@$params, "--threads", $cpu);
    }


    if (my $mem = $ENV{P3_ALLOCATED_MEMORY})
    {
	my $bytes;
	my %fac = (k => 1024, m => 1024*1024, g => 1024*1024*1024, t => 1024*1024*1024*1024 );
	my($val, $suffix) = $mem =~ /^(.*)([mkgt])$/i;
	if ($suffix)
	{
	    $bytes = $val * $fac{lc($suffix)};
	}
	else
	{
	    $bytes = $mem;
	}
	$mem = int($bytes / (1024*1024*1024));
	push(@$params, "--memory", $mem);
    }
    
    my @cmd = ($self->spades, @$params);
    my $rc = system(@cmd);
    #my $rc = 0;
    if ($rc != 0)
    {
	die "Error running assembly command: @cmd\n";
    }

    $self->app->workspace->save_file_to_file($self->assembly_dir . "/contigs.fasta", {},
					     $self->output_folder . "/contigs.fasta", 'contigs', 1, 1, $self->token);
    $self->app->workspace->save_file_to_file($self->assembly_dir . "/spades.log", {},
					     $self->output_folder . "/spades.log", 'txt', 1, 1, $self->token);
    $self->app->workspace->save_file_to_file($self->assembly_dir . "/params.txt", {},
					     $self->output_folder . "/params.txt", 'txt', 1, 1, $self->token);

    $self->contigs($self->assembly_dir . "/contigs.fasta");
}

sub assemble_with_megahit
{
    my($self) = @_;

    my $params = $self->assembly_params;
    rmdir($self->assembly_dir);
    push(@$params,
	 "-o", $self->assembly_dir);

    #
    # Prior code looked at binning_spades_threads and binning_spades_ram
    # to set these parameters; use the scheduler-allocated value instead.
    #
    my $cpu = $ENV{P3_ALLOCATED_CPU} // 2;
    push(@$params, "-t", $cpu);

    if (my $mem = $ENV{P3_ALLOCATED_MEMORY})
    {
	my $bytes;
	my %fac = (k => 1024, m => 1024*1024, g => 1024*1024*1024, t => 1024*1024*1024*1024 );
	my($val, $suffix) = $mem =~ /^(.*)([mkgt])$/i;
	if ($suffix)
	{
	    $bytes = $val * $fac{lc($suffix)};
	}
	else
	{
	    $bytes = $mem;
	}
	$mem = int($bytes / (1024*1024*1024));
	push(@$params, "--memory", $bytes);
    }
    
    my @cmd = ("megahit", @$params);
    print STDERR "@cmd\n";
    my $rc = system(@cmd);
    #my $rc = 0;
    if ($rc != 0)
    {
	die "Error running assembly command: @cmd\n";
    }

    $self->app->workspace->save_file_to_file($self->assembly_dir . "/final.contigs.fa", {},
					     $self->output_folder . "/contigs.fasta", 'contigs', 1, 1, $self->token);
    $self->app->workspace->save_file_to_file($self->assembly_dir . "/log", {},
					     $self->output_folder . "/megahit.log", 'txt', 1, 1, $self->token);
    $self->app->workspace->save_file_to_file($self->assembly_dir . "/opts.txt", {},
					     $self->output_folder . "/opts.txt", 'txt', 1, 1, $self->token) if -f $self->assembly_dir . "/opts.txt";
    $self->app->workspace->save_file_to_file($self->assembly_dir . "/options.json", {},
					     $self->output_folder . "/options.json", 'json', 1, 1, $self->token) if -f $self->assembly_dir . "/options.json";

    $self->contigs($self->assembly_dir . "/final.contigs.fa");
}

#
# Use bins_coverage to compute coverage. This has a side effect of
# copying the input fasta data to the work directory.
sub compute_coverage
{
    my($self) = @_;

#    local $ENV{PATH} = seedtk . "/bin:$ENV{PATH}";

    my @cmd = ("bins_coverage",
	       "--statistics-file", "coverage.stats.txt",
	       $self->contigs, $self->work_dir);
    my $rc = system(@cmd);

    $rc == 0 or die "Error $rc running coverage: @cmd";
	
    $self->app->workspace->save_file_to_file("coverage.stats.txt", {},
					     $self->output_folder . "/coverage.stats.txt", 'txt', 1, 1, $self->token);
}

sub compute_bins
{
    my($self) = @_;

#    local $ENV{PATH} = seedtk . "/bin:$ENV{PATH}";

    #
    # Work around the change in name of reference data
    my $seedprot1 = "$FIG_Config::p3data/seedprot.fa";
    my $seedprot2 = "$FIG_Config::p3data/seedProt.fa";
    my $seedprot = $seedprot1;
    if (! -f $seedprot)
    {
	$seedprot = $seedprot2;
	if (! -f $seedprot)
	{
	    die "Cannot find either $seedprot1 or $seedprot2\n";
	}
    }

    my @cmd = ("bins_generate",
	       "--dataAPIUrl", binning_data_api_url,
	       "--statistics-file", "bins.stats.txt",
	       "--seedProtFasta", $seedprot,
	       $self->work_dir);
    my $rc = system(@cmd);

    $rc == 0 or die "Error $rc computing bins: @cmd";

    $self->app->workspace->save_file_to_file($self->work_dir . "/unbinned.fasta", {},
					     $self->output_folder . "/unbinned.fasta", 'contigs', 1, 1, $self->token);
    $self->app->workspace->save_file_to_file($self->work_dir . "/unplaced.fasta", {},
					     $self->output_folder . "/unplaced.fasta", 'contigs', 1, 1, $self->token);
    $self->app->workspace->save_file_to_file("bins.stats.txt", {},
					     $self->output_folder . "/bins.stats.txt", 'txt', 1, 1, $self->token);
    $self->app->workspace->save_file_to_file("bins.stats.txt", {},
					     $self->output_folder . "/bins.stats.txt", 'txt', 1, 1, $self->token);
}

sub extract_fasta
{
    my($self) = @_;

    #my @cmd = ("bins_fasta", $self->work_dir);
    #my $rc = system(@cmd);
    # $rc == 0 or die "Error $rc computing bins: @cmd";

    #
    # We essentially inline the bins_fasta code here since we want
    # to extract the metadata from the bins.json file as we go.
    #

    local $/ = "//\n";
    open(BINS, "<", $self->work_dir . "/bins.json") or die "Cannot open " . $self->work_dir . "/bins.json: $!" ;
    open(SAMPLE, "<", $self->work_dir . "/sample.fasta") or die "Cannot open " . $self->work_dir . "/sample.fasta: $!" ;

    #
    # App params exemplar for annotation submission
    #
    # {
    # "contigs": "/olson@patricbrc.org/home/buchnera.fa",
    #     "scientific_name": "Buchnera aphidicola",
    #     "code": 11,
    #     "domain": "Bacteria",
    #     "taxonomy_id": 107806,
    #     "output_path": "/olson@patricbrc.org/home/output",
    #     "output_file": "buch32"
    #     }

    my $app_list = $self->app_params;

    my $api = P3DataAPI->new(binning_data_api_url);

    my $idx = 1;

    my $all_bins = [];
    while (defined(my $bin_txt = <BINS>))
    {
	chomp $bin_txt;
	my $bin;
	eval { $bin = decode_json($bin_txt); };
	if ($@)
	{
	    warn "Bad parse on '$bin_txt'\n";
	    last;
	}
	push(@$all_bins, $bin);

	my $taxon_id = $bin->{taxonID};
	print "$bin->{name} $taxon_id\n";
	my %want;
	for my $c (@{$bin->{contigs}})
	{
	    $want{$c->[0]} = 1;
	}

	my $bin_base_name = "bin.$idx.$taxon_id";
	my $bin_name = "$bin_base_name.fa";
	$bin->{binFastaFile} = $bin_name;
	$bin->{binIndex} = $idx;
	$idx++;
	my $bin_fa = $self->work_dir . "/$bin_name";
	open(BIN, ">", $bin_fa) or die "Cannot write $bin_fa: $!";
	seek(SAMPLE, 0, 0);

	local $/ = "\n";
	while (my($id, $def, $seq) = read_next_fasta_seq(\*SAMPLE))
	{
	    if ($want{$id})
	    {
		write_fasta(\*BIN, [[$id, $def, $seq]]);
	    }
	}
	close(BIN);

	my $ws_path = $self->output_folder . "/$bin_name";
	$self->app->workspace->save_file_to_file($bin_fa, $bin, $ws_path, 'contigs', 1, 1, $self->token);
	$bin->{binFastaPath} = $ws_path;

	my $code = 11;
	my $domain = 'Bacteria';
	my @res = $api->query("taxonomy", ["eq", 'taxon_id', $taxon_id], ["select", "genetic_code,lineage_names"]);
	if (@res)
	{
	    my $lineage;
	    ($code, $lineage) = @{$res[0]}{'genetic_code', 'lineage_names'};
	    shift @$lineage if ($lineage->[0] =~ /cellular organisms/);
	    $domain = $lineage->[0];
	}

	$bin->{domain} = $domain;
	$bin->{geneticCode} = $code;

	my $descr = {
	    contigs => $ws_path,
	    code => $code,
	    domain => $domain,
	    scientific_name=> $bin->{name},
	    taxonomy_id => $taxon_id,
	    reference_genome_id => $bin->{refGenomes}->[0],
	    output_path => $self->output_folder,
	    output_file => $bin_base_name,
#	    _parent_job => $self->app->task_id,
	    queue_nowait => 1,
	    analyze_quality => 1,
	    ($self->params->{skip_indexing} ? (skip_indexing => 1) : ()),
	    recipe => $self->params->{recipe},
	    (binning_genome_annotation_clientgroup ? (_clientgroup => binning_genome_annotation_clientgroup) : ()),
	};
	push(@$app_list, $descr);
    }
    my $json = JSON::XS->new->pretty(1)->canonical(1);
    $self->app->workspace->save_data_to_file($json->encode($all_bins), {},
					     $self->output_folder . '/bins.json', 'json', 1, 1, $self->token);
    print "SAVE to " . $self->output_folder . "/bins.json\n";

    close(SAMPLE);
    close(BINS);
    print STDERR Dumper($self->app_params);

    #
    # Return the bins so that we can cleanly terminate the job if no bins
    # were found. Also used later for reporting.
    #
    return $all_bins;
}

sub write_empty_bin_report
{
    my($self) = @_;

    #
    # No bins found. Write a simple HTML stating that.
    #
    my $report = "<h1>No bins found</h1>\n<p>No bins were found in this sample.\n";
    $self->app->workspace->save_data_to_file($report, {},
					     $self->output_folder . "/BinningReport.html", 'html', 1, 0, $self->app->token);
}

sub annotate_bins
{
    my($self, $all_bins) = @_;

    my @good_results;

    my $annotations_inline = $ENV{P3_BINNING_ANNOTATIONS_INLINE} || $self->params->{force_inline_annotation};
    if ($annotations_inline)
    {
	@good_results = $self->compute_annotations_local();
    }
    else
    {
	@good_results = $self->compute_annotations_cluster();
    }


}

#
# Use vbins_generate to bin the viruses.
# vbins_generate /disks/patric-common/runtime/checkv-db/checkv-db-v1.0 `pwd`
#
sub bin_viruses
{
    my($self) = @_;

    if (! -s $self->work_dir . "/unbinned.fasta")
    {
	warn "No unbinned data for viral binning\n";
	return ();
    }

    my $cmd = "vbins_generate";
    my @params;
    if (my $cpu = $ENV{P3_ALLOCATED_CPU})
    {
	push(@params, "--threads", $cpu);
    }
    push(@params, $ENV{CHECKVDB}, $self->work_dir);

    my @cmd = ($cmd, @params);
    print STDERR "@cmd\n";
    my $rc = system(@cmd);
    if ($rc != 0)
    {
	warn "Viral binning failed with $rc: @cmd\n";
	return ();
    }

    my $vbins = eval { csv(in => $self->work_dir . "/vbins.tsv", headers => "auto", sep_char => "\t") };
    if (@$vbins == 0)
    {
	warn "No viral bins created\n";
	$self->write_empty_vbin_report();
	return ();
    }

    my @ret_vbins;
    for my $vbin (@$vbins)
    {
	#
	# Force conversion to numeric for the numeric fields.
	#
	for my $f (qw(pct_error length completeness coverage))
	{
	    $vbin->{$f} += 0 if exists $vbin->{$f};
	}
	
	my $bin_name = "vBin$vbin->{bin}.fa";
	my $fa_file = $self->work_dir . "/$bin_name";
	if (! -f $fa_file)
	{
	    print STDERR "Could not find $fa_file\n";
	    next;
	}
	my $ws_path = $self->output_folder . "/$bin_name";
	$self->app->workspace->save_file_to_file($fa_file, $vbin, $ws_path, 'contigs', 1, 1, $self->token);
	$vbin->{ws_path} = $ws_path;
	$vbin->{bin_name} = $bin_name;
	push @ret_vbins, $vbin;
    }
    # Don't need this with the other report that looks better.
    # $self->app->workspace->save_file_to_file($self->work_dir . "/vbins.html",
    #					     {}, $self->output_folder . "/ViralBins.html", 'html', 1, 1, $self->token);
    return \@ret_vbins;
}

sub annotate_viruses
{
    my($self, $bins) = @_;

    my @good_results;

    my $api = P3DataAPI->new(data_api_url);
    my $app_spec = $self->find_app_spec("GenomeAnnotation");

    my $sub_time = time;
    my $n = 1;
    my $recipe = $self->params->{viral_recipe} // "viral";
    for my $bin (@$bins)
    {
    	my $code = 1;
	my $domain = 'Viruses';
	my $taxon_id = $bin->{taxon_id};
	my @res = $api->query("taxonomy", ["eq", 'taxon_id', $taxon_id], ["select", "genetic_code,lineage_names"]);
	if (@res)
	{
	    my $lineage;
	    ($code, $lineage) = @{$res[0]}{'genetic_code', 'lineage_names'};
	    shift @$lineage if ($lineage->[0] =~ /cellular organisms/);
	    $domain = $lineage->[0];
	}

	$bin->{domain} = $domain;
	$bin->{genetic_code} = $code;

	(my $bin_base_name = $bin->{bin_name}) =~ s/\.[^.]+$//;

	my $descr = {
	    contigs => $bin->{ws_path},
	    code => $code,
	    domain => $domain,
	    scientific_name=> $bin->{name},
	    taxonomy_id => $taxon_id,
	    output_path => $self->output_folder,
	    output_file => $bin_base_name,
#	    _parent_job => $self->app->task_id,
	    queue_nowait => 1,
	    analyze_quality => 0,
	    ($self->params->{skip_indexing} ? (skip_indexing => 1) : ()),
	    recipe => $recipe,
	};

	my $json = JSON::XS->new->pretty->canonical();

	my $tmp = File::Temp->new;
	print $tmp $json->encode($descr);
	close($tmp);
	my @cmd = ("App-GenomeAnnotation", "xx", $app_spec, "$tmp");
	print STDERR "Run annotation: @cmd\n";
	my $start = time;
	my $rc = system(@cmd);
	my $end = time;
	if ($rc != 0)
	{
	    warn "Annotation failed with rc=$rc\n";
	    next;
	}

	my $gto;
	my $genome_id;
	eval {
	    $gto = $self->app->workspace->download_json("$descr->{output_path}/.$descr->{output_file}/$descr->{output_file}.genome");
	    $genome_id = $gto->{id};
	};
	
	push(@good_results, {
	    id => $n++,
	    genome_id => $genome_id,
	    app => "App-GenomeAnnotation",
	    parameters => $descr,
	    user_id => 'immediate-user',
	    submit_time => strftime("%Y-%m-%dT%H:%M:%SZ", gmtime $sub_time),
	    start_time => strftime("%Y-%m-%dT%H:%M:%SZ", gmtime $start),
	    completed_time => strftime("%Y-%m-%dT%H:%M:%SZ", gmtime $end),
	    vbin => $bin,
	});
    }
    return @good_results;
}

sub write_empty_vbin_report
{
    my($self) = @_;

    #
    # No bins found. Write a simple HTML stating that.
    #
    my $report = "<h1>No viral bins found</h1>\n<p>No viral bins were found in this sample.\n";
    $self->app->workspace->save_data_to_file($report, {},
					     $self->output_folder . "/ViralBinningReport.html", 'html', 1, 0, $self->app->token);
}

sub write_viral_summary_report
{
    my($self, $good, $all_bins) = @_;

    eval {

	#
	# Find template.
	#
	
	my $mpath = Module::Metadata->find_module_by_name(__PACKAGE__);
	$mpath =~ s/\.pm$//;
	
	my $summary_tt = dirname($mpath) . "/ViralBinSummary.tt";
	-f $summary_tt or die "Summary not found at $summary_tt\n";

	my $url_base = 'https://alpha.bv-brc.org/view/Genome';

	my $bins = [];
	my %vars = (bins => $bins,
		    params => $self->params,
		    job_id => $self->task_id,
		    );
	for my $bin (@$good)
	{
	    my $genome_url = "$url_base/$bin->{genome_id}";
	    my $val = {
		vbin => $bin->{vbin},
		genome_id => $bin->{genome_id},
		genome_url => $genome_url,
	    };
	    push(@$bins, $val);
	}

	my $templ = Template->new(ABSOLUTE => 1);
	my $html;

	print STDERR Dumper(\%vars);
	$templ->process($summary_tt, \%vars, \$html);

	my $output_path = $self->params->{output_path} . "/." . $self->params->{output_file};
	$self->app->workspace->save_data_to_file($html, {},
						 "$output_path/ViralBinningReport.html", 'html', 1, 0);
    };
    if ($@)
    {
	warn "Error creating final viral binning report: $@";
    }
}



sub write_db_record
{
    my($self, $n_children) = @_;
    
    my $dsn = "DBI:mysql:database=" . db_name . ";host=" . db_host;
    my $dbh = DBI->connect($dsn, db_user, db_pass, { RaiseError => 1, AutoCommit => 0 });

    my $json = JSON::XS->new->pretty(1)->canonical(1);

    $dbh->do(qq(INSERT INTO JobGroup (parent_job, children_created, parent_app, app_spec, app_params)
		VALUES (?, ?, ?, ?, ?)), undef,
	     $self->app->task_id, $n_children, "MetagenomeBinning",
	     $json->encode($self->app_def), $json->encode($self->params));
    $dbh->commit();
}

sub compute_annotations_cluster
{
    my($self) = @_;

    my $app_service = Bio::KBase::AppService::ClientExt->new();

    my @tasks = $self->submit_annotations($app_service);

    print STDERR "Awaiting completion of " . scalar(@tasks) . " annotations\n";
    my $results = $app_service->await_task_completion(\@tasks, 10, 0);
    print STDERR "Tasks completed\n";
	
    #
    # Examine task output to ensure all succeeded
    #
    my @good_results;
    my $fail = 0;
    for my $res (@$results)
    {
	if ($res->{status} eq 'completed')
	{
	    push(@good_results, $res);
	}
	else
	{
	    warn "Task $res->{id} resulted with unsuccessful status $res->{status}\n" . Dumper($res);
	    $fail++;
	}
    }
    
    if ($fail > 0)
    {
	if ($fail == @$results)
	{
	    die "Annotation failed on all $fail bins\n";
	}
	else
	{
	    my $n = @$results;
	    warn "Annotation failed on $fail of $n bins, continuing\n";
	}
    }
    return @good_results;
}	

#
# Compute annotations inline by invoking the annotation script.
# Mostly used for testing, but may be useful for standalone implementation.
#
# Returns a list of Task hashes.
#
sub compute_annotations_local
{
    my($self) = @_;

    #
    # Need to find our app specs. Different locations if we are in production or development.
    #

    my $app_spec = $self->find_app_spec("GenomeAnnotation");

    my @good_results;

    my $sub_time = time;
    my $n = 1;
    for my $task (@{$self->app_params})
    {
	my $tmp = File::Temp->new;
	print $tmp encode_json($task);
	close($tmp);
	my @cmd = ("App-GenomeAnnotation", "xx", $app_spec, "$tmp");
	# my @cmd = ("bash", "-c", "source /vol/patric3/production/P3Slurm2/tst_deployment/user-env.sh; App-GenomeAnnotation xx $app_spec $tmp");
	print STDERR "Run annotation: @cmd\n";
	my $start = time;
	my $rc = system(@cmd);
	my $end = time;
	if ($rc != 0)
	{
	    warn "Annotation failed with rc=$rc\n";
	    next;
	}
	push(@good_results, {
	    id => $n++,
	    app => "App-GenomeAnnotation",
	    parameters => $task,
	    user_id => 'immediate-user',
	    submit_time => strftime("%Y-%m-%dT%H:%M:%SZ", gmtime $sub_time),
	    start_time => strftime("%Y-%m-%dT%H:%M:%SZ", gmtime $start),
	    completed_time => strftime("%Y-%m-%dT%H:%M:%SZ", gmtime $end),
	});
    }
print Dumper(GOOD => \@good_results);
    return @good_results;
}
    

sub submit_annotations
{
    my($self, $client) = @_;

    if ($ENV{BINNING_PRINT_SUBS})
    {
	my $json = JSON::XS->new->pretty(1)->canonical(1);
	my $i = 1;
	for my $task (@{$self->app_params})
	{
	    open(OUT, ">", "binning.job.$i") or die;
	    print OUT $json->encode($task);
	    close(OUT);
	    $i++;
	}
	exit;
    }

    #
    # No longer needed with new code that waits for annos.
    # $self->write_db_record(scalar @{$self->app_params});

    my $start_params = {
	parent_id => $self->app->task_id,
	workspace => $self->output_folder,
    };
    my @tasks;
    for my $task (@{$self->app_params})
    {
	my $submitted = $client->start_app2("GenomeAnnotation", $task, $start_params);
	push(@tasks, $submitted);
    }
    return @tasks;
}
    

#
# Write the summary. $tasks is the list of annotation tasks from the app service; included therein
# are the parameters to the annotation runs which includes the output locations. Use those to
# pull the genome objects.
#
sub write_summary_report
{
    my($self, $tasks, $bins_report, $ws, $token) = @_;

    my @genomes;
    my @geos;
    my %report_url_map;

    my %geo_opts = (
	detail => 2,
	p3 => P3DataAPI->new(data_api_url, $token->token),
    );
    
    local $FIG_Config::global = seedtk . "/data";
    for my $task (@$tasks)
    {
	my $params = $task->{parameters};
	my $name = $params->{output_file};
	my $genome_path = $params->{output_path} . "/.$name";
	my $gto_path = "$genome_path/$name.genome";

	#
	# we need to convert to GEOs for the binning reports.
	#
	my $temp = File::Temp->new(UNLINK => 1);
	my $qual_temp = File::Temp->new(UNLINK => 1);

	print "$genome_path/genome_quality_details.txt\n";
	eval {
	    $ws->copy_files_to_handles(1, $token,
				       [[$gto_path, $temp],
					]);
	};
	warn "Error copying $gto_path: $@" if $@;
	eval {
	    $ws->copy_files_to_handles(1, $token,
				       [["$genome_path/genome_quality_details.txt", $qual_temp],
					]);
	};
	warn "Error copying $genome_path/genome_quality_details.txt: $@" if $@;
	close($temp);
	close($qual_temp);

	if (! -s "$temp")
	{
	    warn "Could not load $gto_path\n";
	    next;
	}
		
	my $gret = GEO->CreateFromGtoFiles(["$temp"], %geo_opts);
	my($geo) = values %$gret;

	if (-s "$qual_temp")
	{
	    $geo->AddQuality("$qual_temp");
	    write_file("$name.geo", Dumper($geo));
	    push(@geos, $geo);
	}
	else
	{
	    warn "Could not read qual file $genome_path/genome_quality_details.txt\n";
	}
	my $genome_id = $geo->id;
	print "$genome_id: $geo->{name}\n";
	push(@genomes, $genome_id);

	#
	# We assume the report URL is available in the same workspace
	# directory as the genome.
	#
	# The genome is actually in the same subtree as the summary report,
	# so we use a relative path instead of $genome_path.
	#
	my $report_url = ".$name/GenomeReport.html";
	#my $report_url = "https://www.patricbrc.org/workspace$report_path";
	$report_url_map{$genome_id} = $report_url;
    }

    #
    # Write the genome group
    #

    my $params = $self->params;

    my $group_path;
    if (my $group = $params->{genome_group})
    {
	my $home;
	if ($token->token =~ /(^|\|)un=([^|]+)/)
	{
	    my $un = $2;
	    $home = "/$un/home";
	}

	if ($home)
	{
	    $group_path = "$home/Genome Groups/$group";
	    
	    my $group_data = { id_list => { genome_id => \@genomes } };
	    my $group_txt = encode_json($group_data);
	    
	    my $res = $ws->create({
		objects => [[$group_path, "genome_group", {}, $group_txt]],
		permission => "w",
		overwrite => 1,
	    });
	    print STDERR Dumper(group_create => $res);
	}
	else
	{
	    warn "Cannot find home path '$home'\n";
	}
    }
    #
    # Generate the binning report. We need to load the various reports into memory to do this.
    #
    eval {

	#
	# Find template.
	#
	
	my $mpath = Module::Metadata->find_module_by_name("BinningReports");
	$mpath =~ s/\.pm$//;
	
	my $summary_tt = "$mpath/summary.tt";
	-f $summary_tt or die "Summary not found at $summary_tt\n";

	#
	# Read SEEDtk role map.
	#
	my %role_map;
	if (open(R, "<", seedtk . "/data/roles.in.subsystems"))
	{
	    while (<R>)
	    {
		chomp;
		my($abbr, $hash, $role) = split(/\t/);
		$role_map{$abbr} = $role;
	    }
	    close(R);
	}

	my $html = BinningReports::Summary($self->task_id, $params, $bins_report, $summary_tt,
					   $group_path, \@geos, \%report_url_map, "/view/Genome");

	my $output_path = $params->{output_path} . "/." . $params->{output_file};
	$ws->save_data_to_file($html, {},
			       "$output_path/BinningReport.html", 'html', 1, 0, $token);
    };
    if ($@)
    {
	warn "Error creating final report: $@";
    }
}

sub find_app_spec
{
    my($self, $app) = @_;
    
    my $top = $ENV{KB_TOP};
    my $specs_deploy = "$top/services/app_service/app_specs";

    my $app_spec;
    
    if (-d $specs_deploy)
    {
	$app_spec = "$specs_deploy/$app.json";
    }
    else
    {
	my @specs_dev = glob("$top/modules/*/app_specs");
	for my $path (@specs_dev)
	{
	    my $s = "$path/$app.json";
	    if (-s $s)
	    {
		$app_spec = $s;
		last;
	    }
	}
	-s $app_spec or 
	    die "cannot find specs file for $app in $specs_deploy or @specs_dev\n";
    }
    -s $app_spec or die "Spec file $app_spec does not exist\n";
    return $app_spec;
}


1;
