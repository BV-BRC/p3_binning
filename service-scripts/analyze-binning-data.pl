use Data::Dumper;
use strict;
use JSON::XS;
use Bio::KBase::AppService::SchedulerDB;
use Bio::P3::Workspace::WorkspaceClientExt;

my $ws = Bio::P3::Workspace::WorkspaceClientExt->new;

my $db = Bio::KBase::AppService::SchedulerDB->new;

my $res = $db->dbh->selectall_arrayref(qq(SELECT id, owner, start_time, finish_time, timediff(finish_time, start_time), output_path, output_file, params, maxrss
					  FROM TaskWithActiveJob
					  WHERE application_id = 'MetagenomeBinning' AND
					     cluster_id = 'P3Slurm' AND
					     state_code = 'C'));
for my $e (@$res)
{
    my($id, $owner, $start, $finish, $elap, $path, $file, $params, $mem) = @$e;

    die Dumper($e);
}

__END__

    my $path1 = "$path/.$file/details/${file}_run_details.json";
    my $path2 = "$path/.$file/${file}_run_details.json";
    my $path3 = "$path/.$file/${file}_run_details.txt";
    my $txt;
    open(my $fh, ">", \$txt) or die $!;
    eval { $ws->copy_files_to_handles(1, undef, [[$path1, $fh]], { admin => 1 }); };
    if ($@)
    {
	eval { $ws->copy_files_to_handles(1, undef, [[$path2, $fh]], { admin => 1 }); };
    }
    if ($@)
    {
	eval { $ws->copy_files_to_handles(1, undef, [[$path3, $fh]], { admin => 1 }); };
    }
    
    warn "$path1: $@ " if $@;
    next unless $txt;

    my $dat = eval { decode_json($txt); };
    next unless $dat;

    next if @{$dat->{problem}};

    my $tot_sz = 0;
    while (my($readkey, $rdata) = each (%{$dat->{reads}}))
    {
	my $sz = $rdata->{num_reads} * $rdata->{avg_len};
	$tot_sz += $sz;
    }

    print join("\t", $id, $owner, $tot_sz, $elap, $mem), "\n";
}


