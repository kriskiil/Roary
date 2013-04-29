package Bio::PanGenome::CommandLine::QueryPanGenome;

# ABSTRACT: Take in a groups file and the protein fasta files and output selected data

=head1 SYNOPSIS

Take in a groups file and the protein fasta files and output selected data

=cut

use Moose;
use Getopt::Long qw(GetOptionsFromArray);
use Bio::PanGenome::AnalyseGroups;
use Bio::PanGenome::Output::GroupsMultifastas;
use Bio::PanGenome::Output::OneGenePerGroupFasta;
use Bio::PanGenome::Output::QueryGroups;

has 'args'        => ( is => 'ro', isa => 'ArrayRef', required => 1 );
has 'script_name' => ( is => 'ro', isa => 'Str',      required => 1 );
has 'help'        => ( is => 'rw', isa => 'Bool',     default  => 0 );

has 'fasta_files'     => ( is => 'rw', isa => 'ArrayRef' );
has 'groups_filename' => ( is => 'rw', isa => 'Str' );
has 'group_names'     => ( is => 'rw', isa => 'ArrayRef' );
has 'output_filename' => ( is => 'rw', isa => 'Str', default => 'pan_genome_results' );
has 'action'          => ( is => 'rw', isa => 'Str', default => 'one_gene_per_group' );

has '_error_message' => ( is => 'rw', isa => 'Str' );

sub BUILD {
    my ($self) = @_;

    my ( $fasta_files, $output_filename, $groups_filename, @group_names, $action, $help );

    GetOptionsFromArray(
        $self->args,
        'o|output=s'          => \$output_filename,
        'g|groups_filename=s' => \$groups_filename,
        'n|group_names=s'     => \@group_names,
        'a|action=s'          => \$action,
        'h|help'              => \$help,
    );

    if ( @{ $self->args } == 0 ) {
        $self->_error_message("Error: You need to provide a FASTA file");
    }

    $self->output_filename($output_filename) if ( defined($output_filename) );
    $self->action($action)                   if ( defined($action) );
    if ( defined($groups_filename) && ( -e $groups_filename ) ) {
        $self->groups_filename($groups_filename);
    }
    else {
        $self->_error_message("Error: Cant access the groups file $groups_filename");
    }

    @group_names = split( /,/, join( ',', @group_names ) );
    $self->group_names( \@group_names ) if (@group_names);

    for my $filename ( @{ $self->args } ) {
        if ( !-e $filename ) {
            $self->_error_message("Error: Cant access file $filename");
            last;
        }
    }
    $self->fasta_files( $self->args );

}

sub run {
    my ($self) = @_;

    ( !$self->help ) or die $self->usage_text;
    if ( defined( $self->_error_message ) ) {
        print $self->_error_message . "\n";
        die $self->usage_text;
    }

    my $analyse_groups_obj = Bio::PanGenome::AnalyseGroups->new(
        fasta_files     => $self->fasta_files,
        groups_filename => $self->groups_filename,
    );

    if ( $self->action eq 'one_gene_per_group' ) {
        my $one_gene_per_fasta = Bio::PanGenome::Output::OneGenePerGroupFasta->new(
            analyse_groups  => $analyse_groups_obj,
            output_filename => $self->output_filename
        );
        $one_gene_per_fasta->create_file();
    }
    elsif ( $self->action eq 'union' ) {
        my $query_groups = Bio::PanGenome::Output::QueryGroups->new(
            analyse_groups        => $analyse_groups_obj,
            output_union_filename => $self->output_filename,
            input_filenames       => $self->fasta_files
        );
        $query_groups->groups_union();
    }
    elsif ( $self->action eq 'intersection' ) {
        my $query_groups = Bio::PanGenome::Output::QueryGroups->new(
            analyse_groups               => $analyse_groups_obj,
            output_intersection_filename => $self->output_filename, 
            input_filenames => $self->fasta_files
        );
        $query_groups->groups_intersection();
    }
    elsif ( $self->action eq 'complement' ) {
        my $query_groups = Bio::PanGenome::Output::QueryGroups->new(
            analyse_groups             => $analyse_groups_obj,
            output_complement_filename => $self->output_filename, 
            input_filenames => $self->fasta_files
        );
        $query_groups->groups_complement();
    }
    elsif ( $self->action eq 'gene_multifasta' && defined( $self->group_names ) ) {
        my $group_multi_fastas = Bio::PanGenome::Output::GroupsMultifastas->new(
            group_names          => $self->group_names,
            analyse_groups       => $analyse_groups_obj,
            output_filename_base => $self->output_filename
        );
        $group_multi_fastas->create_files();
    }
    else {
        print "Nothing done\n";
    }
}

sub usage_text {
    my ($self) = @_;

    return <<USAGE;
    Usage: query_pan_genome [options]
    Take in a groups file and the protein fasta files and output selected data
    
    # Create a FASTA file with one gene per group (representative pan genome)
    query_pan_genome -a one_gene_per_group -g groupfile example.faa
    
    # Provide an output filename
    query_pan_genome  -a one_gene_per_group -g groupfile -o results.fa *.faa
    
    # Create multifasta files for each group/gene passed in
    query_pan_genome  -a gene_multifasta -g groupfile -n gryA,mecA,abc *.faa
    
    # Union
    query_pan_genome  -a union -g groupfile *.faa
    
    # Intersection
    query_pan_genome  -a intersection -g groupfile *.faa

    # Complement (Union minus Intersection)
    query_pan_genome  -a complement -g groupfile *.faa

    # This help message
    query_pan_genome -h

USAGE
}

__PACKAGE__->meta->make_immutable;
no Moose;
1;
