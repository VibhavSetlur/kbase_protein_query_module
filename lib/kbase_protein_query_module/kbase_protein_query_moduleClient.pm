package kbase_protein_query_module::kbase_protein_query_moduleClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

kbase_protein_query_module::kbase_protein_query_moduleClient

=head1 DESCRIPTION


A KBase module: kbase_protein_query_module

This module provides comprehensive protein query analysis capabilities using UniProt IDs as the canonical identifier:

COMPREHENSIVE ANALYSIS WORKFLOW:
1. CheckProteinExistence: Verify protein exists using UniProt ID, optionally generate embedding
2. GenerateProteinEmbeddings: Create embeddings from sequence input or protein check results
3. AssignProteinFamily: Assign proteins to families using similarity to centroids
4. FindTopMatches: Perform similarity search within families
5. SummarizeAndVisualize: Generate comprehensive HTML reports with network analysis

ADVANCED CAPABILITIES:
- UniProt ID canonical identifier system (exact match only)
- ESM-2 protein language model for embedding generation
- Efficient FAISS-based similarity search and clustering
- Family assignment using binary centroid similarity
- Comprehensive metadata storage and retrieval
- HTML report generation with network visualization
- Workspace object management for downstream analysis
- Bioinformatics integration with protein databases
- Network analysis and protein relationship mapping
- Advanced similarity metrics and statistical analysis

Authors: Vibhav Setlur
Contact: https://kbase.us/contact-us/


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => kbase_protein_query_module::kbase_protein_query_moduleClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my %arg_hash2 = @args;
	if (exists $arg_hash2{"token"}) {
	    $self->{token} = $arg_hash2{"token"};
	} elsif (exists $arg_hash2{"user_id"}) {
	    my $token = Bio::KBase::AuthToken->new(@args);
	    if (!$token->error_message) {
	        $self->{token} = $token->token;
	    }
	}
	
	if (exists $self->{token})
	{
	    $self->{client}->{token} = $self->{token};
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 check_protein_existence

  $output = $obj->check_protein_existence($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.CheckProteinExistenceResults
CheckProteinExistenceResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	exists has a value which is an int
	family_id has a value which is a string
	metadata has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	summary has a value which is a string
	protein_existence_result_ref has a value which is a string
	embedding_result_ref has a value which is a string

</pre>

=end html

=begin text

$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.CheckProteinExistenceResults
CheckProteinExistenceResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	exists has a value which is an int
	family_id has a value which is a string
	metadata has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	summary has a value which is a string
	protein_existence_result_ref has a value which is a string
	embedding_result_ref has a value which is a string


=end text

=item Description



=back

=cut

 sub check_protein_existence
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function check_protein_existence (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to check_protein_existence:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'check_protein_existence');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kbase_protein_query_module.check_protein_existence",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'check_protein_existence',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method check_protein_existence",
					    status_line => $self->{client}->status_line,
					    method_name => 'check_protein_existence',
				       );
    }
}
 


=head2 generate_protein_embedding

  $output = $obj->generate_protein_embedding($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.GenerateProteinEmbeddingResults
GenerateProteinEmbeddingResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	embedding_result_ref has a value which is a string
	summary has a value which is a string
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	embedding_norm has a value which is a float
	sequence_length has a value which is an int
	embedding_dim has a value which is an int

</pre>

=end html

=begin text

$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.GenerateProteinEmbeddingResults
GenerateProteinEmbeddingResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	embedding_result_ref has a value which is a string
	summary has a value which is a string
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	embedding_norm has a value which is a float
	sequence_length has a value which is an int
	embedding_dim has a value which is an int


=end text

=item Description



=back

=cut

 sub generate_protein_embedding
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function generate_protein_embedding (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to generate_protein_embedding:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'generate_protein_embedding');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kbase_protein_query_module.generate_protein_embedding",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'generate_protein_embedding',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method generate_protein_embedding",
					    status_line => $self->{client}->status_line,
					    method_name => 'generate_protein_embedding',
				       );
    }
}
 


=head2 assign_family_fast

  $output = $obj->assign_family_fast($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.AssignFamilyFastResults
AssignFamilyFastResults is a reference to a hash where the following keys are defined:
	family_id has a value which is a string
	confidence has a value which is a float
	eigenprotein_id has a value which is a string
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	family_assignment_result_ref has a value which is a string

</pre>

=end html

=begin text

$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.AssignFamilyFastResults
AssignFamilyFastResults is a reference to a hash where the following keys are defined:
	family_id has a value which is a string
	confidence has a value which is a float
	eigenprotein_id has a value which is a string
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	family_assignment_result_ref has a value which is a string


=end text

=item Description



=back

=cut

 sub assign_family_fast
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function assign_family_fast (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to assign_family_fast:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'assign_family_fast');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kbase_protein_query_module.assign_family_fast",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'assign_family_fast',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method assign_family_fast",
					    status_line => $self->{client}->status_line,
					    method_name => 'assign_family_fast',
				       );
    }
}
 


=head2 find_top_matches_from_embedding

  $output = $obj->find_top_matches_from_embedding($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.FindTopMatchesFromEmbeddingResults
FindTopMatchesFromEmbeddingResults is a reference to a hash where the following keys are defined:
	matches has a value which is a reference to a list where each element is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	summary has a value which is a string
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	family_id has a value which is a string
	top_n has a value which is an int
	similarity_stats has a value which is a reference to a hash where the key is a string and the value is a float
	similarity_search_result_ref has a value which is a string

</pre>

=end html

=begin text

$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.FindTopMatchesFromEmbeddingResults
FindTopMatchesFromEmbeddingResults is a reference to a hash where the following keys are defined:
	matches has a value which is a reference to a list where each element is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	summary has a value which is a string
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	family_id has a value which is a string
	top_n has a value which is an int
	similarity_stats has a value which is a reference to a hash where the key is a string and the value is a float
	similarity_search_result_ref has a value which is a string


=end text

=item Description



=back

=cut

 sub find_top_matches_from_embedding
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function find_top_matches_from_embedding (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to find_top_matches_from_embedding:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'find_top_matches_from_embedding');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kbase_protein_query_module.find_top_matches_from_embedding",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'find_top_matches_from_embedding',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method find_top_matches_from_embedding",
					    status_line => $self->{client}->status_line,
					    method_name => 'find_top_matches_from_embedding',
				       );
    }
}
 


=head2 summarize_and_visualize_results

  $output = $obj->summarize_and_visualize_results($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.SummarizeAndVisualizeResultsResults
SummarizeAndVisualizeResultsResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	output_dir has a value which is a string
	summary has a value which is a string
	html_report_path has a value which is a string
	sequence_analysis_ref has a value which is a string

</pre>

=end html

=begin text

$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.SummarizeAndVisualizeResultsResults
SummarizeAndVisualizeResultsResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	output_dir has a value which is a string
	summary has a value which is a string
	html_report_path has a value which is a string
	sequence_analysis_ref has a value which is a string


=end text

=item Description



=back

=cut

 sub summarize_and_visualize_results
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function summarize_and_visualize_results (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to summarize_and_visualize_results:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'summarize_and_visualize_results');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kbase_protein_query_module.summarize_and_visualize_results",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'summarize_and_visualize_results',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method summarize_and_visualize_results",
					    status_line => $self->{client}->status_line,
					    method_name => 'summarize_and_visualize_results',
				       );
    }
}
 


=head2 run_protein_query_analysis

  $output = $obj->run_protein_query_analysis($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.ProteinQueryAnalysisResults
ProteinQueryAnalysisResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	analysis_result_ref has a value which is a string
	summary has a value which is a string
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	html_report_path has a value which is a string
	protein_count has a value which is an int
	stages_completed has a value which is a reference to a list where each element is a string

</pre>

=end html

=begin text

$params is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
$output is a kbase_protein_query_module.ProteinQueryAnalysisResults
ProteinQueryAnalysisResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	analysis_result_ref has a value which is a string
	summary has a value which is a string
	input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
	start_time has a value which is a float
	html_report_path has a value which is a string
	protein_count has a value which is an int
	stages_completed has a value which is a reference to a list where each element is a string


=end text

=item Description



=back

=cut

 sub run_protein_query_analysis
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function run_protein_query_analysis (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_protein_query_analysis:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'run_protein_query_analysis');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kbase_protein_query_module.run_protein_query_analysis",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'run_protein_query_analysis',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method run_protein_query_analysis",
					    status_line => $self->{client}->status_line,
					    method_name => 'run_protein_query_analysis',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "kbase_protein_query_module.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kbase_protein_query_module.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'run_protein_query_analysis',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method run_protein_query_analysis",
            status_line => $self->{client}->status_line,
            method_name => 'run_protein_query_analysis',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for kbase_protein_query_module::kbase_protein_query_moduleClient\n";
    }
    if ($sMajor == 0) {
        warn "kbase_protein_query_module::kbase_protein_query_moduleClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 ReportResults

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string


=end text

=back



=head2 CheckProteinExistenceResults

=over 4



=item Description

Check if a protein exists in the storage system using UniProt ID and create a workspace object with the result.
Input: UniProt ID (e.g., P00001, P12345)
Output: Existence status, family assignment, metadata, optional embedding


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
exists has a value which is an int
family_id has a value which is a string
metadata has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
summary has a value which is a string
protein_existence_result_ref has a value which is a string
embedding_result_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
exists has a value which is an int
family_id has a value which is a string
metadata has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
summary has a value which is a string
protein_existence_result_ref has a value which is a string
embedding_result_ref has a value which is a string


=end text

=back



=head2 ProteinExistenceResult

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
protein_id has a value which is a string
exists has a value which is an int
family_id has a value which is a string
metadata has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
embedding_ref has a value which is a string
embedding has a value which is a reference to a list where each element is a float
model_name has a value which is a string
search_timestamp has a value which is a float
summary has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
protein_id has a value which is a string
exists has a value which is an int
family_id has a value which is a string
metadata has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
embedding_ref has a value which is a string
embedding has a value which is a reference to a list where each element is a float
model_name has a value which is a string
search_timestamp has a value which is a float
summary has a value which is a string


=end text

=back



=head2 GenerateProteinEmbeddingResults

=over 4



=item Description

Generate protein embeddings from direct sequence input.
Creates embeddings using ESM-2 model for downstream analysis.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
embedding_result_ref has a value which is a string
summary has a value which is a string
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
embedding_norm has a value which is a float
sequence_length has a value which is an int
embedding_dim has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
embedding_result_ref has a value which is a string
summary has a value which is a string
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
embedding_norm has a value which is a float
sequence_length has a value which is an int
embedding_dim has a value which is an int


=end text

=back



=head2 ProteinEmbeddingResult

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
input_id has a value which is a string
input_source has a value which is a string
embedding_ref has a value which is a string
embedding has a value which is a reference to a list where each element is a float
model_name has a value which is a string
pooling_method has a value which is a string
metadata has a value which is a reference to a hash where the key is a string and the value is a string
sequence_length has a value which is an int
embedding_norm has a value which is a float
embedding_dim has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
input_id has a value which is a string
input_source has a value which is a string
embedding_ref has a value which is a string
embedding has a value which is a reference to a list where each element is a float
model_name has a value which is a string
pooling_method has a value which is a string
metadata has a value which is a reference to a hash where the key is a string and the value is a string
sequence_length has a value which is an int
embedding_norm has a value which is a float
embedding_dim has a value which is an int


=end text

=back



=head2 AssignFamilyFastResults

=over 4



=item Description

Assign a protein embedding to a family using similarity to family centroids.
Uses binary Hamming distance for fast family assignment.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
family_id has a value which is a string
confidence has a value which is a float
eigenprotein_id has a value which is a string
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
family_assignment_result_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
family_id has a value which is a string
confidence has a value which is a float
eigenprotein_id has a value which is a string
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
family_assignment_result_ref has a value which is a string


=end text

=back



=head2 ProteinFamilyAssignmentResult

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
input_id has a value which is a string
input_type has a value which is a string
embedding_ref has a value which is a string
assigned_family_id has a value which is a string
similarity_score has a value which is a float
metadata has a value which is a reference to a hash where the key is a string and the value is a string
eigenprotein_id has a value which is a string
confidence has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
input_id has a value which is a string
input_type has a value which is a string
embedding_ref has a value which is a string
assigned_family_id has a value which is a string
similarity_score has a value which is a float
metadata has a value which is a reference to a hash where the key is a string and the value is a string
eigenprotein_id has a value which is a string
confidence has a value which is a float


=end text

=back



=head2 FindTopMatchesFromEmbeddingResults

=over 4



=item Description

Find top matches for a given protein embedding within a family.
Uses FAISS IVF float index for efficient similarity search.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
matches has a value which is a reference to a list where each element is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
summary has a value which is a string
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
family_id has a value which is a string
top_n has a value which is an int
similarity_stats has a value which is a reference to a hash where the key is a string and the value is a float
similarity_search_result_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
matches has a value which is a reference to a list where each element is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
summary has a value which is a string
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
family_id has a value which is a string
top_n has a value which is an int
similarity_stats has a value which is a reference to a hash where the key is a string and the value is a float
similarity_search_result_ref has a value which is a string


=end text

=back



=head2 ProteinSimilaritySearchResult

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
input_id has a value which is a string
input_type has a value which is a string
embedding_ref has a value which is a string
family_id has a value which is a string
top_n has a value which is an int
matches has a value which is a reference to a list where each element is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
similarity_stats has a value which is a reference to a hash where the key is a string and the value is a float
metadata has a value which is a reference to a hash where the key is a string and the value is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
input_id has a value which is a string
input_type has a value which is a string
embedding_ref has a value which is a string
family_id has a value which is a string
top_n has a value which is an int
matches has a value which is a reference to a list where each element is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
similarity_stats has a value which is a reference to a hash where the key is a string and the value is a float
metadata has a value which is a reference to a hash where the key is a string and the value is a string


=end text

=back



=head2 SummarizeAndVisualizeResultsResults

=over 4



=item Description

Summarize and visualize protein network analysis results.
Generates comprehensive HTML reports with network visualization.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
output_dir has a value which is a string
summary has a value which is a string
html_report_path has a value which is a string
sequence_analysis_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
output_dir has a value which is a string
summary has a value which is a string
html_report_path has a value which is a string
sequence_analysis_ref has a value which is a string


=end text

=back



=head2 SummarizeVisualizeResult

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
input_id has a value which is a string
input_type has a value which is a string
top_matches_result_ref has a value which is a string
summary_html has a value which is a string
metadata has a value which is a reference to a hash where the key is a string and the value is a string
embedding_ref has a value which is a string
family_assignment_ref has a value which is a string
protein_existence_ref has a value which is a string
matches has a value which is a reference to a list where each element is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
family_id has a value which is a string
top_n has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
input_id has a value which is a string
input_type has a value which is a string
top_matches_result_ref has a value which is a string
summary_html has a value which is a string
metadata has a value which is a reference to a hash where the key is a string and the value is a string
embedding_ref has a value which is a string
family_assignment_ref has a value which is a string
protein_existence_ref has a value which is a string
matches has a value which is a reference to a list where each element is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
family_id has a value which is a string
top_n has a value which is an int


=end text

=back



=head2 ProteinQueryAnalysisResults

=over 4



=item Description

Unified Protein Query Analysis Pipeline

This method provides a single entry point for comprehensive protein analysis,
supporting multiple input types and configurable analysis stages.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
analysis_result_ref has a value which is a string
summary has a value which is a string
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
html_report_path has a value which is a string
protein_count has a value which is an int
stages_completed has a value which is a reference to a list where each element is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
analysis_result_ref has a value which is a string
summary has a value which is a string
input_parameters has a value which is a reference to a hash where the key is a string and the value is an UnspecifiedObject, which can hold any non-null object
start_time has a value which is a float
html_report_path has a value which is a string
protein_count has a value which is an int
stages_completed has a value which is a reference to a list where each element is a string


=end text

=back



=cut

package kbase_protein_query_module::kbase_protein_query_moduleClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
