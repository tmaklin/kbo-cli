// kbo-cli:Command-line interface to the kbo local aligner.
//
// Copyright 2024 Tommi Mäklin [tommi@maklin.fi].

// Copyrights in this project are retained by contributors. No copyright assignment
// is required to contribute to this project.

use std::path::PathBuf;

// Except as otherwise noted (below and/or in individual files), this
// project is licensed under the Apache License, Version 2.0
// <LICENSE-APACHE> or <http://www.apache.org/licenses/LICENSE-2.0> or
// the MIT license, <LICENSE-MIT> or <http://opensource.org/licenses/MIT>,
// at your option.
//
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(version)]
#[command(propagate_version = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
    // Build SBWT index
    Build {
        // Input fasta or fastq sequence file(s)
        #[arg(group = "input", required = true, help = "Sequence data file(s).")]
        seq_files: Vec<String>,

        #[arg(short = 'l', long = "input-list", group = "input", required = true, help_heading = "Input", help = "File with paths or tab separated name and path on each line.")]
        input_list: Option<String>,

        // Outputs
        #[arg(short = 'o', long = "output-prefix", required = true, help_heading = "Output", help = "Prefix for output files <prefix>.sbwt and <prefix>.lcs.")]
        output_prefix: Option<String>,

        // Build parameters
        // // k-mer size
        #[arg(short = 'k', default_value_t = 31, help_heading = "Build options", help = "k-mer size, larger values are slower and use more space.")]
        kmer_size: usize,
        // // prefix precalc
        #[arg(short = 'p', long = "prefix-precalc", default_value_t = 8, help_heading = "Build options", help = "Length of precalculated prefixes included in the index.")]
        prefix_precalc: usize,
        // // deduplicate k-mer batches
        #[arg(short = 'd', long = "dedup-batches", default_value_t = false, help_heading = "Build options", help = "Deduplicate k-mer batches to save some memory.")]
        dedup_batches: bool,

        // Resources
        // // Threads
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        num_threads: usize,
        // // Memory in GB
        #[arg(short = 'm', long = "mem-gb", default_value_t = 4, help_heading = "Build options", help = "Memory available when building on temp disk space (in gigabytes).")]
        mem_gb: usize,
        // // Temporary directory
        #[arg(long = "temp-dir", required = false, help_heading = "Build options", help = "Build on temporary disk space at this path instead of in-memory.")]
        temp_dir: Option<String>,

        // Verbosity
        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,
    },


    // Call variants in query relative to a reference
    Call{
        // Inputs
        // // Input fasta or fastq query file(s)
        #[arg(group = "input", required = true, help = "Query file with sequence data.")]
        query_file: PathBuf,
        // // Reference fasta or fastq file
        #[arg(long = "reference", short = 'r', required = true, help_heading = "Input", help = "Reference sequence to call variants in.")]
        ref_file: PathBuf,

        // Outputs
        #[arg(short = 'o', long = "output", required = false, help_heading = "Output", help = "Write output to a file instead of printing.")]
        output_file: Option<String>,

        // Upper bound for random match probability
        #[arg(long = "max-error-prob", default_value_t = 0.00000001, help_heading = "Algorithm", help = "Tolerance for errors in k-mer matching.")]
        max_error_prob: f64,

        // Resources
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        num_threads: usize,

        // SBWT build parameters
        // k-mer size
        #[arg(short = 'k', default_value_t = 51, help_heading = "Build options", help = "k-mer size, larger values are slower and use more space.")]
        kmer_size: usize,
        // prefix precalc
        #[arg(short = 'p', long = "prefix-precalc", default_value_t = 8, help_heading = "Build options", help = "Length of precalculated prefixes included in the index.")]
        prefix_precalc: usize,
        // deduplicate k-mer batches
        #[arg(short = 'd', long = "dedup-batches", default_value_t = false, help_heading = "Build options", help = "Deduplicate k-mer batches to save some memory.")]
        dedup_batches: bool,

        // Verbosity
        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,
    },

    // Find indexed k-mers in a query
    Find {
        // Inputs
        // // Input fasta or fastq query file(s)
        #[arg(group = "input", required = true, help = "Query file(s) with sequence data.")]
        query_files: Vec<String>,
        // // Input list
        #[arg(short = 'l', long = "input-list", group = "input", required = true, help_heading = "Input", help = "File with paths or tab separated name and path on each line.")]
        input_list: Option<String>,

        // Output options
        // // Output file
        #[arg(short = 'o', long = "output", required = false, help_heading = "Output", help = "Write output to a file instead of printing.")]
        output_file: Option<String>,

        // Reference
        // // Sequence file
        #[arg(short = 'r', long = "reference", group = "reference", help_heading = "Input", help = "File with target sequence data (excludes -i).")]
        ref_file: Option<String>,
        // // ... or a prebuilt index
        #[arg(short = 'i', long = "index", group = "reference", help_heading = "Input", help = "Prefix for prebuilt <prefix>.sbwt and <prefix>.lcs (excludes -r).")]
        index_prefix: Option<String>,
        // // Concatenate contigs in reference
        #[arg(long = "detailed", help_heading = "Input", default_value_t = false)]
        detailed: bool,

        // Parameters
        // // Minimum length to report an alignment
        #[arg(long = "min-len", default_value_t = 100, help_heading = "Algorithm", help = "Minimum alignment length to report.")]
        min_len: u64,
        #[arg(long = "max-gap-len", default_value_t = 0, help_heading = "Algorithm", help = "Allow gaps of this length in the alignment.")]
        max_gap_len: u64,

        // // Upper bound for random match probability
        #[arg(long = "max-error-prob", default_value_t = 0.0000001, help_heading = "Algorithm", help = "Tolerance for errors in k-mer matching.")]
        max_error_prob: f64,

        // Resources
        // // Threads
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        num_threads: usize,

        // Build parameters
        // // k-mer size
        #[arg(short = 'k', default_value_t = 31, help_heading = "Build options", help = "k-mer size, larger values are slower and use more space.")]
        kmer_size: usize,
        // // prefix precalc
        #[arg(short = 'p', long = "prefix-precalc", default_value_t = 8, help_heading = "Build options", help = "Length of precalculated prefixes included in the index.")]
        prefix_precalc: usize,
        // // deduplicate k-mer batches
        #[arg(short = 'd', long = "dedup-batches", default_value_t = false, help_heading = "Build options", help = "Deduplicate k-mer batches to save some memory.")]
        dedup_batches: bool,
        // // Memory in GB
        #[arg(short = 'm', long = "mem-gb", default_value_t = 4, help_heading = "Build options", help = "Memory available when building on temp disk space (in gigabytes).")]
        mem_gb: usize,
        // // Temporary directory
        #[arg(long = "temp-dir", required = false, help_heading = "Build options", help = "Build on temporary disk space at this path instead of in-memory.")]
        temp_dir: Option<String>,


        // Verbosity
        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,
    },

    // Map a query or queries to a reference and return the alignment
    Map {
        // Input fasta or fastq query file(s)
        #[arg(group = "input", required = true, help = "Sequence data file(s).")]
        query_files: Vec<String>,

        // Input options
        // // Input list
        #[arg(short = 'l', long = "input-list", group = "input", required = true, help_heading = "Input", help = "File with paths or tab separated name and path on each line.")]
        input_list: Option<String>,
        // // Reference file
        #[arg(short = 'r', long = "reference", required = true, help_heading = "Input", help = "Reference sequence to map against.")]
        ref_file: String,

        // Output options
        // // Output file
        #[arg(short = 'o', long = "output", required = false, help_heading = "Output", help = "Write output to a file instead of printing.")]
        output_file: Option<String>,
        #[arg(long = "raw", default_value_t = false, help_heading = "Output", help = "Output the internal representation instead of an alignment.")]
        skip_formatting: bool,

        // Parameters
        // // Upper bound for random match probability
        #[arg(long = "max-error-prob", default_value_t = 0.0000001, help_heading = "Algorithm", help = "Tolerance for errors in k-mer matching.")]
        max_error_prob: f64,
        // // Skip gap filling
        #[arg(long = "no-gap-filling", default_value_t = false, help_heading = "Algorithm", help = "Skip running the gap filling algorithm.")]
        skip_gap_filling: bool,
        // // Skip variant calling
        #[arg(long = "no-variant-calling", default_value_t = false, help_heading = "Algorithm", help = "Skip using variant calling to improve the alignment.")]
        // // Don't format the results, return the raw translation
        skip_variant_calling: bool,

        // Resources
        // // Threads
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        num_threads: usize,

        // Build parameters
        // // k-mer size
        #[arg(short = 'k', default_value_t = 31, help_heading = "Build options", help = "k-mer size, larger values are slower and use more space.")]
        kmer_size: usize,
        // // prefix precalc
        #[arg(short = 'p', long = "prefix-precalc", default_value_t = 8, help_heading = "Build options", help = "Length of precalculated prefixes included in the index.")]
        prefix_precalc: usize,
        // // deduplicate k-mer batches
        #[arg(short = 'd', long = "dedup-batches", default_value_t = false, help_heading = "Build options", help = "Deduplicate k-mer batches to save some memory.")]
        dedup_batches: bool,
        // // Memory in GB
        #[arg(short = 'm', long = "mem-gb", default_value_t = 4, help_heading = "Build options", help = "Memory available when building on temp disk space (in gigabytes).")]
        mem_gb: usize,
        // // Temporary directory
        #[arg(long = "temp-dir", required = false, help_heading = "Build options", help = "Build on temporary disk space at this path instead of in-memory.")]
        temp_dir: Option<String>,

        // Verbosity
        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,
    },
}
