// kbo-cli: Command-line interface to the kbo local aligner.
//
// Copyright 2024 Tommi Mäklin [tommi@maklin.fi].

// Copyrights in this project are retained by contributors. No copyright assignment
// is required to contribute to this project.

use std::io::stdout;
// Except as otherwise noted (below and/or in individual files), this
// project is licensed under the Apache License, Version 2.0
// <LICENSE-APACHE> or <http://www.apache.org/licenses/LICENSE-2.0> or
// the MIT license, <LICENSE-MIT> or <http://opensource.org/licenses/MIT>,
// at your option.
//
use std::io::Write;

use clap::Parser;
use jseqio::record::Record;
use kbo::variant_calling::call_variants;
use kbo::variant_calling::write_in_vcf_format;
use log::info;
use needletail::Sequence;
use needletail::parser::SequenceRecord;
use rayon::iter::ParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use sbwt::BitPackedKmerSortingMem;
use sbwt::SbwtIndexBuilder;

// Command-line interface
mod cli;

// Given a needletail parser, reads the next contig sequence
fn read_from_fastx_parser(
	reader: &mut dyn needletail::parser::FastxReader,
) -> Option<SequenceRecord> {
	let rec = reader.next();
	match rec {
	    Some(Ok(seqrec)) => {
			Some(seqrec)
	    },
	    None => None,
		Some(Err(_)) => todo!(),
	}
}

// Reads all sequence data from a fastX file
fn read_fastx_file(
    file: &str,
) -> Vec<Vec<u8>> {
    let mut seq_data: Vec<Vec<u8>> = Vec::new();
    let mut reader = needletail::parse_fastx_file(file).unwrap_or_else(|_| panic!("Expected valid fastX file at {}", file));
	while let Some(rec) = read_from_fastx_parser(&mut *reader) {
		let seqrec = rec.normalize(true);
		seq_data.push(seqrec.to_vec());
	}
    seq_data
}

/// Initializes the logger with verbosity given in `log_max_level`.
fn init_log(log_max_level: usize) {
    stderrlog::new()
	.module(module_path!())
	.quiet(false)
	.verbosity(log_max_level)
	.timestamp(stderrlog::Timestamp::Off)
	.init()
	.unwrap();
}

/// Use `kbo` to list the available commands or `kbo <command>` to run.
///
/// # Input format detection
/// The sequence data is read using
/// [needletail::parser::parse_fastx_file](https://docs.rs/needletail/latest/needletail/parser/fn.parse_fastx_file.html).
///
/// Input file format (fasta or fastq) is detected automatically and
/// the files may be compressed in a
/// [DEFLATE-based](https://en.wikipedia.org/wiki/Deflate) format (.gz
/// files).
///
/// [Bzip2](https://sourceware.org/bzip2/) and
/// [liblzma](https://tukaani.org/xz/) compression (.bz2 and .xz
/// files) can be enabled using the needletail features field in
/// kbo Cargo.toml if compiling from source.
///
#[allow(clippy::needless_update)]
fn main() {
    let cli = cli::Cli::parse();

    // Subcommands:
    match &cli.command {
        Some(cli::Commands::Build {
			seq_files,
			output_prefix,
            kmer_size,
			prefix_precalc,
			dedup_batches,
			num_threads,
			mem_gb,
			temp_dir,
			verbose,
        }) => {
			init_log(if *verbose { 2 } else { 1 });

            let mut sbwt_build_options = kbo::index::BuildOpts::default();
			sbwt_build_options.k = *kmer_size;
			sbwt_build_options.num_threads = *num_threads;
			sbwt_build_options.prefix_precalc = *prefix_precalc;
			sbwt_build_options.dedup_batches = *dedup_batches;
			sbwt_build_options.mem_gb = *mem_gb;
			sbwt_build_options.temp_dir = temp_dir.clone();

			info!("Building SBWT index from {} files...", seq_files.len());
			let mut seq_data: Vec<Vec<u8>> = Vec::new();
			seq_files.iter().for_each(|file| {
				seq_data.append(&mut read_fastx_file(file));
			});

			let (sbwt, lcs) = kbo::build(&seq_data, sbwt_build_options);

			info!("Serializing SBWT index to {}.sbwt ...", output_prefix.as_ref().unwrap());
			info!("Serializing LCS array to {}.lcs ...", output_prefix.as_ref().unwrap());
			kbo::index::serialize_sbwt(output_prefix.as_ref().unwrap(), &sbwt, &lcs);

		},
        Some(cli::Commands::Find {
			query_files,
			ref_file,
			index_prefix,
			detailed,
			min_len,
			max_gap_len,
			max_error_prob,
			num_threads,
            kmer_size,
			prefix_precalc,
			dedup_batches,
			mem_gb,
			temp_dir,
			verbose,
        }) => {
			init_log(if *verbose { 2 } else { 1 });
            let mut sbwt_build_options = kbo::index::BuildOpts::default();
			sbwt_build_options.k = *kmer_size;
			sbwt_build_options.num_threads = *num_threads;
			sbwt_build_options.prefix_precalc = *prefix_precalc;
			sbwt_build_options.dedup_batches = *dedup_batches;
			sbwt_build_options.mem_gb = *mem_gb;
			sbwt_build_options.temp_dir = temp_dir.clone();

			let mut find_opts = kbo::FindOpts::default();
			find_opts.max_error_prob = *max_error_prob;
			find_opts.max_gap_len = *max_gap_len as usize;

			let mut indexes: Vec<((sbwt::SbwtIndexVariant, sbwt::LcsArray), String, usize)> = Vec::new();

			if index_prefix.is_some() && !ref_file.is_some() {
				info!("Loading SBWT index...");
				let (sbwt, lcs) = kbo::index::load_sbwt(index_prefix.as_ref().unwrap());
				let n_kmers = match sbwt {
					sbwt::SbwtIndexVariant::SubsetMatrix(ref index) => {
						index.n_kmers()
					},
				};
				indexes.push(((sbwt, lcs), index_prefix.clone().unwrap(), n_kmers + *kmer_size - 1));
			} else if !index_prefix.is_some() && ref_file.is_some() {
				info!("Building SBWT from file {}...", ref_file.as_ref().unwrap());

				if !*detailed {
					let ref_data = read_fastx_file(ref_file.as_ref().unwrap());
					let n_bases = ref_data.iter().map(|x| x.len()).reduce(|a, b| a + b).unwrap();
					indexes.push((kbo::index::build_sbwt_from_vecs(&ref_data, &Some(sbwt_build_options)), ref_file.clone().unwrap(), n_bases));
				} else {
					let mut reader = needletail::parse_fastx_file(ref_file.as_ref().unwrap()).expect("valid path/file");
					while let Some(seqrec) = read_from_fastx_parser(&mut *reader) {
						let contig = seqrec.id();
						let contig_name = std::str::from_utf8(contig).expect("UTF-8");
						let seq = seqrec.normalize(true);
						indexes.push((kbo::index::build_sbwt_from_vecs(&[seq.to_vec()], &Some(sbwt_build_options.clone())), contig_name.to_string(), seq.len()));
					}
				}

			} else {
				panic!("Ambiguous reference, supply only one of `-r/--reference` and `-i/--index`");
			};

			rayon::ThreadPoolBuilder::new()
				.num_threads(*num_threads)
				.thread_name(|i| format!("rayon-thread-{}", i))
				.build()
				.unwrap();

			info!("Querying SBWT index...");
			println!("query\tref\tq.start\tq.end\tstrand\tlength\tmismatches\tgap_opens\tidentity\tcoverage\tquery.contig\tref.contig");
			let stdout = std::io::stdout();
			query_files.iter().for_each(|file| {
				let mut run_lengths: Vec<(kbo::format::RLE, char, String, String, usize)> = indexes.par_iter().map(|((sbwt, lcs), ref_contig, ref_bases)| {
					let mut reader = needletail::parse_fastx_file(file).expect("valid path/file");
					let mut res: Vec<(kbo::format::RLE, char, String, String, usize)> = Vec::new();
					while let Some(seqrec) = read_from_fastx_parser(&mut *reader) {
						let query_contig = std::str::from_utf8(seqrec.id()).expect("UTF-8");
						let seq = seqrec.normalize(true);
						// Get local alignments for forward strand
						res.append(&mut kbo::find(&seq, sbwt, lcs, find_opts)
								   .iter()
								   .map(|x| (*x, '+',
											 ref_contig.clone(), query_contig.to_string().clone(),
											 *ref_bases
								   )).collect());

						// Add local alignments for reverse _complement
						res.append(&mut kbo::find(&seq.reverse_complement(), sbwt, lcs, find_opts)
								   .iter()
								   .map(|x| (*x, '-',
											 ref_contig.clone(), query_contig.to_string().clone(),
											 *ref_bases
								   )).collect());
					}
					res
				}).flatten().collect();

				// Sort by q.start
				run_lengths.sort_by_key(|(aln, _, _, _, _)| aln.start);

				// Print results with query and ref name added
				run_lengths.iter().filter(|x| x.0.end - x.0.start + 1 >= *min_len as usize)
								  .for_each(|(aln, strand, ref_contig, query_contig, ref_bases)| {
									  let aln_len = aln.end - aln.start + 1;
									  let coverage = (aln_len as f64)/(*ref_bases as f64) * 100_f64;
									  let identity = (aln.matches as f64)/(aln_len as f64) * 100_f64;
					let _ = writeln!(&mut stdout.lock(),
									 "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}",
									 file, ref_file.clone().unwrap(),
									 aln.start,
									 aln.end,
									 strand,
									 aln.end - aln.start + 1,
									 aln.mismatches,
									 aln.gap_opens,
									 identity,
									 coverage,
									 query_contig,
									 ref_contig);
				});
			});
		},
        Some(cli::Commands::CallVariants{
			query_file,
			ref_file,
			max_error_prob,
			num_threads,
            kmer_size,
			prefix_precalc,
			dedup_batches,
			verbose,
        }) => {
			init_log(if *verbose { 2 } else { 1 });

			rayon::ThreadPoolBuilder::new()
				.num_threads(*num_threads)
				.thread_name(|i| format!("rayon-thread-{}", i))
				.build()
				.unwrap();

			let ref_seqs = jseqio::reader::DynamicFastXReader::from_file(ref_file).unwrap().into_db().unwrap();
			let ref_seq_slices: Vec<&[u8]> = ref_seqs.iter().map(|rec| rec.seq).collect();

			let query_seqs = jseqio::reader::DynamicFastXReader::from_file(query_file).unwrap().into_db().unwrap();
			if query_seqs.sequence_count() != 1 {
				eprintln!("Error: query sequence file must contain just one sequence");
				return;
			}
			let query_seq = query_seqs.get(0).seq;
			let query_name = query_seqs.get(0).name().to_vec();

			let builder = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
				.add_rev_comp(true)
				.build_lcs(true)
				.build_select_support(true)
				.k(*kmer_size)
				.n_threads(*num_threads)
				.precalc_length(*prefix_precalc)
				.algorithm(
					BitPackedKmerSortingMem::new()
					.dedup_batches(*dedup_batches)
				);

			
			let (sbwt_ref, lcs_ref) = builder.clone().run_from_slices(&ref_seq_slices);
			let (sbwt_query, lcs_query) = builder.run_from_slices(&[query_seq]);

			let calls = call_variants(
						  &sbwt_ref, 
						  lcs_ref.as_ref().unwrap(), 
						  &sbwt_query, 
						  lcs_query.as_ref().unwrap(), 
						  query_seq, 
						  *max_error_prob);

			// Todo: buffer writing?
			write_in_vcf_format(&mut std::io::stdout(), &calls, &String::from_utf8_lossy(&query_name)).unwrap();

		}

        Some(cli::Commands::Map {
			query_files,
			ref_file,
			max_error_prob,
			num_threads,
            kmer_size,
			prefix_precalc,
			dedup_batches,
			mem_gb,
			temp_dir,
			verbose,
        }) => {
			init_log(if *verbose { 2 } else { 1 });
            let mut sbwt_build_options = kbo::index::BuildOpts::default();
			// These are required for the subcommand to work correctly
			sbwt_build_options.add_revcomp = true;
			sbwt_build_options.build_select = true;
			// These can be adjusted
			sbwt_build_options.k = *kmer_size;
			sbwt_build_options.num_threads = *num_threads;
			sbwt_build_options.prefix_precalc = *prefix_precalc;
			sbwt_build_options.dedup_batches = *dedup_batches;
			sbwt_build_options.mem_gb = *mem_gb;
			sbwt_build_options.temp_dir = temp_dir.clone();

			let mut map_opts = kbo::MapOpts::default();
			map_opts.max_error_prob = *max_error_prob;

			rayon::ThreadPoolBuilder::new()
				.num_threads(*num_threads)
				.thread_name(|i| format!("rayon-thread-{}", i))
				.build()
				.unwrap();

			let ref_data = read_fastx_file(ref_file);

			let stdout = std::io::stdout();
			query_files.par_iter().for_each(|query_file| {
				let query_data = read_fastx_file(query_file);
				let (sbwt, lcs) = kbo::index::build_sbwt_from_vecs(&query_data, &Some(sbwt_build_options.clone()));

				let mut res: Vec<u8> = Vec::new();
				ref_data.iter().for_each(|ref_contig| {
					res.append(&mut kbo::map(ref_contig, &sbwt, &lcs, map_opts));
				});

				let _ = writeln!(&mut stdout.lock(),
								 ">{}\n{}", query_file, std::str::from_utf8(&res).expect("UTF-8"));
			});
		},
		None => {}
    }
}
