// kbo-cli: Command-line interface to the kbo local aligner.
//
// Copyright 2024 Tommi Mäklin [tommi@maklin.fi].

// Copyrights in this project are retained by contributors. No copyright assignment
// is required to contribute to this project.

// Except as otherwise noted (below and/or in individual files), this
// project is licensed under the Apache License, Version 2.0
// <LICENSE-APACHE> or <http://www.apache.org/licenses/LICENSE-2.0> or
// the MIT license, <LICENSE-MIT> or <http://opensource.org/licenses/MIT>,
// at your option.
//
use std::io::Write;

use clap::Parser;
use chrono::offset::Local;
use log::info;
use needletail::Sequence;
use needletail::parser::SequenceRecord;
use rayon::iter::ParallelIterator;
use rayon::iter::IntoParallelRefIterator;

// TODO clean up
use noodles_vcf::{
    self as vcf,
    header::record::value::{map::Contig, Map},
    header::record::value::Collection,
    record::{
        alternate_bases::Allele,
        genotypes::{keys::key, sample::Value, Keys},
        reference_bases::Base,
        AlternateBases, Genotypes, Position,
    },
};

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

/// [`u8`] representation used elsewhere to [`noodles_vcf::record::reference_bases::Base`]
#[inline]
fn u8_to_base(ref_base: u8) -> Base {
    match ref_base {
        b'A' => Base::A,
        b'C' => Base::C,
        b'G' => Base::G,
        b'T' => Base::T,
        _ => Base::N,
    }
}

/// Write the header of a .vcf file
fn write_vcf_header<W: Write>(f: &mut W,
                              ref_name: &str,
                              contig_info: &[(&String, &usize)],
) -> Result<vcf::Header, std::io::Error> {
    let mut writer = vcf::Writer::new(f);
    let mut header_builder = vcf::Header::builder();
    for (name, length) in contig_info.iter() {
        let record = Map::<Contig>::builder()
            .set_length(**length)
            .build();

        header_builder = header_builder.add_contig(
            name.parse().expect("Could not add contig to header"),
            record.expect("Valid contig"),
        );

    };
    header_builder = header_builder.add_sample_name("unknown");
    let mut header = header_builder.build();

    let current_date = Local::now().format("%Y%m%d").to_string();
    let vcf_date = Collection::Unstructured(vec![current_date]);
    header.other_records_mut().insert("fileDate".parse().expect("Valid string"), vcf_date.clone());

    let vcf_source = Collection::Unstructured(vec![format!("kbo-cli v{}", env!("CARGO_PKG_VERSION"))]);
    header.other_records_mut().insert("source".parse().expect("Valid string"), vcf_source.clone());

    let vcf_reference = Collection::Unstructured(vec![ref_name.to_string()]);
    header.other_records_mut().insert("reference".parse().expect("Valid string"), vcf_reference.clone());

    let vcf_phasing = Collection::Unstructured(vec!["none".to_string()]);
    header.other_records_mut().insert("phasing".parse().expect("Valid string"), vcf_phasing.clone());

    writer.write_header(&header)?;

    Ok(header)
}

/// Write the contents of a .vcf file
fn write_vcf<W: Write>(f: &mut W,
                       header: &vcf::Header,
                       ref_seq: &[&u8],
                       mapped_seq: &[u8],
                       contig_name: &str
) -> Result<(), std::io::Error> {
    let mut writer = vcf::Writer::new(f);

    // Write each record (column)
    let keys = Keys::try_from(vec![key::GENOTYPE]).unwrap();

    for (mapped_pos, ref_base) in ref_seq.iter().enumerate() {
        let alt_base: Base;
        let mut variant = false;

        let mapped_base = mapped_seq[mapped_pos];

        let (genotype, alt_base) = if mapped_base == **ref_base {
            (String::from("0"), u8_to_base(mapped_base))
        } else if mapped_base == b'-' {
            variant = true;
            (String::from("."), u8_to_base(b'-'))
        } else {
            variant = true;
            alt_base = u8_to_base(mapped_base);
            (mapped_base.to_string(), alt_base)
        };

        if variant {
            let ref_allele = u8_to_base(**ref_base);
            let genotypes = Genotypes::new(keys.clone(), vec![vec![Some(Value::String(genotype))]]);
            let alt_allele = vec![Allele::Bases(vec![alt_base])];
            let record = vcf::Record::builder()
                .set_chromosome(
                    contig_name
                        .parse()
                        .expect("Invalid chromosome name"),
                )
                .set_position(Position::from(mapped_pos + 1))
                .add_reference_base(ref_allele)
                .set_alternate_bases(AlternateBases::from(alt_allele))
                .set_genotypes(genotypes)
                .build()
                .expect("Could not construct record");
            writer.write_record(header, &record)?;
        }
    }
    Ok(())
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

        Some(cli::Commands::Map {
			query_files,
			ref_file,
            detailed,
            out_format,
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

            let mut indexes: Vec<((sbwt::SbwtIndexVariant, sbwt::LcsArray), String, usize)> = Vec::new();

            query_files.iter().for_each(|query_file| {
                if *detailed {
                    let mut reader = needletail::parse_fastx_file(query_file).expect("valid path/file");
                    while let Some(contig) = read_from_fastx_parser(&mut *reader) {
                        let name = std::str::from_utf8(contig.id()).expect("UTF-8");
                        let seq = contig.normalize(true);
                        indexes.push((kbo::index::build_sbwt_from_vecs(&[seq.to_vec()], &Some(sbwt_build_options.clone())), name.to_string(), seq.len()));
                    }
                } else {
                    let contigs = read_fastx_file(query_file);
                    let n_bases = ref_data.iter().map(|x| x.len()).reduce(|a, b| a + b).unwrap();
                    let name = query_file.clone();
                    indexes.push((kbo::index::build_sbwt_from_vecs(&contigs, &Some(sbwt_build_options.clone())), name, n_bases));
                }

                if out_format == "aln" {
                    ref_data.iter().for_each(|ref_contig| {
                        indexes.iter().for_each(|((sbwt, lcs), contig_name, _)| {
                            let res = kbo::map(ref_contig, &sbwt, &lcs, map_opts);
                            let _ = writeln!(&mut stdout.lock(),
                                             ">{}\n{}", contig_name, std::str::from_utf8(&res).expect("UTF-8"));
                        });
                    });
                } else if out_format == "vcf" {
                    let vcf_header = write_vcf_header(&mut stdout.lock(), ref_file,
                                     &indexes.iter().map(|((_, _), contig_name, contig_len)| {
                                         (contig_name, contig_len)
                                     }).collect::<Vec<(&String, &usize)>>())
                        .expect("I/O error");

                        ref_data.iter().for_each(|ref_contig| {
                            indexes.iter().for_each(|((sbwt, lcs), contig_name, _)| {
                                let res = kbo::map(ref_contig, &sbwt, &lcs, map_opts);
                                let _ = write_vcf(&mut stdout.lock(), &vcf_header, &ref_data.iter().flatten().collect::<Vec<&u8>>(), &res, &contig_name);
                            });
                        });
                };
            });
        },
        None => {}
    }
}
