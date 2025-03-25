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
use std::path::PathBuf;

use clap::Parser;
use log::info;
use needletail::Sequence;
use needletail::parser::SequenceRecord;
use rayon::iter::ParallelIterator;
use rayon::iter::IntoParallelRefIterator;

// Command-line interface
mod cli;
mod vcf_writer;

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
) -> Vec<(String, Vec<u8>)> {
    let mut seq_data: Vec<(String, Vec<u8>)> = Vec::new();
    let mut reader = needletail::parse_fastx_file(file).unwrap_or_else(|_| panic!("Expected valid fastX file at {}", file));
    while let Some(rec) = read_from_fastx_parser(&mut *reader) {
        let seqrec = rec.normalize(true);
        let seqname = String::from_utf8(rec.id().to_vec()).expect("UTF-8");
        seq_data.push((seqname, seqrec.to_vec()));
    }
    seq_data
}

fn read_input_list(
    input_list_file: &String,
    delimiter: u8
) -> Vec<(String, PathBuf)> {
    let fs = match std::fs::File::open(input_list_file) {
        Ok(fs) => fs,
        Err(e) => panic!("  Error in reading --input-list: {}", e),
    };

    let mut reader = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(false)
        .from_reader(fs);

    reader.records().map(|line| {
        if let Ok(record) = line {
            if record.len() > 1 {
                (record[0].to_string(), PathBuf::from(record[1].to_string()))
            } else {
                (record[0].to_string(), PathBuf::from(record[0].to_string()))
            }
        } else {
            panic!("  Error in reading --input-list: {}", input_list_file);
        }
    }).collect::<Vec<(String, PathBuf)>>()
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
            input_list,
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

            let mut sbwt_build_options = kbo::BuildOpts::default();
            sbwt_build_options.k = *kmer_size;
            sbwt_build_options.num_threads = *num_threads;
            sbwt_build_options.prefix_precalc = *prefix_precalc;
            sbwt_build_options.dedup_batches = *dedup_batches;
            sbwt_build_options.mem_gb = *mem_gb;
            sbwt_build_options.temp_dir = temp_dir.clone();

            let mut in_files = seq_files.clone();
            if let Some(list) = input_list {
                let contents = read_input_list(list, b'\t');
                let contents_iter = contents.iter().map(|(_, path)| path.to_str().unwrap().to_string());
                in_files.extend(contents_iter);
            }

            info!("Building SBWT index from {} files...", in_files.len());
            let mut seq_data: Vec<Vec<u8>> = Vec::new();
            in_files.iter().for_each(|file| {
                seq_data.append(&mut read_fastx_file(file).into_iter().map(|(_, seq)| seq).collect::<Vec<Vec<u8>>>());
            });

            let (sbwt, lcs) = kbo::build(&seq_data, sbwt_build_options);

            info!("Serializing SBWT index to {}.sbwt ...", output_prefix.as_ref().unwrap());
            info!("Serializing LCS array to {}.lcs ...", output_prefix.as_ref().unwrap());
            kbo::index::serialize_sbwt(output_prefix.as_ref().unwrap(), &sbwt, &lcs);

        },

        Some(cli::Commands::Find {
            query_files,
            input_list,
            ref_file,
            output_file,
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
            let mut sbwt_build_options = kbo::BuildOpts::default();
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
                    let ref_data = read_fastx_file(ref_file.as_ref().unwrap()).into_iter().map(|(_, seq)| seq).collect::<Vec<Vec<u8>>>();
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


            let mut in_files: Vec<(String, PathBuf)> = query_files.iter().map(|file| (file.clone(), PathBuf::from(file))).collect();
            if let Some(list) = input_list {
                let mut contents = read_input_list(list, b'\t');
                in_files.append(&mut contents);
            }

            let mut ofs = if output_file.is_some() {
                let ofs = match std::fs::File::create(output_file.as_ref().unwrap()) {
                    Ok(file) => file,
                    Err(e) => panic!("  Error in opening --output: {}", e),
                };
                Some(ofs)
            } else {
                None
            };

            info!("Querying SBWT index...");
            if let Some(ofs) = &mut ofs {
                let _ = ofs.write(b"query\tref\tq.start\tq.end\tstrand\tlength\tmismatches\tgap_bases\tgap_opens\tidentity\tcoverage\tquery.contig\tref.contig");
            } else {
                let stdout = std::io::stdout();
                let _ = stdout.lock().write(b"query\tref\tq.start\tq.end\tstrand\tlength\tmismatches\tgap_bases\tgap_opens\tidentity\tcoverage\tquery.contig\tref.contig");
            }
            in_files.iter().for_each(|(file, path)| {
                let mut run_lengths: Vec<(kbo::format::RLE, char, String, String, usize, usize)> = indexes.par_iter().map(|((sbwt, lcs), ref_contig, ref_bases)| {
                    let mut reader = needletail::parse_fastx_file(path).ok().unwrap();
                    let mut res: Vec<(kbo::format::RLE, char, String, String, usize, usize)> = Vec::new();
                    while let Some(seqrec) = read_from_fastx_parser(&mut *reader) {
                        let query_contig = std::str::from_utf8(seqrec.id()).expect("UTF-8");
                        let seq = seqrec.normalize(true);
                        // Get local alignments for forward strand
                        res.append(&mut kbo::find(&seq, sbwt, lcs, find_opts)
                                   .iter()
                                   .map(|x| (*x, '+',
                                             ref_contig.clone(), query_contig.to_string().clone(),
                                             *ref_bases, seq.len()
                                   )).collect());

                        // Add local alignments for reverse _complement
                        res.append(&mut kbo::find(&seq.reverse_complement(), sbwt, lcs, find_opts)
                                   .iter()
                                   .map(|x| (*x, '-',
                                             ref_contig.clone(), query_contig.to_string().clone(),
                                             *ref_bases, seq.len()
                                   )).collect());
                    }
                    res
                }).flatten().collect();

                // Sort by q.start
                run_lengths.sort_by_key(|(aln, _, _, _, _, _)| aln.start);

                // Print results with query and ref name added
                run_lengths.iter().filter(|x| x.0.end - x.0.start + 1 >= *min_len as usize)
                                  .for_each(|(aln, strand, ref_contig, query_contig, ref_bases, query_bases)| {
                                      let aln_len = aln.end - aln.start;
                                      let aln_start = if *strand == '+' { aln.start } else { query_bases - aln.end } + 1;
                                      let aln_end = if *strand == '+' { aln.end } else { query_bases - aln.start };
                                      let coverage = (aln.matches as f64 + aln.mismatches as f64)/(*ref_bases as f64) * 100_f64;
                                      let identity = (aln.matches as f64)/(aln_len as f64) * 100_f64;
                                      if ofs.is_some() {
                                          let _ = writeln!(&mut ofs.as_ref().unwrap(),
                                                           "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}",
                                                           file, ref_file.clone().unwrap(),
                                                           aln_start,
                                                           aln_end,
                                                           strand,
                                                           aln.end - aln.start,
                                                           aln.mismatches,
                                                           aln.gap_bases,
                                                           aln.gap_opens,
                                                           identity,
                                                           coverage,
                                                           query_contig,
                                                           ref_contig);
                                      } else {
                                          let stdout = std::io::stdout();
                                          let _ = writeln!(&mut stdout.lock(),
                                                           "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}",
                                                           file, ref_file.clone().unwrap(),
                                                           aln_start,
                                                           aln_end,
                                                           strand,
                                                           aln.end - aln.start,
                                                           aln.mismatches,
                                                           aln.gap_bases,
                                                           aln.gap_opens,
                                                           identity,
                                                           coverage,
                                                           query_contig,
                                                           ref_contig);
                                      }
                });
            });
        },
        Some(cli::Commands::Call{
            query_file,
            ref_file,
            output_file,
            max_error_prob,
            num_threads,
            kmer_size,
            prefix_precalc,
            dedup_batches,
            verbose,
        }) => {
            init_log(if *verbose { 2 } else { 1 });

            let mut sbwt_build_options = kbo::BuildOpts::default();
            // These are required for the subcommand to work correctly
            sbwt_build_options.add_revcomp = true;
            sbwt_build_options.build_select = true;
            // These can be adjusted
            sbwt_build_options.k = *kmer_size;
            sbwt_build_options.num_threads = *num_threads;
            sbwt_build_options.prefix_precalc = *prefix_precalc;
            sbwt_build_options.dedup_batches = *dedup_batches;

            let mut call_opts = kbo::CallOpts::default();
            call_opts.sbwt_build_opts = sbwt_build_options.clone();
            call_opts.max_error_prob = *max_error_prob;

            rayon::ThreadPoolBuilder::new()
                .num_threads(*num_threads)
                .thread_name(|i| format!("rayon-thread-{}", i))
                .build()
                .unwrap();

            let ref_data: Vec<(String, Vec<u8>)> = read_fastx_file(ref_file.to_str().unwrap());

            let query_data: Vec<Vec<u8>> = read_fastx_file(query_file.to_str().unwrap()).iter().map(|(_, seq)| seq.clone()).collect();
            let (sbwt_query, lcs_query) = kbo::index::build_sbwt_from_vecs(&query_data, &Some(sbwt_build_options.clone()));

            // Will map separately against each contig in the reference
            if output_file.is_some() {
                // this is dumb, surely there is a way to pass stdout or the file?
                let mut ofs = match std::fs::File::create(output_file.as_ref().unwrap()) {
                    Ok(file) => file,
                    Err(e) => panic!("  Error in opening --output: {}", e),
                };

                let vcf_header = vcf_writer::write_vcf_header(&mut ofs,
                                                              ref_file.to_str().unwrap(),
                                                              &ref_data.iter().map(|(contig_name, contig_seq)| {
                                                                  (contig_name.clone(), contig_seq.len())
                                                              }).collect::<Vec<(String, usize)>>())
                    .expect("Write header to .vcf file");

                ref_data.iter().for_each(|(ref_contig_header, ref_seq)| {
                    let calls = kbo::call(&sbwt_query, &lcs_query, ref_seq, call_opts.clone());
                    vcf_writer::write_vcf_contents(
                        &mut ofs,
                        &vcf_header,
                        &calls,
                        ref_seq,
                        ref_contig_header,
                    ).expect("Wrote .vcf record");
                });
            } else {
                let stdout = std::io::stdout();
                let vcf_header = vcf_writer::write_vcf_header(&mut stdout.lock(),
                                                              ref_file.to_str().unwrap(),
                                                              &ref_data.iter().map(|(contig_name, contig_seq)| {
                                                                  (contig_name.clone(), contig_seq.len())
                                                              }).collect::<Vec<(String, usize)>>())
                    .expect("Write header to .vcf file");

                ref_data.iter().for_each(|(ref_contig_header, ref_seq)| {
                    let calls = kbo::call(&sbwt_query, &lcs_query, ref_seq, call_opts.clone());
                    vcf_writer::write_vcf_contents(
                        &mut stdout.lock(),
                        &vcf_header,
                        &calls,
                        ref_seq,
                        ref_contig_header,
                    ).expect("Wrote .vcf record");
                });
            }
        }

        Some(cli::Commands::Map {
            query_files,
            ref_file,
            input_list,
            output_file,
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
            let mut sbwt_build_options = kbo::BuildOpts::default();
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
            map_opts.sbwt_build_opts = sbwt_build_options.clone();

            rayon::ThreadPoolBuilder::new()
                .num_threads(*num_threads)
                .thread_name(|i| format!("rayon-thread-{}", i))
                .build()
                .unwrap();

            let mut in_files: Vec<(String, PathBuf)> = query_files.iter().map(|file| (file.clone(), PathBuf::from(file))).collect();
            if let Some(list) = input_list {
                let mut contents = read_input_list(list, b'\t');
                in_files.append(&mut contents);
            }

            let ref_data: Vec<(String, Vec<u8>)> = read_fastx_file(ref_file);

            let ofs = if output_file.is_some() {
                let ofs = match std::fs::File::create(output_file.as_ref().unwrap()) {
                    Ok(file) => file,
                    Err(e) => panic!("  Error in opening --output: {}", e),
                };
                Some(ofs)
            } else {
                None
            };

            let mut first_write = true;
            in_files.iter().for_each(|(file, path)| {
                let contigs: Vec<Vec<u8>> = read_fastx_file(path.to_str().unwrap()).iter().map(|(_, seq)| seq.clone()).collect();
                let (sbwt, lcs) = kbo::index::build_sbwt_from_vecs(&contigs, &Some(sbwt_build_options.clone()));

                let mut res: Vec<u8> = Vec::new();
                ref_data.iter().for_each(|(_, ref_seq)| {
                    res.append(&mut kbo::map(ref_seq, &sbwt, &lcs, map_opts.clone()));
                });

                if ofs.is_some() {

                    if first_write {
                        let _ = writeln!(&mut ofs.as_ref().unwrap(),
                                         ">{}\n{}", ref_file, std::str::from_utf8(&ref_data.iter().flat_map(|x| x.1.clone()).collect::<Vec<u8>>()).expect("UTF-8"));
                        first_write = false;
                    }
                    let _ = writeln!(&mut ofs.as_ref().unwrap(),
                                     ">{}\n{}", file, std::str::from_utf8(&res).expect("UTF-8"));
                } else {
                    let stdout = std::io::stdout();

                    if first_write {
                        let _ = writeln!(&mut stdout.lock(),
                                         ">{}\n{}", ref_file, std::str::from_utf8(&ref_data.iter().flat_map(|x| x.1.clone()).collect::<Vec<u8>>()).expect("UTF-8"));
                        first_write = false;
                    }
                    let _ = writeln!(&mut stdout.lock(),
                                     ">{}\n{}", file, std::str::from_utf8(&res).expect("UTF-8"));
                }
            });
        },
        None => {}
    }
}
