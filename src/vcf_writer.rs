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

use chrono::offset::Local;
use kbo::variant_calling::Variant;
use noodles_vcf::{
    header::record::value::{map::Contig, Map},
    header::record::value::Collection,
    record::{
        alternate_bases::Allele,
        genotypes::{keys::key, sample::Value, Keys},
        reference_bases::Base,
        AlternateBases, Genotypes, Position,
    },
};

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

fn split_flanking_variants(
    ref_var: &[u8],
    query_var: &[u8],
    query_pos: usize,
) -> Option<(Variant, Variant)> {
    let ref_len = ref_var.len();
    if ref_len != query_var.len() || ref_len == 1 {
        return None
    }

    let first_mismatch = ref_var[0] != query_var[0];
    let last_mismatch = ref_var[ref_len - 1] != query_var[ref_len - 1];

    let mut middle_match = true;
    for pos in 1..(ref_len - 1) {
        middle_match &= ref_var[pos] == query_var[pos];
    }

    if first_mismatch && last_mismatch && middle_match {
        Some(
            (Variant{query_chars: vec![query_var[0]], ref_chars: vec![ref_var[0]], query_pos},
             Variant{query_chars: vec![query_var[ref_len - 1]], ref_chars: vec![ref_var[ref_len - 1]], query_pos: query_pos + ref_len - 1})
        )
    } else {
        None
    }
}

fn build_vcf_record(
    variant: &Variant,
    ref_seq: &[u8],
    keys: &Keys,
    contig_name: &str,
) -> noodles_vcf::Record {
    let is_indel = variant.ref_chars.len() != variant.query_chars.len();

    let (alt_bases, ref_bases) = if is_indel {
        // Add nucleotide preceding an indel to the output
        // (.vcf does not like empty bases in REF or ALT)
        //
        // Note need to decrement query_pos by 1 later
        let alt_bases = [vec![u8_to_base(ref_seq[variant.query_pos - 1])], variant.ref_chars.iter().map(|nt| u8_to_base(*nt)).collect::<Vec<Base>>()].concat();
        let ref_bases = [vec![ref_seq[variant.query_pos - 1] as char], variant.query_chars.iter().map(|nt| *nt as char).collect::<Vec<char>>()].concat().iter().collect::<String>();
        (alt_bases, ref_bases)
    } else {
        let alt_bases = variant.ref_chars.iter().map(|nt| u8_to_base(*nt)).collect::<Vec<Base>>();
        let ref_bases = variant.query_chars.iter().map(|nt| *nt as char).collect::<String>();
        (alt_bases, ref_bases)
    };
    let alt_allele = vec![Allele::Bases(alt_bases)];
    let genotypes = Genotypes::new(keys.clone(), vec![vec![Some(Value::String("1".to_string()))]]);

    let mut record_builder = noodles_vcf::Record::builder()
        .set_chromosome(
            contig_name
                .parse()
                .expect("Invalid chromosome name"),
        )
        .set_position(Position::from(variant.query_pos + 1 - is_indel as usize * (1 + (variant.ref_chars.len() as isize - variant.query_chars.len() as isize).unsigned_abs())))
        .set_reference_bases(ref_bases.parse().expect("Reference bases"))
        .set_alternate_bases(AlternateBases::from(alt_allele))
        .set_genotypes(genotypes);

    if variant.ref_chars.len() != 1 || variant.query_chars.len() != 1 {
        record_builder = record_builder.set_info("INDEL".parse().expect("Parsed string"));
    }

    record_builder.build().expect("Build .vcf record")
}

/// Write the header of a .vcf file
pub fn write_vcf_header<W: Write>(f: &mut W,
                              ref_name: &str,
                              contig_info: &[(String, usize)],
) -> Result<noodles_vcf::Header, std::io::Error> {
    let mut writer = noodles_vcf::Writer::new(f);
    let mut header_builder = noodles_vcf::Header::builder();
    for (contig_header, length) in contig_info.iter() {
        let record = Map::<Contig>::builder()
            .set_length(*length)
            .build();

        let mut header_contents = contig_header.split_whitespace();
        let contig_name = header_contents.next().expect("Contig name");
        header_builder = header_builder.add_contig(
            contig_name.parse().expect("Query contig name in header"),
            record.expect("Record of type noodles_vcf::header::record::value::map::Contig"),
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
pub fn write_vcf_contents<W: Write>(f: &mut W,
                                    header: &noodles_vcf::Header,
                                    variants: &[Variant],
                                    ref_seq: &[u8],
                                    contig_header: &str,
) -> Result<(), std::io::Error> {
    let mut writer = noodles_vcf::Writer::new(f);

    // Write each record (column)
    let keys = Keys::try_from(vec![key::GENOTYPE]).unwrap();

    let mut header_contents = contig_header.split_whitespace();
    let contig_name = header_contents.next().expect("Contig name");

    variants.iter().for_each(|variant| {
        let flanking = split_flanking_variants(&variant.ref_chars, &variant.query_chars, variant.query_pos);
        if flanking.is_some() {
            let (var1, var2) = flanking.unwrap();
            let record1 = build_vcf_record(&var1, ref_seq, &keys, contig_name);
            let record2 = build_vcf_record(&var2, ref_seq, &keys, contig_name);
            writer.write_record(header, &record1).expect("Wrote .vcf record");
            writer.write_record(header, &record2).expect("Wrote .vcf record");
        } else {
            let record = build_vcf_record(variant, ref_seq, &keys, contig_name);
            writer.write_record(header, &record).expect("Wrote .vcf record");
        }
    });
    Ok(())
}
