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
                       ref_seq: &[u8],
                       mapped_seq: &[u8],
                       contig_header: &str
) -> Result<(), std::io::Error> {
    let mut writer = noodles_vcf::Writer::new(f);

    // Write each record (column)
    let keys = Keys::try_from(vec![key::GENOTYPE]).unwrap();

    for (mapped_pos, ref_base) in ref_seq.iter().enumerate() {
        let alt_base: Base;
        let mut variant = false;

        let mapped_base = mapped_seq[mapped_pos];

        let (genotype, alt_base) = if mapped_base == *ref_base {
            (String::from("0"), u8_to_base(mapped_base))
        } else if mapped_base == b'-' {
            // Only output changes that can be resolved
            variant = false;
            (String::from("."), u8_to_base(b'-'))
        } else {
            variant = true;
            alt_base = u8_to_base(mapped_base);
            (mapped_base.to_string(), alt_base)
        };

        if variant {
            let ref_allele = u8_to_base(*ref_base);
            let genotypes = Genotypes::new(keys.clone(), vec![vec![Some(Value::String(genotype))]]);
            let alt_allele = vec![Allele::Bases(vec![alt_base])];

            let mut header_contents = contig_header.split_whitespace();
            let contig_name = header_contents.next().expect("Contig name");

            let record = noodles_vcf::Record::builder()
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
